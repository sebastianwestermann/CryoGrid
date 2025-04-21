%========================================================================
% CryoGrid TILE class TILE_ML
% TILE class designed for machine learning 

% S. Westermann, April 2025
%========================================================================

classdef TILE_ML_ENSEMBLE < matlab.mixin.Copyable
    
    properties
        
        PARA
        RUN_INFO
        FORCING
        CONST
        OUT        
        STORE
        TEMP
        STATVAR
        ML
        ML_STORE
        t
        FEATURE_CLASS
    end
    
    
    methods
        
       

        function tile = provide_PARA(tile)

            tile.PARA.number_of_MLs = [];
            tile.PARA.variables = [];
            
            tile.PARA.ml_class = [];
            tile.PARA.ml_class_index = [];

            tile.PARA.training_class = [];
            tile.PARA.training_class_index = [];

            tile.PARA.target_data_class = [];
            tile.PARA.target_data_class_index = [];

            tile.PARA.feature_data_class = [];
            tile.PARA.feature_data_class_index = [];

            tile.PARA.out_class = [];
            tile.PARA.out_class_index = [];

            tile.PARA.read_from_store = 0;
            tile.PARA.strip4store = 0;
            
            tile.PARA.worker_number = 1;
        end

        function tile = provide_CONST(tile)

        end
        
        function tile = provide_STATVAR(tile)

        end

        function tile = finalize_init(tile)

            tile.PARA.run_name =  tile.RUN_INFO.PPROVIDER.PARA.run_name;
            tile.PARA.result_path =  tile.RUN_INFO.PPROVIDER.PARA.result_path;
            
            disp('getting training data')

            rng(1234+tile.PARA.worker_number)
            target_data_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.target_data_class){tile.PARA.target_data_class_index,1}); %read out
            target_data_class.PARA.variables = tile.PARA.variables;
            target_data_class = finalize_init(target_data_class, tile);

            feature_data_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.feature_data_class){tile.PARA.feature_data_class_index,1}); %read terrain and forcing
            feature_data_class = finalize_init(feature_data_class, tile);
            tile.STATVAR.in = generate_feature_data_training(feature_data_class, tile);
            tile.FEATURE_CLASS = feature_data_class;
            target_ML = generate_target_data(target_data_class, tile);

            if tile.PARA.read_from_store
                tile.ML_STORE = tile.RUN_INFO.TILE.ML_STORE;
                tile.ML_STORE.in = [tile.ML_STORE.in; tile.STATVAR.in];
                tile.ML_STORE.out = [tile.ML_STORE.out; target_ML];
            else
                tile.ML_STORE.in = tile.STATVAR.in;
                tile.ML_STORE.out = target_ML;
            end

            %divide in and out by std

            disp('training neural net')

            tile.STATVAR.out_mean = mean(target_ML,1);
            tile.STATVAR.out_std = std(target_ML,[],1);
            tile.STATVAR.in_mean = mean(tile.STATVAR.in,1);
            tile.STATVAR.in_std = std(tile.STATVAR.in,[],1);     

            target_ML = (target_ML - tile.STATVAR.out_mean) ./ tile.STATVAR.out_std;
            tile.STATVAR.in = (tile.STATVAR.in - tile.STATVAR.in_mean) ./ tile.STATVAR.in_std;
            
            for i=1:size(tile.PARA.variables,1)
                tile.TEMP.var_ID = i;
                tile.STATVAR.out = target_ML(:,i);
                for j=1:tile.PARA.number_of_MLs
                    tile.ML = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.ml_class){tile.PARA.ml_class_index,1});
                    tile.ML = finalize_init(tile.ML, tile);

                    training_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.training_class){tile.PARA.training_class_index,1});
                    training_class = finalize_init(training_class, tile);
                    tile.ML = train_ML(training_class, tile);
                    tile.ML_STORE.ML{i,j} = tile.ML;
                end
            end
            tile.STATVAR.out = target_ML;

            tile = assign_timestamp_prediction(feature_data_class, tile);

            tile.OUT = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.out_class){tile.PARA.out_class_index,1});
            tile.OUT = finalize_init(tile.OUT, tile);
        end

        function tile = run_model(tile)
            disp('running neural net forward')
            for var_ID = 1:size(tile.PARA.variables,1)
                for i = 1:size(tile.STATVAR.timestamp_prediction,1)
                    tile.t = tile.STATVAR.timestamp_prediction(i,1);
                    in = generate_feature_data_prediction_single_timestamp(tile.FEATURE_CLASS, tile);
                    in = (in - tile.STATVAR.in_mean) ./ tile.STATVAR.in_std;
                    tile.STATVAR.(tile.PARA.variables{var_ID,1}) = [];
                    for j=1:tile.PARA.number_of_MLs
                        tile.ML = tile.ML_STORE.ML{var_ID,j};
                        [out_j, ~] = progapagate_ML(tile.ML, in);
                        out_j = mean(out_j,2) .* tile.STATVAR.out_std(1,var_ID) + tile.STATVAR.out_mean(1,var_ID);
                        tile.STATVAR.(tile.PARA.variables{var_ID,1}) = [tile.STATVAR.(tile.PARA.variables{var_ID,1}) out_j];
                    end
                    tile = store_OUT_tile(tile, var_ID);
                end
            end

            tile = write_OUT_tile(tile); %write to SPATIAL and excahnge info between tiles if parallel
            
            if tile.PARA.strip4store
                tile.OUT = [];
                tile.ML = [];
                tile.FEATURE_CLASS = [];
                tile.PARA = [];
                tile.RUN_INFO = [];
                tile.FORCING = [];
                tile.CONST = [];
                tile.STORE = [];
                tile.TEMP = [];
                tile.STATVAR = [];
                tile.t = [];
            end
        end
        

        function tile = store_OUT_tile(tile, var_ID)
            tile.OUT = store_OUT(tile.OUT, tile, var_ID);
        end        

        function tile = write_OUT_tile(tile)
            tile.OUT = write_OUT(tile.OUT, tile);
        end  
        
        
%         function tile = run_model(tile)
% 
% 
%             TOP = tile.TOP;
%             BOTTOM = tile.BOTTOM;
%             TOP.LATERAL = tile.LATERAL;
% 
%             %=========================================================================
%             %TIME INTEGRATION
%             %=========================================================================
%             while tile.t < tile.FORCING.PARA.end_time
% 
%                 %interpolate focing data to time t
%                 tile = interpolate_forcing_tile(tile);
% 
%                 %upper boundar condition (uppermost class only)
%                 TOP.NEXT = get_boundary_condition_u(TOP.NEXT, tile);
% 
%                 %set fluxes between classes in the stratigrapht
%                 CURRENT = TOP.NEXT;
%                 while ~isequal(CURRENT.NEXT, BOTTOM)
%                     get_boundary_condition_m(CURRENT.IA_NEXT, tile); %call interaction class function
%                     CURRENT = CURRENT.NEXT;
%                 end
% 
%                 %lower boundary condition (lowermost class)
%                 CURRENT = get_boundary_condition_l(CURRENT,  tile);  %At this point, CURRENT is equal to BOTTOM_CLASS
% 
%                 %calculate spatial derivatives
%                 CURRENT = TOP.NEXT;
%                 while ~isequal(CURRENT, BOTTOM)
%                     CURRENT = get_derivatives_prognostic(CURRENT, tile);
%                     CURRENT = CURRENT.NEXT;
%                 end
% 
%                 %calculate timestep [second]
%                 CURRENT = TOP.NEXT;
%                 tile.timestep = 1e8;
%                 while ~isequal(CURRENT, BOTTOM)
%                     tile.timestep = min(tile.timestep, get_timestep(CURRENT, tile));
%                     CURRENT = CURRENT.NEXT;
%                 end   
%                 tile.next_break_time = min(tile.LATERAL.IA_TIME, tile.OUT.OUTPUT_TIME);
%                 tile.timestep = min(tile.timestep, (tile.next_break_time - tile.t).*tile.CONST.day_sec);
% 
% 
%                 %prognostic step - integrate prognostic variables in time
%                 CURRENT = TOP.NEXT;
%                 while ~isequal(CURRENT, BOTTOM)
%                     CURRENT = advance_prognostic(CURRENT, tile);
%                     CURRENT = CURRENT.NEXT;
%                 end
% 
%                 %diagnostic step - compute diagnostic variables
%                 TOP.NEXT = compute_diagnostic_first_cell(TOP.NEXT, tile); %calculate Lstar, only uppermost class
%                 CURRENT = BOTTOM.PREVIOUS;
%                 while ~isequal(CURRENT, TOP)
%                     CURRENT = compute_diagnostic(CURRENT, tile);
%                     CURRENT = CURRENT.PREVIOUS;
%                 end
% 
% 
%                 %triggers
%                 CURRENT = TOP.NEXT;
%                 while ~isequal(CURRENT, BOTTOM)
%                     CURRENT = check_trigger(CURRENT, tile);
%                     CURRENT = CURRENT.NEXT;
%                 end
% 
%                 %lateral interactions
%                 tile = interact_lateral(tile);
% 
%                 %set TOP_CLASS and BOTTOM_CLASS for convenient access
%                 tile.TOP_CLASS = TOP.NEXT;
%                 tile.BOTTOM_CLASS = BOTTOM.PREVIOUS;
% 
%                 %update time variable t
%                 tile.t = tile.t + tile.timestep./tile.CONST.day_sec;
% 
%                 %model
%                 tile = store_OUT_tile(tile);
%             end
% 
%             % extra carriage return needed for OUT classes that use
%             % fprintf to print and overwrite step date output to 
%             % console window.
%             fprintf('\n\n')
% 
%         end
% 
% 
%         %---BUILDER functions--------------
% 
%         function check_if_PARA_assigned(tile)
%             if size(tile.PARA.builder,1) == 0 && size(tile.PARA.builder,2) == 0
%                 disp(['PARA builder in class ' class(tile) ' not assigned'])
%             end
%             if strcmp(tile.PARA.builder, 'new_init') 
%                 parameters = {'domain_depth'; 'forcing_class'; 'forcing_class_index'; 'grid_class'; 'grid_class_index'; 'out_class'; ...
%                             'out_class_index'; 'strat_classes_class'; 'strat_classes_class_index'; 'strat_statvar_class'; 'strat_statvar_class_index'; 'lateral_class'; ...
%                             'lateral_class_index'; 'lateral_IA_classes'; 'lateral_IA_classes_index'};
%             elseif strcmp(tile.PARA.builder, 'new_init_steady_state')
%                 parameters = {'domain_depth'; 'forcing_class'; 'forcing_class_index'; 'grid_class'; 'grid_class_index'; 'out_class'; ...
%                     'out_class_index'; 'strat_classes_class'; 'strat_classes_class_index'; 'strat_statvar_class'; 'strat_statvar_class_index'; 'lateral_class'; ...
%                     'lateral_class_index'; 'lateral_IA_classes'; 'lateral_IA_classes_index'; 'init_steady_state_class'; 'init_steady_state_class_index';...
%                     'T_first_cell'; 'start_depth_steady_state'};
%             elseif strcmp(tile.PARA.builder, 'update_forcing_out')
%                 parameters = { 'forcing_class'; 'forcing_class_index';  'out_class'; 'out_class_index'};
%             elseif strcmp(tile.PARA.builder, 'restart_OUT_last_timestep')
%                 parameters = {'restart_file_path'; 'restart_file_name'};
%             else
%                 parameters = fieldnames(tile.PARA);
%             end
%             for i=1:size(parameters,1)
%                 if size(tile.PARA.(parameters{i,1}),1) == 0 && size(tile.PARA.(parameters{i,1}),2) == 0
%                     disp(['Warning: PARA ' parameters{i,1} ' in class ' class(tile) ' not assigned'])
%                 end
%             end
%         end
% 
%         function tile = build_tile_new_init(tile)
% 
%             tile.PARA.run_name =  tile.RUN_INFO.PPROVIDER.PARA.run_name;
%             tile.PARA.result_path =  tile.RUN_INFO.PPROVIDER.PARA.result_path;
% 
%             %1. forcing
%             %tile.FORCING = copy(tile.RUN_INFO.PPROVIDER.FUNCTIONAL_CLASSES.FORCING{tile.PARA.forcing_index,1});
%             tile.FORCING = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.forcing_class){tile.PARA.forcing_class_index,1});
%             tile.FORCING = finalize_init(tile.FORCING, tile);
% 
%             %2. grid
%             %tile.GRID = copy(tile.RUN_INFO.PPROVIDER.FUNCTIONAL_CLASSES.GRID{tile.PARA.grid_index,1});
%             tile.GRID = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.grid_class){tile.PARA.grid_class_index,1});
%             tile.GRID = finalize_init(tile.GRID, tile);
% 
%             %3. map statvar to grid, using STRATIGRAPHY_STATVAR classes 
%             for i=1:size(tile.PARA.strat_statvar_class,1)
%                 strat_statvar_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.strat_statvar_class{i,1}){tile.PARA.strat_statvar_class_index(i,1),1});
%                 strat_statvar_class = finalize_init(strat_statvar_class, tile);
%             end
% 
%             %4. build stratigraphy
%             strat_classes_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.strat_classes_class){tile.PARA.strat_classes_class_index,1});
%             strat_classes_class = finalize_init(strat_classes_class, tile);
%             class_list = strat_classes_class.PARA.classes.class_name;
%             class_index = strat_classes_class.PARA.classes.class_index;
%             tile.TOP = Top();
%             CURRENT = tile.TOP;
%             for i=1:size(class_list,1)
%                 CURRENT.NEXT = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(class_list{i,1}){class_index(i,1)});
%                 CURRENT.NEXT.PREVIOUS = CURRENT;
%                 CURRENT = CURRENT.NEXT;
%             end
%             tile.BOTTOM = Bottom();
%             CURRENT.NEXT = tile.BOTTOM;
%             tile.BOTTOM.PREVIOUS = CURRENT;
% 
%             tile.TOP_CLASS = tile.TOP.NEXT;
%             tile.BOTTOM_CLASS = tile.BOTTOM.PREVIOUS;
% 
%             %5. assign STATVAR using STRATGRAPHY_STATVAR classes
%             class_depths = strat_classes_class.PARA.classes.depth;
%             class_depths = [class_depths; tile.GRID.STATVAR.GRID(end,1)];
% 
%             CURRENT = tile.TOP_CLASS;
%             for i=1:size(class_list,1)
%                 variables = fieldnames(CURRENT.STATVAR);
%                 range = (tile.GRID.STATVAR.MIDPOINTS > class_depths(i,1) & tile.GRID.STATVAR.MIDPOINTS <= class_depths(i+1,1));
%                 CURRENT.STATVAR.upperPos = tile.PARA.altitude - class_depths(i,1);
%                 CURRENT.STATVAR.lowerPos = tile.PARA.altitude - class_depths(i+1,1);
%                 for j=1:size(variables,1)
%                     if isfield(tile.GRID.STATVAR, variables{j,1})
%                         CURRENT.STATVAR.(variables{j,1}) = tile.GRID.STATVAR.(variables{j,1})(range);
%                     end
%                 end
%                 CURRENT=CURRENT.NEXT;
%             end
% 
%             %6. set top depths relative to surface and finalize initialization for
%             %subsurface classes
% 
%             CURRENT = tile.TOP_CLASS;
%             CURRENT.STATVAR.top_depth_rel2groundSurface = 0; %set initial surface to zero
% 
%             CURRENT = convert_units(CURRENT, tile);
%             CURRENT = finalize_init(CURRENT, tile);
% 
%             CURRENT.PARA.target_grid = tile.GRID.STATVAR.GRID;
%             CURRENT.PARA.target_layerThick = tile.GRID.STATVAR.layerThick;
%             while ~isequal(CURRENT.NEXT, tile.BOTTOM_CLASS.NEXT)
%                 CURRENT.NEXT.STATVAR.top_depth_rel2groundSurface = CURRENT.STATVAR.top_depth_rel2groundSurface + sum(CURRENT.STATVAR.layerThick,1);
% 
%                 CURRENT.NEXT = convert_units(CURRENT.NEXT, tile);
%                 CURRENT.NEXT = finalize_init(CURRENT.NEXT, tile);
% 
%                 CURRENT.NEXT.PARA.target_grid = tile.GRID.STATVAR.GRID;
%                 CURRENT.NEXT.PARA.target_layerThick = tile.GRID.STATVAR.layerThick;
% 
%                 CURRENT = CURRENT.NEXT;
%             end
% 
%             %7. assign interaction classes 
%             CURRENT = tile.TOP_CLASS;
% 
%             while ~isequal(CURRENT.NEXT, tile.BOTTOM)
%                 ia_class = get_IA_class(class(CURRENT), class(CURRENT.NEXT));
%                 CURRENT.IA_NEXT = ia_class;
%                 CURRENT.IA_NEXT.PREVIOUS = CURRENT;
%                 CURRENT.IA_NEXT.NEXT = CURRENT.NEXT;
%                 CURRENT.NEXT.IA_PREVIOUS = CURRENT.IA_NEXT;
% 
%                 finalize_init(CURRENT.IA_NEXT, tile); %Added Sebastian, defined in IA_BASE as empty
% 
%                 CURRENT = CURRENT.NEXT;
%             end
% 
%             %8. assign SNOW class
%             snow_class_name = strat_classes_class.PARA.snow_class_name;
%             snow_class_index = strat_classes_class.PARA.snow_class_index;
% 
%             if ~isempty(snow_class_name) && sum(isnan(snow_class_name))==0
%                 snow_class =  tile.RUN_INFO.PPROVIDER.CLASSES.(snow_class_name);
%                 snow_class = snow_class{snow_class_index,1};
% 
%                 %tile.TOP.STORE.SNOW = copy(snow_class);
%                 tile.STORE.SNOW = copy(snow_class);
%                 %no convert_units for SNOW classes needed at this point!
%                 %tile.TOP.STORE.SNOW = finalize_init(tile.TOP.STORE.SNOW, tile); %make this dependent on TILE!
%                 tile.STORE.SNOW = finalize_init(tile.STORE.SNOW, tile); %make this dependent on TILE!
%             end
% 
%             %9. assign sleeping classes
%             sleeping_classes = strat_classes_class.PARA.sleeping_classes_name;
%             sleeping_classes_index = strat_classes_class.PARA.sleeping_classes_index; 
% 
% %             for i=1:size(sleeping_classes,1)
% %                 sc = tile.RUN_INFO.PPROVIDER.CLASSES.(sleeping_classes{i,1});
% %                 sc = sc{sleeping_classes_index(i,1),1};
% %                 tile.TOP.STORE.SLEEPING{i,1} = copy(sc);
% %                 tile.TOP.STORE.SLEEPING{i,1} = convert_units(tile.TOP.STORE.SLEEPING{i,1}, tile);
% %                 tile.TOP.STORE.SLEEPING{i,1} = finalize_init(tile.TOP.STORE.SLEEPING{i,1}, tile);
% %                 tile.TOP.STORE.SLEEPING{i,2} = sleeping_classes_index(i,1);
% %             end
%             for i=1:size(sleeping_classes,1)
%                 sc = tile.RUN_INFO.PPROVIDER.CLASSES.(sleeping_classes{i,1});
%                 sc = sc{sleeping_classes_index(i,1),1};
%                 tile.STORE.SLEEPING{i,1} = copy(sc);
%                 tile.STORE.SLEEPING{i,1} = convert_units(tile.STORE.SLEEPING{i,1}, tile);
%                 tile.STORE.SLEEPING{i,1} = finalize_init(tile.STORE.SLEEPING{i,1}, tile);
%                 tile.STORE.SLEEPING{i,2} = sleeping_classes_index(i,1);
%             end
% 
%             %10. assign time, etc.
%             tile.t = tile.FORCING.PARA.start_time;
% 
%             %11. assign LATERAL classes 
%             tile.LATERAL = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.lateral_class){tile.PARA.lateral_class_index,1});
%             tile.LATERAL = finalize_init(tile.LATERAL, tile);
% 
%             %12. assign OUT classes
%             tile.OUT = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.out_class){tile.PARA.out_class_index,1});
%             tile.OUT = finalize_init(tile.OUT, tile);
% 
%         end
% 
%         function tile = build_tile_new_init_steady_state(tile)
% 
%             tile.PARA.run_name =  tile.RUN_INFO.PPROVIDER.PARA.run_name;
%             tile.PARA.result_path =  tile.RUN_INFO.PPROVIDER.PARA.result_path;
% 
%             %1. forcing
%             %tile.FORCING = copy(tile.RUN_INFO.PPROVIDER.FUNCTIONAL_CLASSES.FORCING{tile.PARA.forcing_index,1});
%             tile.FORCING = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.forcing_class){tile.PARA.forcing_class_index,1});
%             tile.FORCING = finalize_init(tile.FORCING, tile);
% 
%             %2. grid
%             %tile.GRID = copy(tile.RUN_INFO.PPROVIDER.FUNCTIONAL_CLASSES.GRID{tile.PARA.grid_index,1});
%             tile.GRID = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.grid_class){tile.PARA.grid_class_index,1});
%             tile.GRID = finalize_init(tile.GRID, tile);
% 
%             %3. map statvar to grid, using STRATIGRAPHY_STATVAR classes 
%             for i=1:size(tile.PARA.strat_statvar_class,1)
%                 strat_statvar_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.strat_statvar_class{i,1}){tile.PARA.strat_statvar_class_index(i,1),1});
%                 strat_statvar_class = finalize_init(strat_statvar_class, tile);
%             end
% 
%             %4. build stratigraphy
%             strat_classes_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.strat_classes_class){tile.PARA.strat_classes_class_index,1});
%             strat_classes_class = finalize_init(strat_classes_class, tile);
%             class_list = strat_classes_class.PARA.classes.class_name;
%             class_index = strat_classes_class.PARA.classes.class_index;
%             tile.TOP = Top();
%             CURRENT = tile.TOP;
%             for i=1:size(class_list,1)
%                 CURRENT.NEXT = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(class_list{i,1}){class_index(i,1)});
%                 CURRENT.NEXT.PREVIOUS = CURRENT;
%                 CURRENT = CURRENT.NEXT;
%             end
%             tile.BOTTOM = Bottom();
%             CURRENT.NEXT = tile.BOTTOM;
%             tile.BOTTOM.PREVIOUS = CURRENT;
% 
%             tile.TOP_CLASS = tile.TOP.NEXT;
%             tile.BOTTOM_CLASS = tile.BOTTOM.PREVIOUS;
% 
%             %5. assign STATVAR using STRATGRAPHY_STATVAR classes
%             class_depths = strat_classes_class.PARA.classes.depth;
%             class_depths = [class_depths; tile.GRID.STATVAR.GRID(end,1)];
% 
%             CURRENT = tile.TOP_CLASS;
%             for i=1:size(class_list,1)
%                 variables = fieldnames(CURRENT.STATVAR);
%                 range = (tile.GRID.STATVAR.MIDPOINTS > class_depths(i,1) & tile.GRID.STATVAR.MIDPOINTS <= class_depths(i+1,1));
%                 %CURRENT.STATVAR.layerThick = tile.GRID.STATVAR.LAYERTHICK(range,1);
%                 CURRENT.STATVAR.upperPos = tile.PARA.altitude - class_depths(i,1);
%                 CURRENT.STATVAR.lowerPos = tile.PARA.altitude - class_depths(i+1,1);
%                 for j=1:size(variables,1)
%                     if isfield(tile.GRID.STATVAR, variables{j,1})
%                         CURRENT.STATVAR.(variables{j,1}) = tile.GRID.STATVAR.(variables{j,1})(range);
%                     end
%                 end
%                 CURRENT=CURRENT.NEXT;
%             end
% 
%             %6. set top depths relative to surface and finalize initialization for
%             %subsurface classes
% 
%             %determine T_first_cell from class 
%             if ~isempty(tile.PARA.init_steady_state_class) && sum(isnan(tile.PARA.init_steady_state_class)) == 0
%                 init_steady_state_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.init_steady_state_class){tile.PARA.init_steady_state_class_index,1});
%                 init_steady_state_class = finalize_init(init_steady_state_class, tile); %assigns T_first_cell and potentially start_depth_steady_state
%             end
% 
%             CURRENT = tile.TOP_CLASS;
%             CURRENT.STATVAR.top_depth_rel2groundSurface = 0; %set initial surface to zero
% 
%             CURRENT.STATVAR.T = CURRENT.STATVAR.layerThick .*0 + tile.PARA.T_first_cell;            
%             CURRENT = convert_units(CURRENT, tile);
%             CURRENT = finalize_init(CURRENT, tile);
%             [CURRENT, T_end, thermCond_end, layerThick_end, start_depth_steady_state] = ...
%                 init_T_steady_state_TOP_CLASS(CURRENT, tile.PARA.T_first_cell, tile.PARA.start_depth_steady_state, tile.FORCING.PARA.heatFlux_lb);
% 
%             CURRENT.PARA.target_grid = tile.GRID.STATVAR.GRID;
%             CURRENT.PARA.target_layerThick = tile.GRID.STATVAR.layerThick;
%             while ~isequal(CURRENT.NEXT, tile.BOTTOM_CLASS.NEXT)
%                 CURRENT.NEXT.STATVAR.top_depth_rel2groundSurface = CURRENT.STATVAR.top_depth_rel2groundSurface + sum(CURRENT.STATVAR.layerThick,1);
% 
%                 CURRENT.NEXT.STATVAR.T = CURRENT.NEXT.STATVAR.layerThick .*0 + tile.PARA.T_first_cell;            
%                 CURRENT.NEXT = convert_units(CURRENT.NEXT, tile);
%                 CURRENT.NEXT = finalize_init(CURRENT.NEXT, tile);
%                 [CURRENT.NEXT, T_end, thermCond_end, layerThick_end, start_depth_steady_state] = ...
%                     init_T_steady_state(CURRENT.NEXT, T_end, start_depth_steady_state, thermCond_end, layerThick_end, tile.FORCING.PARA.heatFlux_lb);
% 
%                 CURRENT.NEXT.PARA.target_grid = tile.GRID.STATVAR.GRID;
%                 CURRENT.NEXT.PARA.target_layerThick = tile.GRID.STATVAR.layerThick;
% 
%                 CURRENT = CURRENT.NEXT;
%             end
% 
%             %7. assign interaction classes 
%             CURRENT = tile.TOP_CLASS;
% 
%             while ~isequal(CURRENT.NEXT, tile.BOTTOM)
%                 ia_class = get_IA_class(class(CURRENT), class(CURRENT.NEXT));
%                 CURRENT.IA_NEXT = ia_class;
%                 CURRENT.IA_NEXT.PREVIOUS = CURRENT;
%                 CURRENT.IA_NEXT.NEXT = CURRENT.NEXT;
%                 CURRENT.NEXT.IA_PREVIOUS = CURRENT.IA_NEXT;
% 
%                 finalize_init(CURRENT.IA_NEXT, tile); %Added Sebastian, defined in IA_BASE as empty
% 
%                 CURRENT = CURRENT.NEXT;
%             end
% 
%             %8. assign SNOW class
%             snow_class_name = strat_classes_class.PARA.snow_class_name;
%             snow_class_index = strat_classes_class.PARA.snow_class_index;
% 
%             if ~isempty(snow_class_name) && sum(isnan(snow_class_name))==0
%                 snow_class =  tile.RUN_INFO.PPROVIDER.CLASSES.(snow_class_name);
%                 snow_class = snow_class{snow_class_index,1};
% 
%                 tile.STORE.SNOW = copy(snow_class);
%                 %no convert_units for SNOW classes needed at this point!
%                 tile.STORE.SNOW = finalize_init(tile.STORE.SNOW, tile); %make this dependent on TILE!
%             end
% 
%             %9. assign sleeping classes
%             sleeping_classes = strat_classes_class.PARA.sleeping_classes_name;
%             sleeping_classes_index = strat_classes_class.PARA.sleeping_classes_index; 
% 
%             for i=1:size(sleeping_classes,1)
%                 sc = tile.RUN_INFO.PPROVIDER.CLASSES.(sleeping_classes{i,1});
%                 sc = sc{sleeping_classes_index(i,1),1};
%                 tile.STORE.SLEEPING{i,1} = copy(sc);
%                 tile.STORE.SLEEPING{i,1} = convert_units(tile.STORE.SLEEPING{i,1}, tile);
%                 tile.STORE.SLEEPING{i,1} = finalize_init(tile.STORE.SLEEPING{i,1}, tile);
%                 tile.STORE.SLEEPING{i,2} = sleeping_classes_index(i,1);
%             end
% 
%             %10. assign time, etc.
%             tile.t = tile.FORCING.PARA.start_time;
% 
%             %11. assign LATERAL classes 
%             tile.LATERAL = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.lateral_class){tile.PARA.lateral_class_index,1});
%             tile.LATERAL = finalize_init(tile.LATERAL, tile);
% 
%             %12. assign OUT classes
%             tile.OUT = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.out_class){tile.PARA.out_class_index,1});
%             tile.OUT = finalize_init(tile.OUT, tile);
% 
%         end
% 
%         function tile = build_tile_update_stratigraphy(tile)
% 
%             old_TILE = tile.RUN_INFO.TILE;   
%             old_TILE = truncate_stratigraphy(old_TILE, tile.PARA.truncate_depth);
% 
%             tile = build_tile_new_init(tile);
%             tile = merge_stratigraphies(tile, old_TILE);
% 
%             if tile.PARA.keep_LATERAL
%                 tile.LATERAL = old_TILE.LATERAL;
%             end
% 
%             tile.TEMP.time_difference = tile.RUN_INFO.TILE.t - tile.FORCING.PARA.start_time; %Used to correct time variables in subsurface classes 
%             tile.t = tile.FORCING.PARA.start_time;
%             tile.LATERAL.IA_TIME = tile.t + tile.LATERAL.IA_TIME_INCREMENT;
%             tile.LATERAL = update_lateral(tile.LATERAL, tile);
%             tile.next_break_time =  tile.LATERAL.IA_TIME; 
%            % tile.next_break_time = old_TILE.next_break_time;
% 
%             CURRENT = tile.TOP.NEXT;
%             while ~isequal(CURRENT.NEXT, tile.BOTTOM)
%                 CURRENT = reset_timestamps(CURRENT, tile);
%                 CURRENT = CURRENT.NEXT;
%             end
% 
%             tile.LATERAL.TOP = tile.TOP;
%             tile.LATERAL.BOTTOM = tile.BOTTOM;
% 
%             tile.RUN_INFO.TILE = tile;
% 
%             tile.PARA.run_name =  tile.RUN_INFO.PPROVIDER.PARA.run_name;
%             tile.PARA.result_path =  tile.RUN_INFO.PPROVIDER.PARA.result_path;
% 
%         end
% 
% 
%         function tile = build_tile_update_forcing_out(tile)
% 
%             %2. forcing -> special forcing class required
%             tile.FORCING = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.forcing_class){tile.PARA.forcing_class_index,1});
% 
%             %12. assign OUT classes
%             tile.OUT = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.out_class){tile.PARA.out_class_index,1});
% 
%             %use old tile
%             tile.GRID = tile.RUN_INFO.TILE.GRID;   
%             tile.BOTTOM = tile.RUN_INFO.TILE.BOTTOM;
%             tile.TOP = tile.RUN_INFO.TILE.TOP;
%             tile.BOTTOM_CLASS = tile.RUN_INFO.TILE.BOTTOM_CLASS;
%             tile.TOP_CLASS = tile.RUN_INFO.TILE.TOP_CLASS;
%             tile.timestep = tile.RUN_INFO.TILE.timestep;
%             tile.LATERAL = tile.RUN_INFO.TILE.LATERAL; 
%             % tile.LATERAL.IA_TIME = tile.FORCING.PARA.start_time + tile.LATERAL.IA_TIME_INCREMENT;
%            % tile.next_break_time = tile.RUN_INFO.TILE.next_break_time; 
%             tile.STORE = tile.RUN_INFO.TILE.STORE;     
% 
%             %use old PARA, but overwrite all newly set values
%             PARA_new = tile.PARA;
%             tile.PARA = tile.RUN_INFO.TILE.PARA;
%             fn = fieldnames(PARA_new);
%             for i=1:size(fn,1)  %be careful, does not work if empty array (and not NaN) is willingly assigned to a parameter
%                 if ~isempty(PARA_new.(fn{i,1}))
%                     tile.PARA.(fn{i,1}) = PARA_new.(fn{i,1});
%                 end
%             end
% 
%             tile.FORCING = finalize_init(tile.FORCING, tile); 
%             tile.OUT = finalize_init(tile.OUT, tile);           
%             %10. assign time, etc.
%             tile.TEMP.time_difference = tile.RUN_INFO.TILE.t - tile.FORCING.PARA.start_time; %Used to correct time variables in subsurface classes 
%             tile.t = tile.FORCING.PARA.start_time;
% 
%             %reset IA time
%             tile.LATERAL.IA_TIME = tile.t + tile.LATERAL.IA_TIME_INCREMENT;
%             tile.LATERAL = update_lateral(tile.LATERAL, tile);
%             tile.next_break_time =  tile.LATERAL.IA_TIME; 
% 
%             %reset time for BGC class (do mothing if no BGC class exists)
%             %-> MAKE THIS A GENERAL RESET_TIME OR ADJUST_TIME FUNCTION THAT
%             %IS DEFINED IN BASE AND OVERWRITTEN IN ALL FUNCTIONS THAT
%             %ACTUALLY HAVE A TIME VARIABLE - CALCULATE TIME OFFSET BETWEEN
%             %OLD AND NEW FORCING, I.E LAST TIMESTAMP OF OLD RUN AND FIRST TIMESTAMP OF NEW RUN 
%             CURRENT = tile.TOP.NEXT;
%             while ~isequal(CURRENT.NEXT, tile.BOTTOM)
%                 CURRENT = reset_timestamps(CURRENT, tile);
%                 CURRENT = CURRENT.NEXT;
%             end
% 
% 
%             tile.RUN_INFO.TILE = tile;
% 
%             tile.PARA.run_name =  tile.RUN_INFO.PPROVIDER.PARA.run_name;
%             tile.PARA.result_path =  tile.RUN_INFO.PPROVIDER.PARA.result_path;
% 
% 
%             if ~isempty(tile.PARA.modify_class) && sum(isnan(tile.PARA.modify_class_index))==0
%                 for i=1:size(tile.PARA.modify_class,1)
%                     mod = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.modify_class{i,1}){tile.PARA.modify_class_index(i,1),1});
%                     mod = finalize_init(mod, tile);
%                     tile = modify(mod, tile);
%                 end
%             end
% 
%         end
% 
% 
%         function tile = build_tile_restart_OUT_all(tile)
% 
%             tile.PARA.run_name =  tile.RUN_INFO.PPROVIDER.PARA.run_name;
%             tile.PARA.result_path =  tile.RUN_INFO.PPROVIDER.PARA.result_path;
% 
%             %1. forcing
%             %tile.FORCING = copy(tile.RUN_INFO.PPROVIDER.FUNCTIONAL_CLASSES.FORCING{tile.PARA.forcing_index,1});
%             tile.FORCING = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.forcing_class){tile.PARA.forcing_class_index,1});
%             tile.FORCING = finalize_init(tile.FORCING, tile);
% 
%             %4. build stratigraphy
%             strat_classes_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.strat_classes_class){tile.PARA.strat_classes_class_index,1});
%             strat_classes_class = finalize_init(strat_classes_class, tile);
% 
%             tile.TOP = Top();
%             tile.BOTTOM = Bottom();
%             OUT = load([tile.PARA.restart_file_path tile.PARA.restart_file_name]);
%             OUT = OUT.out;
%             if ~isempty(tile.PARA.restart_time) && sum(isnan(tile.PARA.restart_time))==0
%                 tile.PARA.restart_time = datenum(tile.PARA.restart_time(1,1), tile.PARA.restart_time(2,1), tile.PARA.restart_time(3,1));
%             else
%                 tile.PARA.restart_time = [];
%             end
%             [uppermost_class, lowermost_class] = restore_stratigraphy_from_OUT(OUT, tile.PARA.restart_time, tile);
%             tile.TOP.NEXT = uppermost_class;
%             tile.TOP.NEXT.PREVIOUS = tile.TOP;
%             tile.BOTTOM.PREVIOUS = lowermost_class;
%             tile.BOTTOM.PREVIOUS.NEXT = tile.BOTTOM;
% 
%             tile.TOP_CLASS = tile.TOP.NEXT;
%             tile.BOTTOM_CLASS = tile.BOTTOM.PREVIOUS;
% 
%             %6. set top depths relative to surface and finalize initialization for
%             %subsurface classes
%             CURRENT = tile.TOP_CLASS;
%             CURRENT.STATVAR.top_depth_rel2groundSurface = 0; %set initial surface to zero
% 
%             %7. assign interaction classes 
%             CURRENT = tile.TOP_CLASS;
% 
%             while ~isequal(CURRENT.NEXT, tile.BOTTOM)
%                 ia_class = get_IA_class(class(CURRENT), class(CURRENT.NEXT));
%                 CURRENT.IA_NEXT = ia_class;
%                 CURRENT.IA_NEXT.PREVIOUS = CURRENT;
%                 CURRENT.IA_NEXT.NEXT = CURRENT.NEXT;
%                 CURRENT.NEXT.IA_PREVIOUS = CURRENT.IA_NEXT;
% 
%                 finalize_init(CURRENT.IA_NEXT, tile); %Added Sebastian, defined in IA_BASE as empty
% 
%                 CURRENT = CURRENT.NEXT;
%             end
% 
%             %8. assign SNOW class
%             snow_class_name = strat_classes_class.PARA.snow_class_name;
%             snow_class_index = strat_classes_class.PARA.snow_class_index;
% 
%             if ~isempty(snow_class_name) && sum(isnan(snow_class_name))==0
%                 snow_class =  tile.RUN_INFO.PPROVIDER.CLASSES.(snow_class_name);
%                 snow_class = snow_class{snow_class_index,1};
%                 tile.STORE.SNOW = copy(snow_class);
%                 tile.STORE.SNOW = finalize_init(tile.STORE.SNOW, tile); %make this dependent on TILE!
%             end
% 
%             %9. assign sleeping classes
%             sleeping_classes = strat_classes_class.PARA.sleeping_classes_name;
%             sleeping_classes_index = strat_classes_class.PARA.sleeping_classes_index; 
% 
%             for i=1:size(sleeping_classes,1)
%                 sc = tile.RUN_INFO.PPROVIDER.CLASSES.(sleeping_classes{i,1});
%                 sc = sc{sleeping_classes_index(i,1),1};
%                 tile.STORE.SLEEPING{i,1} = copy(sc);
%                 tile.STORE.SLEEPING{i,1} = convert_units(tile.STORE.SLEEPING{i,1}, tile);
%                 tile.STORE.SLEEPING{i,1} = finalize_init(tile.STORE.SLEEPING{i,1}, tile);
%                 tile.STORE.SLEEPING{i,2} = sleeping_classes_index(i,1);
%             end
% 
%             %10. assign time, etc.
%             tile.t = tile.FORCING.PARA.start_time;
% 
%             %11. assign LATERAL classes 
%             tile.LATERAL = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.lateral_class){tile.PARA.lateral_class_index,1});
%             tile.LATERAL = finalize_init(tile.LATERAL, tile);
% 
%             %12. assign OUT classes
%             tile.OUT = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.out_class){tile.PARA.out_class_index,1});
%             tile.OUT = finalize_init(tile.OUT, tile);
% 
%         end
% 
%         function tile = build_tile_restart_OUT_last_timestep(tile)
%             new_run_name = tile.RUN_INFO.PPROVIDER.PARA.run_name;
%             new_result_path = tile.RUN_INFO.PPROVIDER.PARA.result_path;
%             adapt_run_name = tile.PARA.adapt_run_name;
%             new_end_time = tile.PARA.new_end_time;
%             if ~isempty(new_end_time) && sum(isnan(new_end_time))==0
%                 new_end_time = datenum(new_end_time(1,1), new_end_time(2,1), new_end_time(3,1));
%             else
%                 new_end_time = [];
%             end
% 
%             temp=load([tile.PARA.restart_file_path tile.PARA.restart_file_name]);
%             variables = fieldnames(temp.out.STRATIGRAPHY);
%             for i=1:size(variables,1)
%                 tile.(variables{i,1}) = temp.out.STRATIGRAPHY.(variables{i,1});
%             end
%             if ~isempty(adapt_run_name) && ~isnan(adapt_run_name)
%                 if adapt_run_name
%                     tile.PARA.run_name = new_run_name;
%                     tile.PARA.result_path = new_result_path;
%                 end
%             end
%             if ~isempty(new_end_time)
%                 tile.FORCING.PARA.end_time = new_end_time;
%             end
% 
% %             tile.LATERAL.IA_CLASSES = {};
% %             tile.LATERAL.PARA.num_realizations = 1;
% %             tile.LATERAL.PARA.worker_number = 1;
% %             tile.OUT.OUTPUT_TIME = tile.OUT.OUTPUT_TIME+100;
%             %tile.LATERAL.IA_TIME = tile.FORCING.PARA.end_time;
%         end
% 
% 
%         %-----------------------------------
%         %service functions at TILE level
% 
%         function tile = truncate_stratigraphy(tile, truncate_depth)
%             CURRENT = tile.TOP.NEXT;
%             while ~is_ground_surface(CURRENT)
%                 CURRENT = CURRENT.NEXT;
%             end
%             %CURRENT is now class that has ground surface
%             while ~isequal(CURRENT, tile.BOTTOM)
% 
%                 i=1;
%                 while truncate_depth > 0 && i<=size(CURRENT.STATVAR.layerThick,1)
%                     truncate_depth = truncate_depth - CURRENT.STATVAR.layerThick(i,1);
%                     i=i+1;
%                 end
%                 i=i-1;
% 
%                 if truncate_depth <= 1e-6 %truncate inside this class
%                     if i==0
%                         %do nothing, keep entire class
%                     elseif i==size(CURRENT.STATVAR.layerThick,1)
%                         %remove entire class
%                         CURRENT = CURRENT.NEXT;
%                     else
%                         %call function at stratgraphy class level and
%                         %truncate the class
%                         CURRENT= truncate_STATVAR(CURRENT, i, truncate_depth);
%                         CURRENT = compute_diagnostic(CURRENT, tile);
%                     end
%                     break
%                 end
%                 CURRENT = CURRENT.NEXT;
%             end
%             tile.TOP.NEXT = CURRENT;
%             tile.TOP.NEXT.PREVIOUS = tile.TOP;
%             tile.TOP.NEXT.IA_PREVIOUS = [];
%         end
% 
% 
%         function tile = merge_stratigraphies(tile, tile_below)
%              if tile.PARA.domain_depth == 0 %delete empty stratigraphy
%                 tile.TOP.NEXT = tile.BOTTOM;
%                 tile.BOTTOM.PREVIOUS = tile.TOP;
%              end
% 
%             if strcmp(class(tile.BOTTOM.PREVIOUS), class(tile_below.TOP.NEXT))
%                 %call function at stratgraphy class level and merge classes
%                 tile.BOTTOM.PREVIOUS = merge_STATVAR(tile.BOTTOM.PREVIOUS, tile_below.TOP.NEXT);
%                 tile.BOTTOM.PREVIOUS = compute_diagnostic(tile.BOTTOM.PREVIOUS, tile);
%                 CURRENT_BELOW = tile_below.TOP.NEXT.NEXT;
%             else
%                 CURRENT_BELOW = tile_below.TOP.NEXT;
%             end
%             %connect classes in between
%             CURRENT_BELOW.PREVIOUS = tile.BOTTOM.PREVIOUS;
%             CURRENT_BELOW.PREVIOUS.NEXT = CURRENT_BELOW;
%             ia_class = get_IA_class(class(CURRENT_BELOW.PREVIOUS), class(CURRENT_BELOW));
%             CURRENT_BELOW.PREVIOUS.IA_NEXT = ia_class;
%             CURRENT_BELOW.IA_PREVIOUS = ia_class;
%             CURRENT_BELOW.PREVIOUS.IA_NEXT.NEXT = CURRENT_BELOW;
%             CURRENT_BELOW.IA_PREVIOUS.PREVIOUS = CURRENT_BELOW.PREVIOUS;
% 
%             finalize_init(CURRENT_BELOW.IA_PREVIOUS, tile); 
% 
%             tile.BOTTOM.PREVIOUS = tile_below.BOTTOM.PREVIOUS;
%             tile.BOTTOM.PREVIOUS.NEXT = tile.BOTTOM;
% 
%         end
% 
% 
% 
%         %-------------param file generation-----
%        function tile = param_file_info(varargin)
% 
%             if nargin==2
%                 tile = varargin{1};
%                 option = varargin{2};
%             else
%                 tile = varargin{1};
%                 option = 0;
%             end
%             tile.PARA.STATVAR = [];
%             tile.PARA.class_category = 'TILE';
%             tile.PARA.default_value=[];
%             tile.PARA.options=[];
%             tile.PARA.comment=[];
% 
%             if strcmp(option, 'new_init')
%                 tile.PARA.builder = [];
%                 tile.PARA.default_value.builder = {'new_init'};
%                 parameters = { 'domain_depth';  'forcing_class'; 'forcing_class_index'; 'grid_class'; 'grid_class_index'; 'out_class'; ...
%                     'out_class_index'; 'strat_classes_class'; 'strat_classes_class_index'; 'strat_statvar_class'; 'strat_statvar_class_index'; 'lateral_class'; ...
%                     'lateral_class_index'; 'lateral_IA_classes'; 'lateral_IA_classes_index'};  %'latitude'; 'longitude'; 'altitude'; 'area';
% 
%                 for i=1:size(parameters,1)
%                     tile.PARA.(parameters{i,1})=[];
%                 end
% 
% 
%                 tile.PARA.default_value.domain_depth = {100};
%                 tile.PARA.comment.domain_depth = {'vertical depth of the model domain [m]'};
% 
%                 tile.PARA.default_value.forcing_class = {'FORCING_seb'};
%                 tile.PARA.default_value.forcing_class_index = {1};
%                 tile.PARA.default_value.grid_class = {'GRID_user_defined'};
%                 tile.PARA.default_value.grid_class_index = {1};
%                 tile.PARA.default_value.out_class = {'OUT_all_lateral'};
%                 tile.PARA.default_value.out_class_index = {1};
%                 tile.PARA.default_value.strat_classes_class = {'STRAT_classes'};
%                 tile.PARA.default_value.strat_classes_class_index = {1};
% 
%                 tile.PARA.comment.strat_statvar_class = {'list of STRATIGRAPHY_STATVAR classes that provide initial state of state variables'};
%                 tile.PARA.options.strat_statvar_class.name =  'H_LIST'; %
%                 tile.PARA.options.strat_statvar_class.entries_x = {'STRAT_layers' 'STRAT_linear'};
%                 tile.PARA.options.strat_statvar_class_index.name =  'H_LIST'; 
%                 tile.PARA.options.strat_statvar_class_index.entries_x = {1 1};
% 
%                 tile.PARA.comment.lateral_class = {'lateral class, e.g. LATERAL1D or LATERAL3D'};
%                 tile.PARA.default_value.lateral_class = {'LATERAL_1D'};
%                 tile.PARA.default_value.lateral_class_index = {1};
% 
%                 tile.PARA.comment.lateral_IA_classes = {'list of lateral interaction classes'};
%                 tile.PARA.options.lateral_IA_classes.name =  'H_LIST'; %
%                 tile.PARA.options.lateral_IA_classes_index.name =  'H_LIST';
% 
%             elseif strcmp(option, 'new_init_steady_state')
%                 tile.PARA.builder = [];
%                 tile.PARA.default_value.builder = {'new_init_steady_state'};
%                 parameters = {'domain_depth'; 'forcing_class'; 'forcing_class_index'; 'grid_class'; 'grid_class_index'; 'out_class'; ...
%                     'out_class_index'; 'strat_classes_class'; 'strat_classes_class_index'; 'strat_statvar_class'; 'strat_statvar_class_index'; 'lateral_class'; ...
%                     'lateral_class_index'; 'lateral_IA_classes'; 'lateral_IA_classes_index'; 'init_steady_state_class'; 'init_steady_state_class_index';...
%                     'T_first_cell'; 'start_depth_steady_state'};
%                 for i=1:size(parameters,1)
%                     tile.PARA.(parameters{i,1})=[];
%                 end
% 
%                 tile.PARA.default_value.domain_depth = {100};
%                 tile.PARA.comment.domain_depth = {'vertical depth of the model domain [m]'};
% 
%                 tile.PARA.default_value.forcing_class = {'FORCING_seb'};
%                 tile.PARA.default_value.forcing_class_index = {1};
%                 tile.PARA.default_value.grid_class = {'GRID_user_defined'};
%                 tile.PARA.default_value.grid_class_index = {1};
%                 tile.PARA.default_value.out_class = {'OUT_all_lateral'};
%                 tile.PARA.default_value.out_class_index = {1};
%                 tile.PARA.default_value.strat_classes_class = {'STRAT_classes'};
%                 tile.PARA.default_value.strat_classes_class_index = {1};
% 
%                 tile.PARA.comment.strat_statvar_class = {'list of STRATIGRAPHY_STATVAR classes that provide initial state of state variables'};
%                 tile.PARA.options.strat_statvar_class.name =  'H_LIST'; %
%                 tile.PARA.options.strat_statvar_class.entries_x = {'STRAT_layers'};
%                 tile.PARA.options.strat_statvar_class_index.name =  'H_LIST'; 
%                 tile.PARA.options.strat_statvar_class_index.entries_x = {1};
% 
%                 tile.PARA.comment.lateral_class = {'lateral class, e.g. LATERAL1D or LATERAL3D'};
%                 tile.PARA.default_value.lateral_class = {'LATERAL_1D'};
%                 tile.PARA.default_value.lateral_class_index = {1};
% 
%                 tile.PARA.comment.lateral_IA_classes = {'list of lateral interaction classes'};
%                 tile.PARA.options.lateral_IA_classes.name =  'H_LIST'; %
%                 tile.PARA.options.lateral_IA_classes_index.name =  'H_LIST';
% 
%                 tile.PARA.comment.init_steady_state_class = {'init_steady_state class to compute temperature of first grid cell, leave empty when using T_first_grid_cell'};
%                 tile.PARA.comment.T_first_cell = {'temperature of first grid cell used to compute the temperature gradient, leave empty when using init_steady_state class'};
%                 tile.PARA.comment.start_depth_steady_state = {'depth [m] where temperature gradient starts, constant above, leave empty when using init_steady_state class'};
% 
%             elseif strcmp(option, 'update_forcing_out')
%                 tile.PARA.builder = [];
%                 tile.PARA.default_value.builder = {'update_forcing_out'};
%                 parameters = { 'forcing_class'; 'forcing_class_index';  'out_class'; 'out_class_index'};
% 
%                 for i=1:size(parameters,1)
%                     tile.PARA.(parameters{i,1})=[];
%                 end
% 
%                 tile.PARA.default_value.forcing_class = {'FORCING_seb'};
%                 tile.PARA.default_value.forcing_class_index = {1};
%                 tile.PARA.default_value.out_class = {'OUT_all_lateral'};
%                 tile.PARA.default_value.out_class_index = {1};
% 
%             elseif strcmp(option, 'restart_OUT_last_timestep')
%                 tile.PARA.builder = [];
%                 tile.PARA.default_value.builder = {'restart_OUT_last_timestep'};
%                 parameters = {'restart_file_path'; 'restart_file_name'};
%                 for i=1:size(parameters,1)
%                     tile.PARA.(parameters{i,1})=[];
%                 end
% 
%                 tile.PARA.comment.restart_file_path = {'path and filename of restart file'};
%             else
%                 tile = provide_PARA(tile);
%             end
% 
%         end
% 

    end
end



