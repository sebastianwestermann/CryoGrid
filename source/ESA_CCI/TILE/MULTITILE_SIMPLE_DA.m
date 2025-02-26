
classdef MULTITILE_SIMPLE_DA < matlab.mixin.Copyable
    
    properties
        PARA
        RUN_INFO
        BUILDER
        FORCING
        CONST
        GRID
        OUT        
        SUBSURFACE_CLASS
        ENSEMBLE
        STRATIGRAPHY
        STORE
        
        t        
        timestep
        
        TOP
        BOTTOM
        DA
    end
    
    
    methods
        
       

        function tile = provide_PARA(tile)
            
            tile.PARA.builder = [];
            tile.PARA.number_of_realizations = [];
            tile.PARA.ensemble_size = [];
            tile.PARA.ensemble_class =[];
            tile.PARA.ensemble_class_index =[];
            tile.PARA.subsurface_class =[];
            tile.PARA.subsurface_class_index =[];
            tile.PARA.forcing_class = [];
            tile.PARA.forcing_class_index = [];
            tile.PARA.out_class = [];
            tile.PARA.out_class_index = [];
            
            tile.PARA.DA_classes = [];
            tile.PARA.DA_classes_index = [];
        end
        
        function tile = provide_CONST(tile)

            tile.CONST.day_sec = [];

        end
        
        function tile = provide_STATVAR(tile)

        end
        
        
        %assemble the stratigraphy
        function tile = finalize_init(tile)
            
            builder = str2func(tile.PARA.builder);
            tile.BUILDER = builder();
            tile.BUILDER.TILE = tile;
            
            build_tile(tile.BUILDER);
           
        end
        

        function tile = interpolate_forcing_tile(tile)
             tile.FORCING = interpolate_forcing(tile.FORCING, tile);
        end

        
        function tile = store_OUT_tile(tile)
            tile.OUT = store_OUT(tile.OUT, tile);
        end        
        
        function tile = data_assimilation(tile)
            for i = 1:size(tile.DA, 1)
                tile.DA{i,1} = DA_step(tile.DA{i,1}, tile);
                %recalculate stratigraphy
%                 if tile.DA{i,1}.TEMP.recalculate_stratigraphy_now == 1 %set in DA class
%                     tile.DA{i,1}.TEMP.recalculate_stratigraphy_now = 0;
%                     tile = recalculate_stratigraphy(tile);
%                     %reset observation operators if necessary
%                     for j=1:size(tile.DA, 1)
%                         for k=1:size(tile.DA{j,1}.OBS_OP,1)
%                             tile.DA{j,1}.OBS_OP{k,1} = reset_new_stratigraphy(tile.DA{j,1}.OBS_OP{k,1}, tile);
%                         end
%                     end
%                 end
%                 if tile.DA{i,1}.TEMP.recalculate_forcing_now == 1 %set in DA class
%                     tile.DA{i,1}.TEMP.recalculate_forcing_now = 0;
%                     tile = recalculate_forcing(tile);
%                 end
            end
        end
        
        
        function tile = recalculate_forcing(tile)
                tile.FORCING = adjust_forcing(tile.FORCING, tile);
        end
        
        
        function tile = run_model(tile)
            
            %=========================================================================
            %TIME INTEGRATION
            %=========================================================================
            while tile.t < tile.FORCING.PARA.end_time
                
                CURRENT = tile.SUBSURFACE_CLASS;
                
                %interpolate focing data to time t
                tile = interpolate_forcing_tile(tile);
                
                %upper boundar condition (uppermost class only)
                CURRENT = get_boundary_condition_u(CURRENT, tile);
                
                %lower boundary condition (lowermost class)
                CURRENT = get_boundary_condition_l(CURRENT,  tile);
                
                %calculate spatial derivatives
                CURRENT = get_derivatives_prognostic(CURRENT, tile);

                %prognostic step - integrate prognostic variables in time
                CURRENT = advance_prognostic(CURRENT, tile);
                
                %diagnostic step - compute diagnostic variables
                CURRENT = compute_diagnostic(CURRENT, tile);

                %triggers
                CURRENT = check_trigger(CURRENT, tile);

                %update time variable t
                tile.t = tile.t + tile.timestep./tile.CONST.day_sec;
                
                %model output
                tile = store_OUT_tile(tile);
                
                %data assimilation
                tile = data_assimilation(tile);
            end
            
        end
        
        
        %---BUILDER functions--------------
        
        function tile = build_tile_new_init(tile)
            
            tile.PARA.run_name =  tile.RUN_INFO.PPROVIDER.PARA.run_name;
            tile.PARA.result_path =  tile.RUN_INFO.PPROVIDER.PARA.result_path;
            tile.PARA.tile_size = tile.PARA.number_of_realizations .* tile.PARA.ensemble_size;

            tile.ENSEMBLE = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.ensemble_class){tile.PARA.ensemble_class_index});
            tile.ENSEMBLE = finalize_init(tile.ENSEMBLE, tile);
            
                        
            %2. forcing -> special forcing class required
            tile.timestep = tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.subsurface_class){tile.PARA.subsurface_class_index,1}.PARA.timestep;
            tile.FORCING = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.forcing_class){tile.PARA.forcing_class_index,1});
            tile.FORCING = finalize_init(tile.FORCING, tile);


            %1. build stratigraphy
            tile.SUBSURFACE_CLASS = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.subsurface_class){tile.PARA.subsurface_class_index,1});
            tile.SUBSURFACE_CLASS = finalize_init(tile.SUBSURFACE_CLASS, tile);
            
            %set this for being able to use existing OUT classes
            tile.TOP = Top();
            tile.BOTTOM = Bottom();
            tile.TOP.NEXT = tile.SUBSURFACE_CLASS;
            tile.TOP.NEXT.PREVIOUS = tile.TOP;
            tile.SUBSURFACE_CLASS.NEXT = tile.BOTTOM();
            tile.BOTTOM.PREVIOUS = tile.SUBSURFACE_CLASS;            
            
            %set fixed timestep!
            tile.timestep = tile.SUBSURFACE_CLASS.PARA.timestep;

             %10. assign time, etc.
            tile.t = tile.FORCING.PARA.start_time;

            %12. assign OUT classes
            tile.OUT = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.out_class){tile.PARA.out_class_index,1});
            tile.OUT = finalize_init(tile.OUT, tile);
            
            tile.RUN_INFO.TILE = tile;
            
            %assign DA classes
            for i = 1:size(tile.PARA.DA_classes,1)
              tile.DA{i,1} = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.DA_classes{i,1}){tile.PARA.DA_classes_index(i,1),1});  
              tile.DA{i,1} = finalize_init(tile.DA{i,1}, tile);
            end
            
        end

        function tile = build_tile_update_forcing_out(tile)
                        
            %2. forcing -> special forcing class required
            tile.FORCING = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.forcing_class){tile.PARA.forcing_class_index,1});

            %12. assign OUT classes
            tile.OUT = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.out_class){tile.PARA.out_class_index,1});
            
            %use old tile
            tile.ENSEMBLE = tile.RUN_INFO.TILE.ENSEMBLE;
            tile.SUBSURFACE_CLASS = tile.RUN_INFO.TILE.SUBSURFACE_CLASS;
            tile.BOTTOM = tile.RUN_INFO.TILE.BOTTOM;
            tile.TOP = tile.RUN_INFO.TILE.TOP;
            tile.timestep = tile.RUN_INFO.TILE.timestep;
            
            %use old PARA, but overwrite all newly set values
            PARA_new = tile.PARA;
            tile.PARA = tile.RUN_INFO.TILE.PARA;
            fn = fieldnames(PARA_new);
            for i=1:size(fn,1)
                if ~isempty(PARA_new.(fn{i,1}))
                    tile.PARA.(fn{i,1}) = PARA_new.(fn{i,1});
                end
            end

            
            tile.FORCING = finalize_init(tile.FORCING, tile); 
            tile.OUT = finalize_init(tile.OUT, tile);           
            %10. assign time, etc.
            tile.t = tile.FORCING.PARA.start_time;
            
            tile.RUN_INFO.TILE = tile;
            
        end

    end
end



