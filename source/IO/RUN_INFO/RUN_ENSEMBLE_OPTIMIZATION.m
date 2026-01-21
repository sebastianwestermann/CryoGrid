%========================================================================
% CryoGrid RUN_INFO class RUN_ENSEMBLE_OPTIMIZATION
% RUN_INFO class for spatially distributed runs (using an appropriate 
% SPATIAL_REFERENCE class, DATA_PROVIDER classes and FORCING classes)
% which can run several TILE classes per point sequentially for model spin-up 
%
% S. Westermann, Dec 2025
%========================================================================

classdef RUN_ENSEMBLE_OPTIMIZATION < matlab.mixin.Copyable
    
    properties
        PPROVIDER
        PARA
        CONST
        STATVAR
        TILE
        TEMP
        SPATIAL
        ENSEMBLE
        OPT
    end
    
    %SPATIAL contains ensemble parameters
    %OPT contains the classes used for optimization, e.g. DA
    methods
        
        function run_info = provide_PARA(run_info)
            
            run_info.PARA.run_mode = 'TILE';
            run_info.PARA.max_number_of_gridcells = Inf;
            run_info.PARA.parallel_pool_open = 0;

            run_info.PARA.number_of_cores = [];
            run_info.PARA.number_of_OPT_threads = [];
            
            run_info.PARA.tile_class = [];
            run_info.PARA.tile_class_index = [];
            %run_info.PARA.number_of_runs_per_tile = []; %vector

            run_info.PARA.projection_class = [];
            run_info.PARA.projection_class_index = [];
            
            run_info.PARA.ensemble_class = []; %same as projection_class, must set the spatial attributes as well 
            run_info.PARA.ensemble_class_index = [];

            run_info.PARA.optimization_class = [];
            run_info.PARA.optimization_class_index = [];                        
        end
        
        function run_info = provide_CONST(run_info)

        end
        
        function run_info = provide_STATVAR(run_info)

        end
        
        
        function run_info = finalize_init(run_info)
        
            if ~isempty(run_info.PARA.projection_class) && ~(sum(isnan(run_info.PARA.projection_class))>0)
                disp('get spatial data')
                spatial_class = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.projection_class){run_info.PARA.projection_class_index,1});
                if ~spatial_class.PARA.new_reference
                    spatial_class.STATVAR = run_info.SPATIAL.STATVAR;
                end
                run_info.SPATIAL = spatial_class;
                run_info.SPATIAL.RUN_INFO = run_info;
                run_info.SPATIAL = finalize_init(run_info.SPATIAL);
            end

            if sum(isnan(run_info.PARA.max_number_of_gridcells)) > 1 || ischar(run_info.PARA.max_number_of_gridcells)
                run_info.PARA.max_number_of_gridcells = Inf;
            end
        end

        
        function [run_info, tile] = run_model(run_info)
            if strcmp(run_info.PARA.run_mode, 'TILE')
                [run_info, tile] = run_model_TILE(run_info);
            elseif strcmp(run_info.PARA.run_mode, 'MULTITILE')
                [run_info, tile] = run_model_MULTITILE(run_info);
            end
        end

        function [run_info, tile] = run_model_TILE(run_info)
            tile = 0;
            if run_info.PARA.number_of_cores > 1
                if ~run_info.PARA.parallel_pool_open
                    parpool(run_info.PARA.number_of_cores)
                    spmd
                        [run_info, tile] = run_model_TILE_parallel(run_info);
                    end
                    delete(gcp('nocreate'));
                else
                    [run_info, tile] = run_model_TILE_parallel(run_info);
                end
            else
                [run_info, tile] = run_model_TILE_sequential(run_info);
            end
        end

        function [run_info, tile] = run_model_MULTITILE(run_info)
            if run_info.PARA.number_of_cores > 1
                if ~run_info.PARA.parallel_pool_open
                    parpool(run_info.PARA.number_of_cores)
                    spmd
                        [run_info, tile] = run_model_MULTITILE_parallel(run_info);
                    end
                    delete(gcp('nocreate'));
                else
                    [run_info, tile] = run_model_MULTITILE_parallel(run_info);
                end
            else
                [run_info, tile] = run_model_MULTITILE_sequential(run_info);
            end
        end

        function [run_info, tile] = run_model_TILE_parallel(run_info)
            tile = 0;

            % number of cores: run_info.PARA.number_of_cores
            %number of realizations: run_info.ENSEMBLE.TEMP.ensemble_size
            % if run_info.ENSEMBLE.TEMP.ensemble_size>=run_info.PARA.number_of_cores
            % -> run the DA on all cores and have sequential runs for the
            % non-da run in space
            %else ->split up the DA runs on several spmds

            worker_number = spmdIndex();
            number_of_cores_per_OPT_thread = []; 
            remaining_cores = run_info.PARA.number_of_cores;
            for i=1:run_info.PARA.number_of_OPT_threads
                number_of_cores_per_OPT_thread = [number_of_cores_per_OPT_thread; round(remaining_cores ./ (run_info.PARA.number_of_OPT_threads-i+1))];
                remaining_cores = remaining_cores - number_of_cores_per_OPT_thread(end,1);
            end

            OPT_thread_number=[]; %for each core gives the number of the OPT thread
            OPT_worker_number = [];
            for i=1:run_info.PARA.number_of_OPT_threads
                OPT_thread_number = [OPT_thread_number; repmat(i,number_of_cores_per_OPT_thread(i,1),1)];
                OPT_worker_number = [OPT_worker_number; [1:number_of_cores_per_OPT_thread(i,1)]'];
            end
            max_number_of_gridcells = min(size(number_of_cores_per_OPT_thread,1), size(run_info.SPATIAL.STATVAR.key,1));
            number_of_runs = size(run_info.SPATIAL.STATVAR.key,1) ./ max_number_of_gridcells;
            run_raster2 = [0; round([number_of_runs:number_of_runs:size(run_info.SPATIAL.STATVAR.key,1)]')];
            run_raster2 = [run_raster2(1:end-1,1)+1 run_raster2(2:end,1)];
            run_raster=[];
            for i=1:size(OPT_thread_number,1)
                run_raster = [run_raster; run_raster2(OPT_thread_number(i,1),:)];
            end
            run_info.TEMP.OPT_thread_number = OPT_thread_number;
            run_info.TEMP.OPT_worker_number = OPT_worker_number;
            run_info.TEMP.run_raster = run_raster;            
            run_info.TEMP.worker_number = worker_number;

            if OPT_worker_number(worker_number) <= size(run_raster,1)./run_info.PARA.number_of_OPT_threads
                for run_number = run_raster(worker_number,1):run_raster(worker_number,2)
                    disp(['running grid cell ' num2str(run_number)])
                    %as normal 1D run

                    %1. run the spatial ACTIONS, this can change the tile->
                    %ensemble paramters
                    for ai=1:size(run_info.SPATIAL.ACTION,1)
                        run_info.SPATIAL.ACTION{ai,1} = assign_tile_properties(run_info.SPATIAL.ACTION{ai,1}, run_number); %writes the provider class, including ENSEMBLE classes
                    end

                    %start the loop
                    % 2. do the ensemble, followed by ACTION
                    if ~isempty(run_info.PARA.ensemble_class) && ~(sum(isnan(run_info.PARA.ensemble_class))>0)
                        disp('get ensemble data')
                        ensemble_class = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.ensemble_class){run_info.PARA.ensemble_class_index,1});
                        run_info.ENSEMBLE = ensemble_class;
                        run_info.ENSEMBLE.RUN_INFO = run_info;
                        run_info.ENSEMBLE = finalize_init(run_info.ENSEMBLE);
                    end
                    %initialization of run_info.OPT needs to go here, with all
                    %fields stored initialized as empty (i.e. the actual values
                    %of the iteration must be assigned later in the loop over
                    %the realizations)
                    number_of_cores_per_DA = size(find(OPT_thread_number == OPT_thread_number(worker_number)),1); %number of cores in this DA (can be different if several threads are in parallel)

                    run_info.OPT = {};
                    for i=1:size(run_info.PARA.optimization_class,1)
                        optimization_class = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.optimization_class{i,1}){run_info.PARA.optimization_class_index(i,1),1});
                        if number_of_cores_per_DA == 1
                            optimization_class.PARA.run_mode = 'OPT_DA_IO_TILE_sequential';
                        else
                            optimization_class.PARA.run_mode = 'OPT_DA_IO_TILE';
                        end
                        optimization_class = finalize_init(optimization_class, run_info); %-> needs to go before realization, so that it can store
                        run_info.OPT = [run_info.OPT; {optimization_class}];
                    end
                    %here the while-loop over time needs to go, ensemble either uses "old
                    %values" when learning coefficient is 0, otherwise taking
                    %values from last round into account
                    %assigns values to ENSEMBLE and sets all stored values to
                    %empty


                    number_of_sequential_DA_runs = run_info.ENSEMBLE.TEMP.ensemble_size ./ number_of_cores_per_DA;
                    run_raster_DA = [0; round([number_of_sequential_DA_runs:number_of_sequential_DA_runs:run_info.ENSEMBLE.TEMP.ensemble_size]')];
                    run_raster_DA = [run_raster_DA(1:end-1,1)+1 run_raster_DA(2:end,1)];

                    ACTIVE = 1;
                    first_init = 1;
                    while ACTIVE == 1
                        for realization_number = run_raster_DA(OPT_worker_number(worker_number),1):run_raster_DA(OPT_worker_number(worker_number),2) %1:run_info.ENSEMBLE.TEMP.ensemble_size %this is parallelized in the real run, no loop (but good to allow for a loop, in case of fewer cores)
                                                        
                            for ai = 1:size(run_info.ENSEMBLE.ACTION,1)
                                run_info.ENSEMBLE.ACTION{ai,1} = assign_tile_properties(run_info.ENSEMBLE.ACTION{ai,1}, realization_number); %writes the provider class
                                %needs to also assign parameters everywhere in
                                %the tree of classes for iteration 2 or
                                %continuation REVISE!!
                            end
                            for ii=1:size(run_info.PARA.optimization_class,1)
                                run_info.OPT{ii,1}.TEMP.realization_number = realization_number;
                            end

                            %3. initialize the DA, especially load the observations ans
                            % establish time raster, add OUT classes to tile with obs time
                            % as proper PARA -> add PARA whether it is part of a DA,
                            % seep point 5
                            %witch 2 and 3, maybe?

                            %tile builder class and PARA must be changed for
                            %iteration 2 and continuation in time, must happen
                            %after DA -> this is were recalcualte_stratigraphy
                            %comes in, in that case no change is made and only
                            %PARAs are changed (works only when there is a
                            %single DA step/optimization for a distinct
                            %framework
                            for ii=1:size(run_info.OPT,1)
                                if first_init == 1 && ii==1
                                    [run_info.OPT{ii,1}, new_tile] = get_tile_class(run_info.OPT{ii,1}, run_info);
                                else
                                    if run_info.OPT{ii,1}.TEMP.ACTIVE == 1
                                        [run_info.OPT{ii,1}, new_tile] = get_tile_class(run_info.OPT{ii,1}, run_info);
                                    end
                                end
                            end
                            % new_tile.RUN_INFO = run_info; see seqential
                            % 
                            % new_tile = finalize_init(new_tile);
                            tile = new_tile;
                            run_info.TILE = tile;
                            tile.PARA.worker_number = worker_number;
                            tile.PARA.range = realization_number;

                            %check if somethig is needed to update obsverations
                            %and the aso reinitialize the Observable OUT
                            %classes if new_interation == 0

                            for ai = 1:size(run_info.ENSEMBLE.ACTION,1)
                                run_info.ENSEMBLE.ACTION{ai,1} = assign_tile_properties(run_info.ENSEMBLE.ACTION{ai,1}, realization_number); %writes the provider class
                                %needs to also assign parameters everywhere in
                                %the tree of classes for iteration 2 or
                                %continuation REVISE!!
                            end

                            for ii=1:size(run_info.OPT,1)
                                if run_info.OPT{ii,1}.TEMP.ACTIVE == 1
                                    run_info.OPT{ii,1} = get_DA_step_time(run_info.OPT{ii,1}, tile);
                                end
                            end
                            %set run_info.OPT{i,1}.TEMP.optimization_time to
                            %get the end-time
                            run_info.STATVAR.next_optimization_time = Inf;
                            for i=1:size(run_info.PARA.optimization_class,1)
                                run_info.STATVAR.next_optimization_time = min(run_info.STATVAR.next_optimization_time, run_info.OPT{i,1}.DA_STEP_TIME);
                            end
                            tile.PARA.start_time = tile.t;
                            tile.PARA.end_time = min(tile.FORCING.PARA.end_time, run_info.STATVAR.next_optimization_time);

                            for ii=1:size(run_info.OPT,1)
                                if run_info.OPT{ii,1}.TEMP.ACTIVE == 1
                                    [run_info.OPT{ii,1}, tile] = add_save_state_classes(run_info.OPT{ii,1}, tile);
                                    [run_info.OPT{ii,1}, tile] = reset_observable_classes(run_info.OPT{ii,1}, tile);
                                end
                            end

                            tile = run_model(tile);  %time integration

                            %move info from OBS/OUT classes to DA
                            for i=1:size(run_info.PARA.optimization_class,1)
                                run_info.OPT{i,1} = move_obs2opt(run_info.OPT{i,1}, tile);
                            end

                        end

                        for ii=1:size(run_info.OPT,1)
                            run_info.OPT{ii,1}.TEMP.ACTIVE = 0;
                        end
                        first_init = 0;

                        %DA step and resampling
                        for i=1:size(run_info.PARA.optimization_class,1)
                            %check which class is triggered through
                            %tile.PARA.end_time
                            run_info.OPT{i,1}.TEMP.number_of_sequential_runs = run_raster_DA(OPT_worker_number(worker_number),2) - run_raster_DA(OPT_worker_number(worker_number),1) + 1;
                            run_info.OPT{i,1} = DA_step(run_info.OPT{i,1}, run_info);
                        end

                        ACTIVE = 0;
                        for i=1:size(run_info.PARA.optimization_class,1)
                            if run_info.OPT{i,1}.TEMP.ACTIVE
                                ACTIVE = 1;
                            end
                        end


                    end
                end
            end
        end

        function [run_info, tile] = run_model_TILE_sequential(run_info)
            tile = 0;
            for run_number = 1:max(run_info.SPATIAL.STATVAR.key) %make dependent on generic variable name so that also DA over multiple grid cells is covered
                %must be loop over all "independent runs", i.e. runs that are ont connected through DA 
                disp(['running grid cell ' num2str(run_number)])
                %as normal 1D run

                %1. run the spatial ACTIONS, this can change the tile->
                %ensemble paramters
                for ai=1:size(run_info.SPATIAL.ACTION,1)
                    run_info.SPATIAL.ACTION{ai,1} = assign_tile_properties(run_info.SPATIAL.ACTION{ai,1}, run_number); %writes the provider class, including ENSEMBLE classes
                end

                %start the loop
                % 2. do the ensemble, followed by ACTION
                if ~isempty(run_info.PARA.ensemble_class) && ~(sum(isnan(run_info.PARA.ensemble_class))>0)
                    disp('get ensemble data')
                    ensemble_class = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.ensemble_class){run_info.PARA.ensemble_class_index,1});
                    run_info.ENSEMBLE = ensemble_class;
                    run_info.ENSEMBLE.RUN_INFO = run_info;
                    run_info.ENSEMBLE = finalize_init(run_info.ENSEMBLE);
                end
                %initialization of run_info.OPT needs to go here, with all
                %fields stored initialized as empty (i.e. the actual values
                %of the iteration must be assigned later in the loop over
                %the realizations)
                run_info.OPT = {};
                for i=1:size(run_info.PARA.optimization_class,1)
                    optimization_class = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.optimization_class{i,1}){run_info.PARA.optimization_class_index(i,1),1});
                    optimization_class = finalize_init(optimization_class, run_info); %-> needs to go before realization, so that it can store
                    run_info.OPT = [run_info.OPT; {optimization_class}];
                end
                %here the while-loop over time needs to go, ensemble either uses "old
                %values" when learning coefficient is 0, otherwise taking
                %values from last round into account
                %assigns values to ENSEMBLE and sets all stored values to
                %empty
                run_info.TEMP.OPT_worker_number = 1;
                ACTIVE = 1;
                first_init = 1;
                while ACTIVE == 1
                    for realization_number = 1:run_info.ENSEMBLE.TEMP.ensemble_size %this is parallelized in the real run, no loop (but good to allow for a loop, in case of fewer cores)
                        %spmd with the correct number of cores, assigned to
                        %this chunk of the spatial problem
                        for ai = 1:size(run_info.ENSEMBLE.ACTION,1)
                            run_info.ENSEMBLE.ACTION{ai,1} = assign_tile_properties(run_info.ENSEMBLE.ACTION{ai,1}, realization_number); %writes the provider class
                            %assign all new PARA's in PROVIDER
                        end
                        for ii=1:size(run_info.PARA.optimization_class,1)
                            run_info.OPT{ii,1}.TEMP.realization_number = realization_number;
                        end

                        %3. initialize the DA, especially load the observations ans
                        % establish time raster, add OUT classes to tile with obs time
                        % as proper PARA -> add PARA whether it is part of a DA,
                        % seep point 5
                        %witch 2 and 3, maybe?

                        %tile builder class and PARA must be changed for
                        %iteration 2 and continuation in time, must happen
                        %after DA -> this is were recalcualte_stratigraphy
                        %comes in, in that case no change is made and only
                        %PARAs are changed (works only when there is a
                        %single DA step/optimization for a distinct
                        %framework
                        for ii=1:size(run_info.OPT,1)
                            if first_init == 1 && ii==1
                                [run_info.OPT{ii,1}, new_tile] = get_tile_class(run_info.OPT{ii,1}, run_info);
                            else
                                if run_info.OPT{ii,1}.TEMP.ACTIVE == 1
                                    [run_info.OPT{ii,1}, new_tile] = get_tile_class(run_info.OPT{ii,1}, run_info);
                                end
                            end
                        end

                        tile = new_tile;
                        run_info.TILE = tile;
                        tile.PARA.worker_number = 1;
                        tile.PARA.range = realization_number;

                        %check if somethig is needed to update obsverations
                        %and the aso reinitialize the Observable OUT
                        %classes if new_interation == 0

                        for ai = 1:size(run_info.ENSEMBLE.ACTION,1)
                            run_info.ENSEMBLE.ACTION{ai,1} = assign_tile_properties(run_info.ENSEMBLE.ACTION{ai,1}, realization_number); %writes the provider class
                            %assign parameters everywhere in the tree of
                            %classes for iteration 2 or continuation
                        end

                        for ii=1:size(run_info.OPT,1)
                            if run_info.OPT{ii,1}.TEMP.ACTIVE == 1 
                                run_info.OPT{ii,1} = get_DA_step_time(run_info.OPT{ii,1}, tile);
                            end
                        end
                        %set run_info.OPT{i,1}.TEMP.optimization_time to
                        %get the end-time
                        run_info.STATVAR.next_optimization_time = Inf;
                        for i=1:size(run_info.PARA.optimization_class,1)
                            run_info.STATVAR.next_optimization_time = min(run_info.STATVAR.next_optimization_time, run_info.OPT{i,1}.DA_STEP_TIME);
                        end
                        tile.PARA.start_time = tile.t;
                        tile.PARA.end_time = min(tile.FORCING.PARA.end_time, run_info.STATVAR.next_optimization_time);

                        for ii=1:size(run_info.OPT,1)
                            if run_info.OPT{ii,1}.TEMP.ACTIVE == 1
                                [run_info.OPT{ii,1}, tile] = add_save_state_classes(run_info.OPT{ii,1}, tile);
                                [run_info.OPT{ii,1}, tile] = reset_observable_classes(run_info.OPT{ii,1}, tile);
                            end
                        end

                        tile = run_model(tile);  %time integration

                        %move info from OBS/OUT classes to DA 
                        for i=1:size(run_info.PARA.optimization_class,1)
                            run_info.OPT{i,1} = move_obs2opt(run_info.OPT{i,1}, tile);
                        end

                    end

                    for ii=1:size(run_info.OPT,1)
                        run_info.OPT{ii,1}.TEMP.ACTIVE = 0;
                    end
                    first_init = 0;

                    %DA step and resampling
                    for i=1:size(run_info.PARA.optimization_class,1)
                        %check which class is triggered through
                        %tile.PARA.end_time
                        run_info.OPT{i,1} = DA_step(run_info.OPT{i,1}, run_info);
                    end

                    ACTIVE = 0;
                    for i=1:size(run_info.PARA.optimization_class,1)
                        if run_info.OPT{i,1}.TEMP.ACTIVE
                            ACTIVE = 1;
                        end
                    end
                end

                %4. initialize TILE, this sets tile.PARA.end_time to
                %forcing endtime
                %5. set tile.PARA.end_time, and start with store_out, then run
                %during store_out, run a function under OPS that loads new
                %obs data if necessary, then add the time raster to OUT
                %initialize optimization classes
                
                %set optimization time (i.e. the DA_step) which becomes tile.PARA.end_time, must
                %be set in the opt classes


                %when done-> collect observations from TILE.OUT
                %do the DA and add the weights to DA.STATVAR 
                %call ensemble class again, with PARA's shifted, and
                %generate new enseble with iteration +1 (DA class does
                %this)
                %assign the correct parameter line with ENSEMBLE->ACTION, specific
                %class needed for this (assigned by DA)
                %depending on DA setting, initialze TILE new (with chnage dparameters), or usethe written files as start 
                %always make sure that tile.t is belwo the true end time

                %while tile.t < ti
                 
                % for i=1:size(run_info.PARA.tile_class,1)
                %     disp(['running tile number ' num2str(i)])
                % 
                %     %add observable OUT classes to TILE, as well as
                %     %OUT_last_timestamp
                % 
                %     %make finalize-init of OPT and ENSEMBLE, with OPT
                %     %setting the MASK value in ENSEMBLE
                % 
                % 
                %     for ai=1:size(run_info.SPATIAL.ACTION,1)
                %         run_info.SPATIAL.ACTION{ai,1} = assign_tile_properties(run_info.SPATIAL.ACTION{ai,1}, run_number); %writes the provider class
                %     end
                % 
                %     new_tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_class{i,1}){run_info.PARA.tile_class_index(i,1),1});
                %     new_tile.RUN_INFO = run_info;
                %     new_tile = finalize_init(new_tile);
                %     tile = new_tile;
                %     run_info.TILE = tile;
                % 
                %     tile.PARA.worker_number = 1;
                %     tile.PARA.range = run_number;
                % 
                %     disp([tile.PARA.latitude tile.PARA.longitude])
                % 
                %     tile = run_model(tile);  %time integration
                % end
            end
        end

        function [run_info, tile] = run_model_MULTITILE_parallel(run_info)
            tile = 0;
            worker_number = spmdIndex();
            max_number_of_gridcells = min(run_info.PARA.number_of_cores, size(run_info.SPATIAL.STATVAR.key,1));
            number_of_runs = size(run_info.SPATIAL.STATVAR.key,1) ./ max_number_of_gridcells;
            run_raster =  [0; round([number_of_runs:number_of_runs:size(run_info.SPATIAL.STATVAR.key,1)]')];
            run_raster = [run_raster(1:end-1,1)+1 run_raster(2:end,1)];
            if worker_number <= size(run_raster,1)
                number_of_grid_cells =  run_raster(worker_number,2)-run_raster(worker_number,1)+1;
                max_number_of_gridcells = min(run_info.PARA.max_number_of_gridcells, number_of_grid_cells);
                number_of_runs = ceil(number_of_grid_cells ./ max_number_of_gridcells);
                run_increment = number_of_grid_cells ./number_of_runs;
                run_raster2 =  [0; round([run_increment:run_increment:number_of_grid_cells]')];
                run_raster2 = [run_raster(worker_number,1)+run_raster2(1:end-1,1) run_raster(worker_number,1)-1+run_raster2(2:end,1)];
                for run_index = 1:size(run_raster2,1)
                    range = [run_raster2(run_index,1):run_raster2(run_index,2)]';
                    disp(['running range index ' num2str(run_index)])

                    for i=1:size(run_info.PARA.tile_class,1)
                        disp(['running tile number ' num2str(i)])
                        for j=1:run_info.PARA.number_of_runs_per_tile(i,1)
                            disp(['running round ' num2str(j)])

                            for ai=1:size(run_info.SPATIAL.ACTION,1)
                                run_info.SPATIAL.ACTION{ai,1} = assign_tile_properties(run_info.SPATIAL.ACTION{ai,1}, range); %writes the provider class
                            end

                            new_tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_class{i,1}){run_info.PARA.tile_class_index(i,1),1});
                            new_tile.RUN_INFO = run_info;
                            new_tile = finalize_init(new_tile);
                            tile = new_tile;
                            run_info.TILE = tile;

                            tile.PARA.worker_number = worker_number;
                            tile.PARA.range = range;

                            tile = run_model(tile);  %time integration
                        end
                    end
                end
            end
        end

        function  [run_info, tile] = run_model_MULTITILE_sequential(run_info)
            tile = 0;
            max_number_of_gridcells = min(run_info.PARA.max_number_of_gridcells, size(run_info.SPATIAL.STATVAR.key,1));
            number_of_runs = ceil(size(run_info.SPATIAL.STATVAR.key,1) ./ max_number_of_gridcells);
            run_increment = size(run_info.SPATIAL.STATVAR.key,1) ./number_of_runs;
            run_raster =  [0; round([run_increment:run_increment:size(run_info.SPATIAL.STATVAR.key,1)]')];
            run_raster = [run_raster(1:end-1,1)+1 run_raster(2:end,1)];

            for run_index = 1:size(run_raster,1)
                disp(['running range index ' num2str(run_index)])
                %as normal 1D run
                for i=1:size(run_info.PARA.tile_class,1)
                    disp(['running tile number ' num2str(i)])
                    for j=1:run_info.PARA.number_of_runs_per_tile(i,1)
                        disp(['running round ' num2str(j)])

                        range = [run_raster(run_index,1):run_raster(run_index,2)]';

                        for ai=1:size(run_info.SPATIAL.ACTION,1)
                            run_info.SPATIAL.ACTION{ai,1} = assign_tile_properties(run_info.SPATIAL.ACTION{ai,1}, range); %writes the provider class
                        end

                        new_tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_class{i,1}){run_info.PARA.tile_class_index(i,1),1});
                        new_tile.RUN_INFO = run_info;
                        new_tile = finalize_init(new_tile);
                        tile = new_tile;
                        run_info.TILE = tile;

                        tile.PARA.worker_number = 1;
                        tile.PARA.range = range;

                        tile = run_model(tile);  %time integration
                    end
                end
            end
        end

        % 
        % function [run_info, tile] = run_model_TILE(run_info)
        % 
        %     tile = 0;
        %     if run_info.PARA.number_of_cores > 1
        %         if isempty(gcp('nocreate'))
        %             parpool(run_info.PARA.number_of_cores)
        %         end
        %         spmd
        %             worker_number = labindex();
        %             max_number_of_gridcells = min(run_info.PARA.number_of_cores, size(run_info.SPATIAL.STATVAR.key,1));
        %             number_of_runs = size(run_info.SPATIAL.STATVAR.key,1) ./ max_number_of_gridcells;
        %             run_raster =  [0; round([number_of_runs:number_of_runs:size(run_info.SPATIAL.STATVAR.key,1)]')];
        %             run_raster = [run_raster(1:end-1,1)+1 run_raster(2:end,1)];
        %             if worker_number <= size(run_raster,1)
        %                 for run_number = run_raster(worker_number,1):run_raster(worker_number,2)
        % 
        %                     disp(['running grid cell ' num2str(run_number)])
        %                     %as normal 1D run
        %                     for i=1:size(run_info.PARA.tile_class,1)
        %                         disp(['running tile number ' num2str(i)])
        %                         for j=1:run_info.PARA.number_of_runs_per_tile(i,1)
        %                             disp(['running round ' num2str(j)])
        % 
        %                             for ai=1:size(run_info.SPATIAL.ACTION,1)
        %                                 run_info.SPATIAL.ACTION{ai,1} = assign_tile_properties(run_info.SPATIAL.ACTION{ai,1}, run_number); %writes the provider class
        %                             end
        % 
        %                             new_tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_class{i,1}){run_info.PARA.tile_class_index(i,1),1});
        %                             new_tile.RUN_INFO = run_info;
        %                             new_tile = finalize_init(new_tile);
        %                             tile = new_tile;
        %                             run_info.TILE = tile;
        % 
        %                             tile.PARA.worker_number = worker_number;
        %                             tile.PARA.range = run_number;
        % 
        %                             tile = run_model(tile);  %time integration
        %                         end
        %                     end
        %                 end
        %             end
        %             if run_info.PARA.close_parallel_pool
        %                 delete(gcp('nocreate'));
        %             end
        %         end
        %     else
        %         for run_number = 1:size(run_info.SPATIAL.STATVAR.key,1)    
        % 
        %             disp(['running grid cell ' num2str(run_number)])
        %             %as normal 1D run
        %             for i=1:size(run_info.PARA.tile_class,1) 
        %                 disp(['running tile number ' num2str(i)])
        %                 for j=1:run_info.PARA.number_of_runs_per_tile(i,1)
        %                     disp(['running round ' num2str(j)])
        % 
        %                     for ai=1:size(run_info.SPATIAL.ACTION,1)
        %                         run_info.SPATIAL.ACTION{ai,1} = assign_tile_properties(run_info.SPATIAL.ACTION{ai,1}, run_number); %writes the provider class
        %                     end
        % 
        %                     new_tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_class{i,1}){run_info.PARA.tile_class_index(i,1),1});
        %                     new_tile.RUN_INFO = run_info;
        %                     new_tile = finalize_init(new_tile);
        %                     tile = new_tile;
        %                     run_info.TILE = tile;
        % 
        %                     tile.PARA.worker_number = 1;
        %                     tile.PARA.range = run_number;
        % 
        %                     [tile.PARA.latitude tile.PARA.longitude]
        % 
        %                     tile = run_model(tile);  %time integration
        %                 end
        %             end
        %         end
        %     end
        % 
        % end
        % 
        % 
        % function [run_info, tile] = run_model_MULTITILE(run_info) %not tested yet
        % 
        %     tile = 0;
        %     if run_info.PARA.number_of_cores > 1
        %         if isempty(gcp('nocreate'))
        %             parpool(run_info.PARA.number_of_cores)
        %         end
        %         spmd
        %             worker_number = labindex();
        %             max_number_of_gridcells = min(run_info.PARA.number_of_cores, size(run_info.SPATIAL.STATVAR.key,1));
        %             number_of_runs = size(run_info.SPATIAL.STATVAR.key,1) ./ max_number_of_gridcells;
        %             run_raster =  [0; round([number_of_runs:number_of_runs:size(run_info.SPATIAL.STATVAR.key,1)]')];
        %             run_raster = [run_raster(1:end-1,1)+1 run_raster(2:end,1)];
        %             if worker_number <= size(run_raster,1)
        %                 number_of_grid_cells =  run_raster(worker_number,2)-run_raster(worker_number,1)+1;
        %                 max_number_of_gridcells = min(run_info.PARA.max_number_of_gridcells, number_of_grid_cells);
        %                 number_of_runs = ceil(number_of_grid_cells ./ max_number_of_gridcells);
        %                 run_increment = number_of_grid_cells ./number_of_runs;
        %                 run_raster2 =  [0; round([run_increment:run_increment:number_of_grid_cells]')];
        %                 run_raster2 = [run_raster(worker_number,1)+run_raster2(1:end-1,1) run_raster(worker_number,1)-1+run_raster2(2:end,1)];
        %                 for run_index = 1:size(run_raster2,1)
        %                     range = [run_raster2(run_index,1):run_raster2(run_index,2)]';
        %                     disp(['running range index ' num2str(run_index)])
        % 
        %                     for i=1:size(run_info.PARA.tile_class,1)
        %                         disp(['running tile number ' num2str(i)])
        %                         for j=1:run_info.PARA.number_of_runs_per_tile(i,1)
        %                             disp(['running round ' num2str(j)])
        % 
        %                             for ai=1:size(run_info.SPATIAL.ACTION,1)
        %                                 run_info.SPATIAL.ACTION{ai,1} = assign_tile_properties(run_info.SPATIAL.ACTION{ai,1}, range); %writes the provider class
        %                             end
        % 
        %                             new_tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_class{i,1}){run_info.PARA.tile_class_index(i,1),1});
        %                             new_tile.RUN_INFO = run_info;
        %                             new_tile = finalize_init(new_tile);
        %                             tile = new_tile;
        %                             run_info.TILE = tile;
        % 
        %                             tile.PARA.worker_number = worker_number;
        %                             tile.PARA.range = range;
        % 
        %                             tile = run_model(tile);  %time integration
        %                         end
        %                     end
        %                 end
        %             end
        %             if run_info.PARA.close_parallel_pool
        %                 delete(gcp('nocreate'));
        %             end
        %         end
        %     else
        % 
        %         max_number_of_gridcells = min(run_info.PARA.max_number_of_gridcells, size(run_info.SPATIAL.STATVAR.key,1));
        %         number_of_runs = ceil(size(run_info.SPATIAL.STATVAR.key,1) ./ max_number_of_gridcells);
        %         run_increment = size(run_info.SPATIAL.STATVAR.key,1) ./number_of_runs;
        %         run_raster =  [0; round([run_increment:run_increment:size(run_info.SPATIAL.STATVAR.key,1)]')];
        %         run_raster = [run_raster(1:end-1,1)+1 run_raster(2:end,1)];
        % 
        %         for run_index = 1:size(run_raster,1)
        %             disp(['running range index ' num2str(run_index)])
        %             %as normal 1D run
        %             for i=1:size(run_info.PARA.tile_class,1) 
        %                 disp(['running tile number ' num2str(i)])
        %                 for j=1:run_info.PARA.number_of_runs_per_tile(i,1)
        %                     disp(['running round ' num2str(j)])
        % 
        %                     range = [run_raster(run_index,1):run_raster(run_index,2)]';
        % 
        %                     for ai=1:size(run_info.SPATIAL.ACTION,1)
        %                         run_info.SPATIAL.ACTION{ai,1} = assign_tile_properties(run_info.SPATIAL.ACTION{ai,1}, range); %writes the provider class
        %                     end
        % 
        %                     new_tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_class{i,1}){run_info.PARA.tile_class_index(i,1),1});
        %                     new_tile.RUN_INFO = run_info;
        %                     new_tile = finalize_init(new_tile);
        %                     tile = new_tile;
        %                     run_info.TILE = tile;
        % 
        %                     tile.PARA.worker_number = 1;
        %                     tile.PARA.range = range;
        % 
        %                     tile = run_model(tile);  %time integration
        %                 end
        %             end
        %         end
        %     end
        % 
        % end


        
        %-------------param file generation-----
        function run_info = param_file_info(run_info)
            run_info = provide_PARA(run_info);

            run_info.PARA.STATVAR = [];
            run_info.PARA.class_category = 'RUN_INFO';
            run_info.PARA.default_value = [];
            run_info.PARA.comment = [];
            
            run_info.PARA.comment.number_of_cores = {'number of cores to be used for calculation'};
            run_info.PARA.default_value.number_of_cores = {2};
            
            run_info.PARA.options.tile_class.name =  'H_LIST';
            run_info.PARA.options.tile_class.entries_x = {'TILE_1D_standard' 'TILE_1D_standard'};
            
            run_info.PARA.options.tile_class_index.name =  'H_LIST'; 
            run_info.PARA.options.tile_class_index.entries_x = {1 2};
            
            run_info.PARA.options.number_of_runs_per_tile.name =  'H_LIST'; % 
            run_info.PARA.options.number_of_runs_per_tile.entries_x = {1 1};
            
            run_info.PARA.comment.ensemble_class = {'projection class providing providing information on the locations and additinal data for each target point'};
            
        end
        
        
            
%             number_of_tiles = ceil(run_info.PARA.total_number_of_cells ./ run_info.PARA.number_of_cells_per_tile);
%             domains_per_worker = max(1, floor(number_of_tiles ./ run_info.PARA.num_ranks - 1e-12));
%             
%             while ~terminate_check(run_info, 'section2_done', number_of_tiles)
%                 %for run_index = 1:number_of_tiles
%                 for run_index=[[(run_info.PARA.my_rank-1).*domains_per_worker+1:number_of_tiles] [1:(run_info.PARA.my_rank-1).*domains_per_worker]]
%                     
%                     if run_slice_yes_no(run_info, 'section2_started', 'section2_done', run_index)
%                     
%                         crap = write_check(run_info, 'section2_started', run_index);
%                         
%                         disp(['running range index ' num2str(run_index)])
%                         
%                         range = [(run_index-1).*run_info.PARA.number_of_cells_per_tile+1:min(run_index.*run_info.PARA.number_of_cells_per_tile, run_info.PARA.total_number_of_cells)]';
%                         
%                         for i=1:size(run_info.PARA.tile_class,1)
%                             disp(['running tile number ' num2str(i)])
%                             for j=1:run_info.PARA.number_of_runs(i,1)
%                                 disp(['running round ' num2str(j)])
%                                 
%                                 %load the next tile from the PROVIDER
%                                 tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_class{i,1}){run_info.PARA.tile_class_index(i,1),1});
%                                 tile.PARA.number_of_realizations = size(range,1);
%                                 tile.PARA.range = range;
%                                 
%                                 tile.PARA.geothermal = run_info.STATVAR.geothermal(range,1);
%                                 
%                                 
%                                 tile.RUN_INFO = run_info;
%                                 
%                                 tile = finalize_init(tile); %here, tile can still access a potentially existing tile through til.RUN_INFO.TILE
%                                 
%                                 tile = run_model(tile);
%                             end
%                         end
%                         crap = write_check(run_info, 'section2_done', run_index);
%                     end
%                 end
%             end
            

%         end
        
        
%         function [run_info, tile] = run_postproc(run_info)
%             
%             number_of_tiles = ceil(run_info.PARA.total_number_of_cells ./ run_info.PARA.number_of_cells_per_tile);
%             %domains_per_worker = max(1, floor(number_of_tiles ./ run_info.PARA.num_ranks - 1e-12))
%             %domains_per_worker = max(1, ceil(number_of_tiles ./ run_info.PARA.num_ranks));
%             domains_per_worker_breaks = round([0:number_of_tiles./run_info.PARA.num_ranks:number_of_tiles]');
%             
%             for i=1:size(run_info.PARA.tile_postproc_class,1)
%                 tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_postproc_class{i,1}){run_info.PARA.tile_postproc_class_index(i,1),1});
%                 
%                 %tile.PARA.run_index = [(run_info.PARA.my_rank-1).*domains_per_worker+1:min(run_info.PARA.my_rank.*domains_per_worker, number_of_tiles)]';
%                 tile.PARA.run_index = [domains_per_worker_breaks(run_info.PARA.my_rank)+1:domains_per_worker_breaks(run_info.PARA.my_rank + 1)]';
%                 if ~isempty(tile.PARA.run_index)
%                     
%                     tile.RUN_INFO = run_info;
%                     
%                     tile = finalize_init(tile); %here, tile can still access a potentially existing tile through til.RUN_INFO.TILE
%                     
%                     tile = run_model(tile);
%                 end
%             end
%        
%             %terminate parallel environment 
%             if run_info.PARA.parallelized == 1
%                 NMPI_Finalize(); % End the MPI communication
%             end
%         end
 
        
%         
%         function [run_info, tile] = run_model(run_info)
%             
%             number_of_tiles = ceil(run_info.PARA.total_number_of_cells ./ run_info.PARA.number_of_cells_per_tile);
%             domains_per_worker = max(1, floor(number_of_tiles ./ run_info.PARA.num_ranks - 1e-12));
%             
%             %parallelize this
%             for run_index = 1:number_of_tiles
%                 
%                 disp(['running range index ' num2str(run_index)])
%                 
%                 range = [(run_index-1).*run_info.PARA.number_of_cells_per_tile+1:min(run_index.*run_info.PARA.number_of_cells_per_tile, run_info.PARA.total_number_of_cells)]';
%                 
%                 for i=1:size(run_info.PARA.tile_class,1)
%                     disp(['running tile number ' num2str(i)])
%                     for j=1:run_info.PARA.number_of_runs(i,1)
%                         disp(['running round ' num2str(j)])
%                         
%                         %load the next tile from the PROVIDER
%                         tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_class{i,1}){run_info.PARA.tile_class_index(i,1),1});
%                         tile.PARA.number_of_realizations = size(range,1);
%                         tile.PARA.range = range;
%                         
%                         tile.PARA.geothermal = run_info.STATVAR.geothermal(range,1);
%                         
%                         
%                         %REMOVE
% %                         [~, pos] = max(run_info.STATVAR.landcover(range,1),[], 2);
% %                         tile.PARA.stratigraphy = pos(1,1) .*0 +1;
%                         %REMOVE
%                         
%                         %tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_class){run_info.PARA.tile_class_index,1});
%                         tile.RUN_INFO = run_info;
%                         
%                         tile = finalize_init(tile); %here, tile can still access a potentially existing tile through til.RUN_INFO.TILE
%                        
%                         tile = run_model(tile);
%                     end
%                 end
%             end
%             
%             if run_info.PARA.parallelized == 1
%                 NMPI_Finalize(); % End the MPI communication
%             end
%         end
 
        
        
        %non-mandatory
        
%         function run_info = customize(run_info)
% 
%             
%         end
%         
%         
%          function run_info = parallelize_preproc(run_info)
%             number_of_years = [];
%             load_index = [1; 1; 1.8];
%             number_of_cores = run_info.PARA.num_ranks;
%             
%             my_core = run_info.PARA.my_rank;
%             
%             number_of_slices = 46;
%             start_year= [];
% 
%             run_info.PARA.tile_active = [0;0;0];
%             
%             number_of_slices = 46;
%             start_year= [];
%             
%             for i=1:size(run_info.PARA.tile_preproc_class_index,1)
%                 tile = run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_preproc_class{i,1}){run_info.PARA.tile_preproc_class_index(i,1),1};
%                 forcing = run_info.PPROVIDER.CLASSES.(tile.PARA.forcing_class){tile.PARA.forcing_class_index,1};
%                 number_of_years = [number_of_years; forcing.PARA.end_time(1)-forcing.PARA.start_time(1)+1];
%                 start_year= [start_year; forcing.PARA.start_time(1)];
%             end
%             
%             load_per_core = sum(number_of_years .* number_of_slices .* load_index) ./ number_of_cores;
%             number_of_slices_per_tile = number_of_years .* number_of_slices;
%             load_per_tile = number_of_years .* number_of_slices .* load_index;
%             load_per_tile_acc = cumsum(load_per_tile);
%             
% 
%             load_start = load_per_core .* (my_core-1)+1;
%             load_end = load_per_core .* my_core;
%             
%             tile_num = 1;
%             slice_count = 0;
%             while tile_num<=3
%                 if load_start > load_per_tile_acc(tile_num)
%                     tile_num=tile_num+1;
%                 else
%                     break
%                 end
%             end
%             tile_num_start = tile_num;
%             
%             tile_num = 1;
%             slice_count = 0;
%             while tile_num<=3
%                 if load_end > load_per_tile_acc(tile_num)
%                     tile_num=tile_num+1;
%                 else
%                     break
%                 end
%             end
%             tile_num_end = tile_num;
%             run_info.PARA.tile_active(tile_num_start:tile_num_end) = run_info.PARA.tile_active(tile_num_start:tile_num_end)+1;
%             
%             load_per_tile_acc = [0; load_per_tile_acc];
%             slice_count_start = round((load_start - load_per_tile_acc(tile_num_start)) ./ load_index(tile_num_start));
%             if slice_count_start > 1
%                 %new_start_time = datenum(start_year(tile_num_start,1) + floor(slice_count_start./46),1,1) + (mod(slice_count_start, 46)-1).*8;  
%                 new_start_time = datenum(start_year(tile_num_start,1) + floor((slice_count_start-1)./46),1,1) + (mod(slice_count_start-1, 46)).*8;
%                 tile = run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_preproc_class{tile_num_start,1}){run_info.PARA.tile_preproc_class_index(tile_num_start,1),1};
%                 [year,month,day,~,~,~] = datevec(new_start_time);
%                 run_info.PPROVIDER.CLASSES.(tile.PARA.forcing_class){tile.PARA.forcing_class_index,1}.PARA.start_time = [year; month; day];
%                 disp([num2str(run_info.PARA.my_rank) ' ' num2str(year) ' ' num2str(month) ' ' num2str(day)])
%             end
%             slice_count_end =  round((load_end - load_per_tile_acc(tile_num_end)) ./ load_index(tile_num_end))+2; %add the 2 for overlap
%             if slice_count_end < number_of_slices_per_tile(tile_num_end)
%                 %new_end_time = datenum(start_year(tile_num_end,1) + floor(slice_count_end./46),1,1) + (mod(slice_count_end, 46)-1).*8;
%                 new_end_time = datenum(start_year(tile_num_end,1) + floor((slice_count_end-1)./46),1,1) + (mod(slice_count_end-1, 46)).*8;
%                 tile = run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_preproc_class{tile_num_end,1}){run_info.PARA.tile_preproc_class_index(tile_num_end,1),1};
%                 [year,month,day,~,~,~] = datevec(new_end_time);
%                 run_info.PPROVIDER.CLASSES.(tile.PARA.forcing_class){tile.PARA.forcing_class_index,1}.PARA.end_time = [year; month; day];
%                 disp([num2str(run_info.PARA.my_rank) ' ' num2str(year) ' ' num2str(month) ' ' num2str(day)])
%             end
% 
%          end
%         
%         
%         
%         function run_info = parallelize_preproc_old(run_info)
%             number_of_years = [];
%             load_index = [1; 1; 1.8];
%             number_of_cores = run_info.PARA.num_ranks;
%             
%             my_core = run_info.PARA.my_rank;
%             
%             number_of_slices = 46;
%             start_year= [];
%             
% 
%             run_info.PARA.tile_active = [0;0;0];
%             
%             number_of_slices = 46;
%             start_year= [];
%             
%             for i=1:size(run_info.PARA.tile_preproc_class_index,1)
%                 tile = run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_preproc_class{i,1}){run_info.PARA.tile_preproc_class_index(i,1),1};
%                 forcing = run_info.PPROVIDER.CLASSES.(tile.PARA.forcing_class){tile.PARA.forcing_class_index,1};
%                 number_of_years = [number_of_years; forcing.PARA.end_time(1)-forcing.PARA.start_time(1)+1];
%                 start_year= [start_year; forcing.PARA.start_time(1)];
%             end
%             
%             load_per_core = sum(number_of_years .* number_of_slices .* load_index) ./ number_of_cores;
%             number_of_slices_per_tile = number_of_years .* number_of_slices;
%             load_per_tile = number_of_years .* number_of_slices .* load_index;
%             load_per_tile_acc = cumsum(load_per_tile);
%             
%             result =[];
%             tile_number = 1;
%             core = 1;
%             slice_count = 0;
%             end_last_slice = 0;
%             load_count = 0;
%             
%             
%             while core <= number_of_cores %tile_number<=3
%                 load_count = load_count + load_per_core;
%                 if load_count < load_per_tile_acc(tile_number)
%                     slice_count = slice_count + load_per_core./load_index(tile_number);
%                     result=[result; [core tile_number round(slice_count)]];
%                 else
%                     while tile_number<3
%                         remaining_load = load_count - load_per_tile_acc(tile_number);
%                         if  remaining_load >0
%                             result=[result; [core tile_number number_of_slices_per_tile(tile_number)]];
%                         else
%                             break
%                         end
%                         tile_number = tile_number+1;
%                         slice_count =  remaining_load./load_index(tile_number);
%                         if round(slice_count) <= number_of_slices_per_tile(tile_number)
%                             result=[result; [core tile_number round(slice_count)]];
%                         end
%                     end
%                 end
%                 
%                 core = core +1;
%             end
%             
%             result2=[];
%             ti=0;
%             for i=1:size(result,1)
%                 
%                 if ti ~= result(i,2)
%                     result2 = [result2; [start_year(result(i,2)) + floor(result(i,3)./46) mod(result(i,3),46).*8 + 1 ...
%                         datenum(start_year(result(i,2)) + floor(result(i,3)./46),1,1) + mod(result(i,3),46).*8-1  ...
%                         datenum(start_year(result(i,2)),1,1)]];
%                     ti = result(i,2);
%                 else
%                     result2 = [result2; [start_year(result(i,2)) + floor(result(i,3)./46) mod(result(i,3),46).*8 + 1 ...
%                         datenum(start_year(result(i,2)) + floor(result(i,3)./46),1,1) + mod(result(i,3),46).*8-1  ...
%                         datenum(start_year(result(i-1,2)) + floor(result(i-1,3)./46),1,1) + mod(result(i-1,3),46).*8]];
%                 end
%             end
%             % result2 = [result2 [datenum(start_year(1,1), 1, 1); result2(1:end-1,3)]];
%             % result2(:,3) = result2(:,3)-1;
%             
%             for i=1:size(result,1)
%                 if my_core == result(i,1)
%                     tile_num =result(i,2);
%                     run_info.PARA.tile_active(tile_num,1) = 1;
% 
%                     tile = run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_preproc_class{tile_num,1}){run_info.PARA.tile_preproc_class_index(tile_num,1),1};
%                     [year,month,day,~,~,~] = datevec(result2(i,4));
%                     run_info.PPROVIDER.CLASSES.(tile.PARA.forcing_class){tile.PARA.forcing_class_index,1}.PARA.start_time = [year; month; day]; 
%                     disp([num2str(run_info.PARA.my_rank) num2str(year) ' ' num2str(month) ' ' num2str(day)])
%                     [year,month,day,~,~,~] = datevec(result2(i,3));
%                     run_info.PPROVIDER.CLASSES.(tile.PARA.forcing_class){tile.PARA.forcing_class_index,1}.PARA.end_time = [year; month; day]; 
%                     disp([num2str(run_info.PARA.my_rank) num2str(year) ' ' num2str(month) ' ' num2str(day)])
% %                     new_start_time = datestr(result2(i,4))
% %                     new_end_time = datestr(result2(i,3))
%                 end
%                 
%             end
%             
%         end
        
        
    end
end



