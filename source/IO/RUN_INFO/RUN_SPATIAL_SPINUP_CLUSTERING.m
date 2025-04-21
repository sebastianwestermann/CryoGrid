%========================================================================
% CryoGrid RUN_INFO class RUN_SPATIAL_SPINUP
% as RUN_SPATIAL, but providing the possibility to apply clustering to
% reduce the number of points, using an appropriate CLUSERING class
%
% S. westermann, Dec 2022
%========================================================================

classdef RUN_SPATIAL_SPINUP_CLUSTERING < matlab.mixin.Copyable
    
    properties
        PPROVIDER
        PARA
        CONST
        STATVAR
        SPATIAL
        CLUSTER
        TILE
    end
    
    
    methods
        
        
        function run_info = provide_PARA(run_info)
            
            run_info.PARA.run_mode = 'TILE';
            run_info.PARA.max_number_of_gridcells = Inf;
            run_info.PARA.parallel_pool_open = 0; %1 if pool is already open

            run_info.PARA.number_of_cores = [];
                        
            run_info.PARA.tile_class = [];
            run_info.PARA.tile_class_index = [];
            run_info.PARA.number_of_runs_per_tile = []; %vector
            
            run_info.PARA.projection_class = [];
            run_info.PARA.projection_class_index = [];
            
            run_info.PARA.clustering_class = [];
            run_info.PARA.clustering_class_index = [];
        end
        
        function run_info = provide_CONST(run_info)
            
        end
        
        function run_info = provide_STATVAR(run_info)
            
        end
        
        
        function run_info = finalize_init(run_info)
            
            disp('get spatial data')
            %provided by coordinate system MODIS LST classes
            if ~isempty(run_info.PARA.projection_class) && ~(sum(isnan(run_info.PARA.projection_class))>0)
                spatial_class = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.projection_class){run_info.PARA.projection_class_index,1});
                if ~spatial_class.PARA.new_reference
                    spatial_class.STATVAR = run_info.SPATIAL.STATVAR;
                end
                run_info.SPATIAL = spatial_class;
                run_info.SPATIAL.RUN_INFO = run_info;
                run_info.SPATIAL = finalize_init(run_info.SPATIAL);
            end
            
            if ~isempty(run_info.PARA.clustering_class) && ~(sum(isnan(run_info.PARA.clustering_class))>0)
                run_info.CLUSTER = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.clustering_class){run_info.PARA.clustering_class_index,1});
                run_info.CLUSTER.RUN_INFO = run_info;
                run_info.CLUSTER = finalize_init(run_info.CLUSTER);
                run_info.CLUSTER = compute_clusters(run_info.CLUSTER);
            end
            if ~run_info.PARA.parallel_pool_open
                save([run_info.PPROVIDER.PARA.result_path run_info.PPROVIDER.PARA.run_name '/run_parameters.mat'], 'run_info', '-v7.3')
            else
                worker_number = spmdIndex;
                if worker_number==1
                    save([run_info.PPROVIDER.PARA.result_path run_info.PPROVIDER.PARA.run_name '/run_parameters.mat'], 'run_info', '-v7.3')
                end
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
            worker_number = spmdIndex;
            max_number_of_gridcells = min(run_info.PARA.number_of_cores, size(run_info.CLUSTER.STATVAR.sample_centroid_index,1));
            number_of_runs = size(run_info.CLUSTER.STATVAR.sample_centroid_index,1) ./ max_number_of_gridcells;
            run_raster =  [0; round([number_of_runs:number_of_runs:size(run_info.CLUSTER.STATVAR.sample_centroid_index,1)]')];
            run_raster = [run_raster(1:end-1,1)+1 run_raster(2:end,1)];
            if worker_number <= size(run_raster,1)
                for sample_number = run_raster(worker_number,1):run_raster(worker_number,2)

                    run_number = run_info.CLUSTER.STATVAR.sample_centroid_index(sample_number,1);

                    disp(['running grid cell ' num2str(run_number)])
                    %as normal 1D run
                    for i=1:size(run_info.PARA.tile_class,1) %can be parallelized
                        disp(['running tile number ' num2str(i)])

                        for j=1:run_info.PARA.number_of_runs_per_tile(i,1)
                            disp(['running round ' num2str(j)])

                            for ai=1:size(run_info.SPATIAL.ACTION,1)
                                run_info.SPATIAL.ACTION{ai,1} = assign_tile_properties(run_info.SPATIAL.ACTION{ai,1}, run_number); %writes the provider class
                            end

                            new_tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_class{i,1}){run_info.PARA.tile_class_index(i,1),1});
                            new_tile.RUN_INFO = run_info;
                            new_tile = finalize_init(new_tile);
                            tile = new_tile;
                            run_info.TILE = tile;

                            tile.PARA.worker_number = worker_number;
                            tile.PARA.range = run_number;

                            tile = run_model(tile);  %time integration
                        end
                    end
                end
            end
        end

        function [run_info, tile] = run_model_TILE_sequential(run_info)
            tile = 0;
            for sample_number = 1:size(run_info.CLUSTER.STATVAR.sample_centroid_index,1)

                run_number = run_info.CLUSTER.STATVAR.sample_centroid_index(sample_number,1);

                disp(['running grid cell ' num2str(run_number)])
                %as normal 1D run
                for i=1:size(run_info.PARA.tile_class,1) %can be parallelized
                    disp(['running tile number ' num2str(i)])

                    for j=1:run_info.PARA.number_of_runs_per_tile(i,1)
                        disp(['running round ' num2str(j)])

                        for ai=1:size(run_info.SPATIAL.ACTION,1)
                            run_info.SPATIAL.ACTION{ai,1} = assign_tile_properties(run_info.SPATIAL.ACTION{ai,1}, run_number); %writes the provider class
                        end
                        new_tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_class{i,1}){run_info.PARA.tile_class_index(i,1),1});
                        new_tile.RUN_INFO = run_info;
                        new_tile = finalize_init(new_tile);
                        tile = new_tile;
                        run_info.TILE = tile;

                        tile.PARA.worker_number = 1;
                        tile.PARA.range = run_number;

                        tile = run_model(tile);  %time integration
                    end
                end
            end
        end

        function [run_info, tile] = run_model_MULTITILE_parallel(run_info)
            tile = 0;
            worker_number = spmdIndex;
            max_number_of_gridcells = min(run_info.PARA.number_of_cores, size(run_info.CLUSTER.STATVAR.sample_centroid_index,1));
            number_of_runs = size(run_info.CLUSTER.STATVAR.sample_centroid_index,1) ./ max_number_of_gridcells;
            run_raster =  [0; round([number_of_runs:number_of_runs:size(run_info.CLUSTER.STATVAR.sample_centroid_index,1)]')];
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

                    range = run_info.CLUSTER.STATVAR.sample_centroid_index(range,1);

                    %as normal 1D run
                    for i=1:size(run_info.PARA.tile_class,1) %can be parallelized
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

        function [run_info, tile] = run_model_MULTITILE_sequential(run_info)
            tile = 0;
            max_number_of_gridcells = min(run_info.PARA.max_number_of_gridcells, size(run_info.CLUSTER.STATVAR.sample_centroid_index,1));
            number_of_runs = ceil(size(run_info.CLUSTER.STATVAR.sample_centroid_index,1) ./ max_number_of_gridcells);
            run_increment = size(run_info.CLUSTER.STATVAR.sample_centroid_index,1) ./number_of_runs;
            run_raster =  [0; round([run_increment:run_increment:size(run_info.SPATIAL.STATVAR.key,1)]')];
            run_raster = [run_raster(1:end-1,1)+1 run_raster(2:end,1)];

            for run_index = 1:size(run_raster,1)

                range = [run_raster(run_index,1):run_raster(run_index,2)]';
                range = run_info.CLUSTER.STATVAR.sample_centroid_index(range,1);

                disp(['running range index ' num2str(run_index)])
                %as normal 1D run
                for i=1:size(run_info.PARA.tile_class,1) %can be parallelized
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
        %     if run_info.PARA.number_of_cores > 1 %parallelized
        %         if isempty(gcp('nocreate'))
        %             parpool(run_info.PARA.number_of_cores)
        %         end
        %         spmd
        %             worker_number = labindex;
        %             max_number_of_gridcells = min(run_info.PARA.number_of_cores, size(run_info.CLUSTER.STATVAR.sample_centroid_index,1));
        %             number_of_runs = size(run_info.CLUSTER.STATVAR.sample_centroid_index,1) ./ max_number_of_gridcells;
        %             run_raster =  [0; round([number_of_runs:number_of_runs:size(run_info.CLUSTER.STATVAR.sample_centroid_index,1)]')];
        %             run_raster = [run_raster(1:end-1,1)+1 run_raster(2:end,1)];
        %             if worker_number <= size(run_raster,1)
        %                 for sample_number = run_raster(worker_number,1):run_raster(worker_number,2)
        % 
        %                     run_number = run_info.CLUSTER.STATVAR.sample_centroid_index(sample_number,1);
        % 
        %                     disp(['running grid cell ' num2str(run_number)])
        %                     %as normal 1D run
        %                     for i=1:size(run_info.PARA.tile_class,1) %can be parallelized
        %                         disp(['running tile number ' num2str(i)])
        % 
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
        %         for sample_number = 1:size(run_info.CLUSTER.STATVAR.sample_centroid_index,1)
        % 
        %             run_number = run_info.CLUSTER.STATVAR.sample_centroid_index(sample_number,1);
        % 
        % 
        %             disp(['running grid cell ' num2str(run_number)])
        %             %as normal 1D run
        %             for i=1:size(run_info.PARA.tile_class,1) %can be parallelized
        %                 disp(['running tile number ' num2str(i)])
        % 
        %                 for j=1:run_info.PARA.number_of_runs_per_tile(i,1)
        %                     disp(['running round ' num2str(j)])
        % 
        %                     for ai=1:size(run_info.SPATIAL.ACTION,1)
        %                         run_info.SPATIAL.ACTION{ai,1} = assign_tile_properties(run_info.SPATIAL.ACTION{ai,1}, run_number); %writes the provider class
        %                     end
        %                     new_tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_class{i,1}){run_info.PARA.tile_class_index(i,1),1});
        %                     new_tile.RUN_INFO = run_info;
        %                     new_tile = finalize_init(new_tile);
        %                     tile = new_tile;
        %                     run_info.TILE = tile;
        % 
        %                     tile.PARA.worker_number = 1;
        %                     tile.PARA.range = run_number;
        % 
        %                     tile = run_model(tile);  %time integration
        %                 end
        %             end
        %         end
        %     end 
        % end
        % 
        % 
        % function [run_info, tile] = run_model_MULTITILE(run_info)
        % 
        %     tile = 0;
        %     if run_info.PARA.number_of_cores > 1 %parallelized
        %         if isempty(gcp('nocreate'))
        %             parpool(run_info.PARA.number_of_cores)
        %         end
        %         spmd
        %             worker_number = labindex;
        %             max_number_of_gridcells = min(run_info.PARA.number_of_cores, size(run_info.CLUSTER.STATVAR.sample_centroid_index,1));
        %             number_of_runs = size(run_info.CLUSTER.STATVAR.sample_centroid_index,1) ./ max_number_of_gridcells;
        %             run_raster =  [0; round([number_of_runs:number_of_runs:size(run_info.CLUSTER.STATVAR.sample_centroid_index,1)]')];
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
        %                     range = run_info.CLUSTER.STATVAR.sample_centroid_index(range,1);
        % 
        %                     %as normal 1D run
        %                     for i=1:size(run_info.PARA.tile_class,1) %can be parallelized
        %                         disp(['running tile number ' num2str(i)])
        % 
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
        %         max_number_of_gridcells = min(run_info.PARA.max_number_of_gridcells, size(run_info.CLUSTER.STATVAR.sample_centroid_index,1));
        %         number_of_runs = ceil(size(run_info.CLUSTER.STATVAR.sample_centroid_index,1) ./ max_number_of_gridcells);
        %         run_increment = size(run_info.CLUSTER.STATVAR.sample_centroid_index,1) ./number_of_runs;
        %         run_raster =  [0; round([run_increment:run_increment:size(run_info.SPATIAL.STATVAR.key,1)]')];
        %         run_raster = [run_raster(1:end-1,1)+1 run_raster(2:end,1)];
        % 
        %         for run_index = 1:size(run_raster,1)
        % 
        %             range = [run_raster(run_index,1):run_raster(run_index,2)]';
        %             range = run_info.CLUSTER.STATVAR.sample_centroid_index(range,1);
        % 
        %             disp(['running range index ' num2str(run_index)])
        %             %as normal 1D run
        %             for i=1:size(run_info.PARA.tile_class,1) %can be parallelized
        %                 disp(['running tile number ' num2str(i)])
        % 
        %                 for j=1:run_info.PARA.number_of_runs_per_tile(i,1)
        %                     disp(['running round ' num2str(j)])
        % 
        %                     for ai=1:size(run_info.SPATIAL.ACTION,1)
        %                         run_info.SPATIAL.ACTION{ai,1} = assign_tile_properties(run_info.SPATIAL.ACTION{ai,1}, range); %writes the provider class
        %                     end
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
            
            run_info.PARA.comment.projection_class = {'projection class providing providing information on the locations and additional data for each target point'};
            
            run_info.PARA.comment.clustering_class = {'clustering class to select representative clusters considering the properties of the different target points'};
                        
        end
        
    end
end



