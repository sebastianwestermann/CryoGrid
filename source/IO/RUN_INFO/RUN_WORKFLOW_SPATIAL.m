%========================================================================
% CryoGrid RUN_INFO class RUN_WORKFLOW_SPATIAL
%
% S. Westermann, April 2025
%========================================================================

classdef RUN_WORKFLOW_SPATIAL < matlab.mixin.Copyable
    
    properties
        PPROVIDER
        PARA
        CONST
        STATVAR
        STORE
    end
    
    
    methods
        
        
        function run_info = provide_PARA(run_info)
            
            run_info.PARA.run_info_class = [];
            run_info.PARA.run_info_class_index = [];
            run_info.PARA.spatial_exchange_class = [];
            run_info.PARA.spatial_exchange_class_index = [];

            run_info.PARA.number_of_cores = [];
            
        end
        
        function run_info = provide_CONST(run_info)
            
        end
        
        function run_info = provide_STATVAR(run_info)
            
        end
        
        
        function run_info = finalize_init(run_info)

        end
        

        function [run_info, tile] = run_model(run_info)
            if run_info.PARA.number_of_cores>1
                parpool(run_info.PARA.number_of_cores)
                spmd
                    tile = 0;
                    for i=1:size(run_info.PARA.run_info_class, 1)
                        current_run_info = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.run_info_class{i,1}){run_info.PARA.run_info_class_index(i,1),1});
                        current_run_info.PARA.parallel_pool_open = 1;
                        current_run_info.PARA.number_of_cores = run_info.PARA.number_of_cores;
                        current_run_info.PPROVIDER = run_info.PPROVIDER;

                        if ~isempty(run_info.PARA.spatial_exchange_class{i,1}) && ~(sum(isnan(run_info.PARA.spatial_exchange_class{i,1}))>0)
                            spatial_exchange = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.spatial_exchange_class{i,1}){run_info.PARA.spatial_exchange_class_index(i,1),1});
                            spatial_exchange = finalize_init(spatial_exchange, run_info);
                            current_run_info = read_spatial_from_master(spatial_exchange, current_run_info, run_info);
                        end
                       disp(['run_info class ' num2str(i)])
                        current_run_info = finalize_init(current_run_info);
                        
                        [current_run_info, tile] = run_model(current_run_info);

                        if ~isempty(run_info.PARA.spatial_exchange_class{i,1}) && ~(sum(isnan(run_info.PARA.spatial_exchange_class{i,1}))>0)
                            run_info = write_spatial_to_master(spatial_exchange, current_run_info, run_info);
                        end
                        spmdBarrier;
                    end
                end
            else
                tile = 0;
                for i=1:size(run_info.PARA.run_info_class, 1)
                    current_run_info = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.run_info_class{i,1}){run_info.PARA.run_info_class_index(i,1),1});
                    current_run_info.PARA.close_parallel_pool = 0;
                    current_run_info.PARA.number_of_cores = run_info.PARA.number_of_cores;
                    current_run_info.PPROVIDER = run_info.PPROVIDER;

                    if ~isempty(run_info.PARA.spatial_exchange_class{i,1}) && ~(sum(isnan(run_info.PARA.spatial_exchange_class{i,1}))>0)
                        spatial_exchange = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.spatial_exchange_class{i,1}){run_info.PARA.spatial_exchange_class_index(i,1),1});
                        spatial_exchange = finalize_init(spatial_exchange, run_info);
                        current_run_info = read_spatial_from_master(spatial_exchange, current_run_info, run_info);
                    end

                    current_run_info = finalize_init(current_run_info);

                    [current_run_info, tile] = run_model(current_run_info);

                    if ~isempty(run_info.PARA.spatial_exchange_class{i,1}) && ~(sum(isnan(run_info.PARA.spatial_exchange_class{i,1}))>0)
                        run_info = write_spatial_to_master(spatial_exchange, current_run_info, run_info);
                    end
                end
            end
        end
        
        
        
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



