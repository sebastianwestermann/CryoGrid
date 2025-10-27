classdef ENSEMBLE_MULTITILE_GRID_AND_SUBGRID < matlab.mixin.Copyable


    properties
        PARA
        CONST
        STATVAR
        TEMP
    end
    
    methods
        
        function ensemble = provide_PARA(ensemble) 
            ensemble.PARA.grid_ensemble_class = [];
            ensemble.PARA.grid_ensemble_class_index = [];
            ensemble.PARA.grid_ensemble_size = [];
            ensemble.PARA.subgrid_ensemble_class = [];
            ensemble.PARA.subgrid_ensemble_class_index = [];
            ensemble.PARA.subgrid_ensemble_size = [];
        end
        
        function ensemble = provide_CONST(ensemble)

        end
        
        function ensemble = provide_STATVAR(ensemble)

        end 
        
        
        function ensemble = finalize_init(ensemble, tile)
            ensemble.TEMP.grid_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(ensemble.PARA.grid_ensemble_class){ensemble.PARA.grid_ensemble_class_index,1});
            tile.PARA.ensemble_size = ensemble.PARA.grid_ensemble_size;
            ensemble.TEMP.grid_class = finalize_init2(ensemble.TEMP.grid_class, tile);
            ensemble.TEMP.subgrid_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(ensemble.PARA.subgrid_ensemble_class){ensemble.PARA.subgrid_ensemble_class_index,1});
            tile.PARA.ensemble_size = ensemble.PARA.subgrid_ensemble_size;
            ensemble.TEMP.subgrid_class = finalize_init2(ensemble.TEMP.subgrid_class, tile);
            tile.PARA.ensemble_size = ensemble.PARA.grid_ensemble_size .* ensemble.PARA.subgrid_ensemble_size;

            
            variables_grid = fieldnames(ensemble.TEMP.grid_class.STATVAR);
            variables_subgrid = fieldnames(ensemble.TEMP.subgrid_class.STATVAR);

            for i=1:size(variables_grid,1)
                var_reshaped = reshape(ensemble.TEMP.grid_class.STATVAR.(variables_grid{i,1}),tile.PARA.number_of_realizations, 1, size(ensemble.TEMP.grid_class.STATVAR.(variables_grid{i,1}),2)./tile.PARA.number_of_realizations);
                var_reshaped = repmat(var_reshaped,1,ensemble.PARA.subgrid_ensemble_size,1);
                ensemble.TEMP.grid_class.STATVAR.(variables_grid{i,1}) = var_reshaped(:)';
                %ensemble.TEMP.grid_class.STATVAR.(variables_grid{i,1}) = repmat(ensemble.TEMP.grid_class.STATVAR.(variables_grid{i,1}), 1, size_of_subgrid_ensemble);
            end
            for i=1:size(variables_subgrid,1)
                var_reshaped = reshape(ensemble.TEMP.subgrid_class.STATVAR.(variables_subgrid{i,1}),tile.PARA.number_of_realizations, size(ensemble.TEMP.subgrid_class.STATVAR.(variables_subgrid{i,1}),2)./tile.PARA.number_of_realizations);
                var_reshaped = repmat(var_reshaped,1,1,ensemble.PARA.grid_ensemble_size);
                ensemble.TEMP.subgrid_class.STATVAR.(variables_subgrid{i,1}) = var_reshaped(:)';
                %ensemble.TEMP.subgrid_class.STATVAR.(variables_subgrid{i,1}) = repmat(ensemble.TEMP.subgrid_class.STATVAR.(variables_subgrid{i,1}), size_of_grid_ensemble, 1);
                %ensemble.TEMP.subgrid_class.STATVAR.(variables_subgrid{i,1}) = reshape(ensemble.TEMP.subgrid_class.STATVAR.(variables_subgrid{i,1}), 1, size_of_grid_ensemble .* size_of_subgrid_ensemble);
            end
            
            ensemble.TEMP.grid_class = write2provider(ensemble.TEMP.grid_class, tile);
            ensemble.TEMP.subgrid_class = write2provider(ensemble.TEMP.subgrid_class, tile);
            
            variables = fieldnames(ensemble.TEMP.grid_class.TEMP);
            for i=1:size(variables,1)
                ensemble.TEMP.(variables{i,1}) = ensemble.TEMP.grid_class.TEMP.(variables{i,1});
            end
        end
        
        %called by DA
        function result = get_best_fitting_parameters(ensemble, tile, variable_list, weights)
            result = get_best_fitting_parameters(ensemble.TEMP.grid_class, tile, variable_list, weights);
        end

    end
end

