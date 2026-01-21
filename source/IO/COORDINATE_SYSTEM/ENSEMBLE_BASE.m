%========================================================================
% CryoGrid COORDINATE_SYSTEM class ENSEMBLE_BASE, based on
% ENSEMBLE_OF_POINTS
% S. Westerman, Dec 2025
%========================================================================


classdef ENSEMBLE_BASE < matlab.mixin.Copyable

    properties
        RUN_INFO
        PARA
        CONST
        STATVAR
        TEMP
        ACTION
    end
    
    methods
        function proj = provide_PARA(proj)

            proj.PARA.parameter_class = [];
            proj.PARA.parameter_class_index = [];
            proj.PARA.parameter_class_additive = [];
            
            proj.PARA.assign_tile_properties_class = [];
            proj.PARA.assign_tile_properties_class_index = [];

            proj.PARA.new_reference = 1;
        end
        
        function proj = provide_STATVAR(proj)

        end
        
        function proj = provide_CONST(proj)
            
        end
        
        function proj = finalize_init(proj)

            proj.TEMP.mean_gaussian = [];
            proj.TEMP.std_gaussian = [];
            proj.TEMP.variable_name = {};
            proj.TEMP.ensemble_class = {};
            proj.TEMP.ensemble_class_index = [];
            proj.TEMP.ensemble_size = 1;

            if ~isempty(proj.PARA.parameter_class)
                for i=1:size(proj.PARA.parameter_class,1)
                    parameter_class = copy(proj.RUN_INFO.PPROVIDER.CLASSES.(proj.PARA.parameter_class{i,1}){proj.PARA.parameter_class_index(i,1),1});
                    parameter_class.PARENT = proj;
                    parameter_class = finalize_init(parameter_class);
                    if proj.PARA.parameter_class_additive(i,1)
                        parameter_class = generate_ensemble_from_existing(parameter_class, proj); %independent ensemble generated, ready to be added to the existing
                        proj.TEMP.update_key = proj.PARA.parameter_class_additive(i,1); 
                        proj = expand_STATVAR(proj, parameter_class);
                    else
                        parameter_class = generate_ensemble_from_existing(parameter_class, proj); %combined ensemble generated
                        proj.STATVAR = parameter_class.STATVAR;
                    end
                end
            end
            
            for i=1:size(proj.PARA.assign_tile_properties_class,1)
                proj.ACTION{i,1} = copy(proj.RUN_INFO.PPROVIDER.CLASSES.(proj.PARA.assign_tile_properties_class{i,1}){proj.PARA.assign_tile_properties_class_index(i,1),1});
                proj.ACTION{i,1} = finalize_init(proj.ACTION{i,1});
                proj.ACTION{i,1}.PROJ = proj;
            end

        end
        

        function proj = expand_STATVAR(proj, parameter_class)
            if ~isempty(proj.STATVAR)
                proj_variables = fieldnames(proj.STATVAR);
                size_of_existing = size(proj.STATVAR.(proj_variables{1,1}), 1);
            else
                proj_variables={};
                size_of_existing = 0;
            end
            ensemble_variables = fieldnames(parameter_class.STATVAR);
            size_of_new = size(parameter_class.STATVAR.(ensemble_variables{1,1}), 1);
            
            %go through existing variables
            for i=1:size(proj_variables,1)
                if any(strcmp(proj_variables{i,1}, ensemble_variables)) %both variables exist
                    if strcmp(proj_variables{i,1}, 'key') &&  ~isempty(proj.STATVAR.key) && proj.TEMP.update_key == 2
                        proj.STATVAR.key = [proj.STATVAR.key; proj.STATVAR.key(end,1) + parameter_class.STATVAR.key];
                    else
                        proj.STATVAR.(proj_variables{i,1}) = [proj.STATVAR.(proj_variables{i,1}); parameter_class.STATVAR.(proj_variables{i,1})];
                    end
                else %variable does not exist in new ensemble
                     proj.STATVAR.(proj_variables{i,1}) = [proj.STATVAR.(proj_variables{i,1}); repmat(NaN, size_of_new,size(proj.STATVAR.(proj_variables{i,1}), 2))];
                end
            end
            %go through new variables and add non-existing ones
            for i=1:size(ensemble_variables,1)
                if ~any(strcmp(ensemble_variables{i,1}, proj_variables)) %variables does not exists in existing ensemble
                     proj.STATVAR.(ensemble_variables{i,1}) = [ repmat(NaN, size_of_existing, size(parameter_class.STATVAR.(ensemble_variables{i,1}),2)) ; parameter_class.STATVAR.(ensemble_variables{i,1})];
                end
            end
        end

        function ensemble = update_ensemble_after_optimization(ensemble, opt_class)
            STATVAR2 = ensemble.STATVAR;

            valid = find(STATVAR2.iteration == opt_class.TEMP.num_iterations);
            variables = fieldnames(STATVAR2);
            for i=1:size(variables,1)
                STATVAR2.(variables{i,1}) = STATVAR2.(variables{i,1})(valid,:);
            end
            STATVAR2.iteration = STATVAR2.iteration+1;
            
            for i=1:size(ensemble.PARA.parameter_class,1)
                if isfield(ensemble.RUN_INFO.PPROVIDER.CLASSES.(ensemble.PARA.parameter_class{i,1}){ensemble.PARA.parameter_class_index(i,1),1}.PARA, 'id_variable') && ...
                        strcmp(ensemble.RUN_INFO.PPROVIDER.CLASSES.(ensemble.PARA.parameter_class{i,1}){ensemble.PARA.parameter_class_index(i,1),1}.PARA.id_variable, opt_class.PARA.ensemble_variable_id)
                    ensemble_gen_class = copy(ensemble.RUN_INFO.PPROVIDER.CLASSES.(ensemble.PARA.parameter_class{i,1}){ensemble.PARA.parameter_class_index(i,1),1});
                end
            end

            for j=1:max(STATVAR2.(opt_class.PARA.ensemble_variable_id))
                for i=1:size(opt_class.PARA.ensemble_variables,1)
                    range = find(STATVAR2.(opt_class.PARA.ensemble_variable_id) == j);
                    STATVAR2.([opt_class.PARA.ensemble_variables{i,1} '_gaussian'])(range,1) = opt_class.TEMP.value_gaussian_resampled(i, j);
                    STATVAR2.(opt_class.PARA.ensemble_variables{i,1})(range,1) = get_value_from_gaussian(ensemble_gen_class, opt_class.PARA.ensemble_variables{i,1}, opt_class.TEMP.value_gaussian_resampled(i, j));
                end
            end

            %append (only when continuing with iterations9
            for i=1:size(variables,1)
                 ensemble.STATVAR.(variables{i,1}) =  [ensemble.STATVAR.(variables{i,1}); STATVAR2.(variables{i,1})];
            end
        end

        function ensemble = update_ensemble_after_optimization_new_timeSlice(ensemble, opt_class)
            STATVAR2 = ensemble.STATVAR;

            valid = find(STATVAR2.iteration == 1);
            variables = fieldnames(STATVAR2);
            for i=1:size(variables,1)
                STATVAR2.(variables{i,1}) = STATVAR2.(variables{i,1})(valid,:);
            end

            for i=1:size(ensemble.PARA.parameter_class,1)
                if isfield(ensemble.RUN_INFO.PPROVIDER.CLASSES.(ensemble.PARA.parameter_class{i,1}){ensemble.PARA.parameter_class_index(i,1),1}.PARA, 'id_variable') && ...
                        strcmp(ensemble.RUN_INFO.PPROVIDER.CLASSES.(ensemble.PARA.parameter_class{i,1}){ensemble.PARA.parameter_class_index(i,1),1}.PARA.id_variable, opt_class.PARA.ensemble_variable_id)
                    ensemble_gen_class = copy(ensemble.RUN_INFO.PPROVIDER.CLASSES.(ensemble.PARA.parameter_class{i,1}){ensemble.PARA.parameter_class_index(i,1),1});
                end
            end

            ensemble_gen_class.PARENT = ensemble;
            ensemble_gen_class = generate_ensemble2(ensemble_gen_class);

            for j=1:max(STATVAR2.(opt_class.PARA.ensemble_variable_id))
                for i=1:size(opt_class.PARA.ensemble_variables,1)
                    range = find(STATVAR2.(opt_class.PARA.ensemble_variable_id) == j);
                    STATVAR2.([opt_class.PARA.ensemble_variables{i,1} '_gaussian'])(range,1) = ensemble.STATVAR.([opt_class.PARA.ensemble_variables{i,1} '_gaussian'])(j,:);
                    STATVAR2.(opt_class.PARA.ensemble_variables{i,1})(range,1) = ensemble.STATVAR.(opt_class.PARA.ensemble_variables{i,1})(j,:);
                end
            end
            
            %could also append here
            for i=1:size(variables,1)
                ensemble.STATVAR.(variables{i,1}) = STATVAR2.(variables{i,1});
                % ensemble.STATVAR.(variables{i,1}) =  [ensemble.STATVAR.(variables{i,1}); STATVAR2.(variables{i,1})];
            end
        end



 
%         %-------------param file generation-----
%         function proj = param_file_info(proj)
%             proj = provide_PARA(proj);
%             
%             proj.PARA.STATVAR = [];
%             proj.PARA.class_category = 'SPATIAL_REFERENCE';
%             proj.PARA.default_value = [];
%             
%             proj.PARA.comment.max_lat = {'maximum latitude of model domain'}; 
%             
%             proj.PARA.comment.min_lat = {'minimum latitude of model domain'};
%             
%             proj.PARA.comment.lat_grid_cell_size = {'latitudinal spacing of grid'};
%             
%             proj.PARA.comment.max_lon = {'maximum longitude of model domain'};
%             
%             proj.PARA.comment.min_lon = {'minimum longitude of model domain'};
%             
%             proj.PARA.comment.lon_grid_cell_size = {'longitudinal spacing of grid'};
%             
%             proj.PARA.comment.mask_class = {'list of mask classes, constains the region of interest based on the coordinates'};
%             proj.PARA.options.mask_class.name = 'H_LIST';
%             proj.PARA.options.mask_class_index.name = 'H_LIST';
%             
%             proj.PARA.comment.data_class = {'list of data provider classes, provide data for each target location'};
%             proj.PARA.options.data_class.name = 'H_LIST';
%             proj.PARA.options.data_class_index.name = 'H_LIST';
%             
%             proj.PARA.comment.data_mask_class = {'list of data mask classes, constrains the region of interest based on the data provided by data provider classes'};
%             proj.PARA.options.data_mask_class.name = 'H_LIST';
%             proj.PARA.options.data_mask_class_index.name = 'H_LIST';
%             
%             proj.PARA.comment.assign_tile_properties_class = {'translates the data to changes to the classes used in the simulations for each point, basically cutomizing the "parameter file" for each taret location'};
%             proj.PARA.options.assign_tile_properties_class.name = 'H_LIST';
%             proj.PARA.options.assign_tile_properties_class.entries_x = {'update_one2one' 'tag_out_w_run_number'};            
%         end
        
    end
end

