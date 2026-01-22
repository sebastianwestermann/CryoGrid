classdef GENERATE_ENSEMBLE_gaussian < matlab.mixin.Copyable


    properties
        PARA
        CONST
        STATVAR
        TEMP
        PARENT
    end
    
    methods
        
 function ensemble = provide_PARA(ensemble) 

            ensemble.PARA.ensemble_size = []; %set in DA class
            ensemble.PARA.variables = [];
            ensemble.PARA.center = [];
            ensemble.PARA.width = [];
            ensemble.PARA.lower_bound = [];
            ensemble.PARA.upper_bound = [];
            ensemble.PARA.sample_std_range = [];

            ensemble.PARA.id_variable = [];
            % ensemble.PARA.tag_variable= [];
            % ensemble.PARA.tag_value = [];

            ensemble.PARA.add_existing = [];
            ensemble.PARA.mask_class = [];
            ensemble.PARA.mask_class_index = [];
            
        end
        
        function ensemble = provide_CONST(ensemble)

        end
        
        function ensemble = provide_STATVAR(ensemble)

        end 
        
        
        function ensemble = finalize_init(ensemble, tile)
            ensemble.PARENT.TEMP.ensemble_size = ensemble.PARENT.TEMP.ensemble_size .* ensemble.PARA.ensemble_size;
            rng(sum(ensemble.PARA.center) + sum(ensemble.PARA.width) + ensemble.PARA.ensemble_size);
        end

        function ensemble = generate_ensemble(ensemble)
           
            %makes initial ensemble          
            for i=1:size(ensemble.PARA.variables,1)
                ensemble.STATVAR.([ensemble.PARA.variables{i,1} '_gaussian']) = ensemble.PARA.center(i,1) + randn(ensemble.PARA.ensemble_size,1) .* ensemble.PARA.width(i,1);
                if ensemble.PARA.sample_std_range(i,1) >= 1  && ensemble.PARA.ensemble_size >=3
                    ensemble.STATVAR.([ensemble.PARA.variables{i,1} '_gaussian'])(1,1) = ensemble.PARA.center(i,1);
                    ensemble.STATVAR.([ensemble.PARA.variables{i,1} '_gaussian'])(2,1) = ensemble.PARA.center(i,1) - ensemble.PARA.width(i,1);
                    ensemble.STATVAR.([ensemble.PARA.variables{i,1} '_gaussian'])(3,1) = ensemble.PARA.center(i,1) + ensemble.PARA.width(i,1);
                end
                if ensemble.PARA.sample_std_range(i,1) >=2  && ensemble.PARA.ensemble_size >=5
                    ensemble.STATVAR.([ensemble.PARA.variables{i,1} '_gaussian'])(4,1) = ensemble.PARA.center(i,1) - 2.* ensemble.PARA.width(i,1);
                    ensemble.STATVAR.([ensemble.PARA.variables{i,1} '_gaussian'])(5,1) = ensemble.PARA.center(i,1) + 2.* ensemble.PARA.width(i,1);
                end
                if ensemble.PARA.sample_std_range(i,1) >=3  && ensemble.PARA.ensemble_size >=7
                    ensemble.STATVAR.([ensemble.PARA.variables{i,1} '_gaussian'])(6,1) = ensemble.PARA.center(i,1) - 3.* ensemble.PARA.width(i,1);
                    ensemble.STATVAR.([ensemble.PARA.variables{i,1} '_gaussian'])(7,1) = ensemble.PARA.center(i,1) + 3.* ensemble.PARA.width(i,1);
                end
                ensemble.STATVAR.(ensemble.PARA.variables{i,1}) = max(ensemble.PARA.lower_bound(i,1), min(ensemble.PARA.upper_bound(i,1), ensemble.STATVAR.([ensemble.PARA.variables{i,1} '_gaussian']) ));
                ensemble.PARENT.TEMP.mean_gaussian = [ensemble.PARENT.TEMP.mean_gaussian; ensemble.PARA.center(i,1)];
                ensemble.PARENT.TEMP.std_gaussian = [ensemble.PARENT.TEMP.std_gaussian; ensemble.PARA.width(i,1)];
                ensemble.PARENT.TEMP.variable_name = [ensemble.PARENT.TEMP.variable_name; ensemble.PARA.variables{i,1}];
                ensemble.PARENT.TEMP.ensemble_class = [ensemble.PARENT.TEMP.ensemble_class; class(ensemble)];
                ensemble.PARENT.TEMP.ensemble_class_index = [ensemble.PARENT.TEMP.ensemble_class_index; ensemble.PARA.class_index];
            end
            if ~isempty(ensemble.PARA.id_variable)
                ensemble.STATVAR.(ensemble.PARA.id_variable) = [1:ensemble.PARA.ensemble_size]';
            end
            % if ~isempty(ensemble.PARA.tag_variable)
            %     ensemble.STATVAR.(ensemble.PARA.tag_variable) = [1:ensemble.PARA.ensemble_size]'.*0+ensemble.PARA.tag_value;
            % end
        end

        %is called to make a new ensemble after DA is successful, this
        %could be used to implement learning based on the history, i.e. the mean and covariances of previous DA steps could be used  
        function ensemble = generate_ensemble2(ensemble)
            for i=1:size(ensemble.PARA.variables,1)
                ensemble.STATVAR.([ensemble.PARA.variables{i,1} '_gaussian']) = ensemble.PARA.center(i,1) + randn(ensemble.PARA.ensemble_size,1) .* ensemble.PARA.width(i,1);
                if ensemble.PARA.sample_std_range(i,1) >= 1  && ensemble.PARA.ensemble_size >=3
                    ensemble.STATVAR.([ensemble.PARA.variables{i,1} '_gaussian'])(1,1) = ensemble.PARA.center(i,1);
                    ensemble.STATVAR.([ensemble.PARA.variables{i,1} '_gaussian'])(2,1) = ensemble.PARA.center(i,1) - ensemble.PARA.width(i,1);
                    ensemble.STATVAR.([ensemble.PARA.variables{i,1} '_gaussian'])(3,1) = ensemble.PARA.center(i,1) + ensemble.PARA.width(i,1);
                end
                if ensemble.PARA.sample_std_range(i,1) >=2  && ensemble.PARA.ensemble_size >=5
                    ensemble.STATVAR.([ensemble.PARA.variables{i,1} '_gaussian'])(4,1) = ensemble.PARA.center(i,1) - 2.* ensemble.PARA.width(i,1);
                    ensemble.STATVAR.([ensemble.PARA.variables{i,1} '_gaussian'])(5,1) = ensemble.PARA.center(i,1) + 2.* ensemble.PARA.width(i,1);
                end
                if ensemble.PARA.sample_std_range(i,1) >=3  && ensemble.PARA.ensemble_size >=7
                    ensemble.STATVAR.([ensemble.PARA.variables{i,1} '_gaussian'])(6,1) = ensemble.PARA.center(i,1) - 3.* ensemble.PARA.width(i,1);
                    ensemble.STATVAR.([ensemble.PARA.variables{i,1} '_gaussian'])(7,1) = ensemble.PARA.center(i,1) + 3.* ensemble.PARA.width(i,1);
                end
                ensemble.STATVAR.(ensemble.PARA.variables{i,1}) = max(ensemble.PARA.lower_bound(i,1), min(ensemble.PARA.upper_bound(i,1), ensemble.STATVAR.([ensemble.PARA.variables{i,1} '_gaussian']) ));
                pos = find(strcmp(ensemble.PARENT.TEMP.variable_name, ensemble.PARA.variables{i,1}) & strcmp(ensemble.PARENT.TEMP.ensemble_class,  class(ensemble)) & ensemble.PARENT.TEMP.ensemble_class_index==ensemble.PARA.class_index)
                ensemble.PARENT.TEMP.mean_gaussian(pos,1) = ensemble.PARA.center(i,1);
                ensemble.PARENT.TEMP.std_gaussian(pos,1) = ensemble.PARA.width(i,1);
            end

            if ~isempty(ensemble.PARA.id_variable)
                ensemble.STATVAR.(ensemble.PARA.id_variable) = [1:ensemble.PARA.ensemble_size]';
            end
            % if ~isempty(ensemble.PARA.tag_variable)
            %     ensemble.STATVAR.(ensemble.PARA.tag_variable) = [1:ensemble.PARA.ensemble_size]'.*0+ensemble.PARA.tag_value;
            % end
        end

        function ensemble_class = generate_ensemble_from_existing(ensemble_class, proj)
            ensemble_class = generate_ensemble(ensemble_class);
            if ensemble_class.PARA.add_existing == 1 %combine with existing

                fn = fieldnames(proj.STATVAR);
                proj.STATVAR.mask = logical(proj.STATVAR.(fn{1,1}).*0+1);
                for i=1:size(ensemble_class.PARA.mask_class_index,1)
                    mask_class = copy(proj.RUN_INFO.PPROVIDER.CLASSES.(ensemble_class.PARA.mask_class{i,1}){ensemble_class.PARA.mask_class_index(i,1),1});
                    mask_class.PARENT = proj;
                    mask_class = finalize_init(mask_class);
                    mask_class = apply_mask(mask_class); %can be additive or subtractive
                end

                %reduce the list to the ones inside the masks
                mask = proj.STATVAR.mask;

                proj_variables = fieldnames(proj.STATVAR);
                ensemble_variables = fieldnames(ensemble_class.STATVAR);
                size_of_existing = sum(double(mask)); %size(proj.STATVAR.(proj_variables{1,1}), 1);
                size_of_new = size(ensemble_class.STATVAR.(ensemble_variables{1,1}), 1);

                %go through new variables
                for i=1:size(ensemble_variables,1)
                    new_param = [];
                    for j=1:size(ensemble_class.STATVAR.(ensemble_variables{i,1}),2)
                        a=ensemble_class.STATVAR.(ensemble_variables{i,1})(:,j);
                        a=repmat(a,1,size_of_existing);
                        new_param = [new_param  a(:)];
                    end
                    ensemble_class.STATVAR.(ensemble_variables{i,1}) = new_param;
                end

                %go through existing variables (have been written above if they
                %also exist in new variables)
                for i=1:size(proj_variables,1)
                    if ~any(strcmp(proj_variables{i,1}, ensemble_variables))
                        ensemble_class.STATVAR.(proj_variables{i,1}) = [];
                        for j=1:size(proj.STATVAR.(proj_variables{i,1}),2)
                            a=proj.STATVAR.(proj_variables{i,1})(mask,j);
                            a=repmat(a,1,size_of_new)';
                            ensemble_class.STATVAR.(proj_variables{i,1}) = [ensemble_class.STATVAR.(proj_variables{i,1}) a(:)];
                        end
                    end
                end
            end
        end


        function new_value = get_value_from_gaussian(ensemble, var_name, value)
            new_value = NaN;
            for i=1:size(ensemble.PARA.variables,1)
                if strcmp(ensemble.PARA.variables{i,1}, var_name)
                    new_value = max(ensemble.PARA.lower_bound(i,1), min(ensemble.PARA.upper_bound(i,1), value));
                end
            end
        end

        % %add possibility to restrict the variables only to the ones which should be changed       
        % function ensemble = recalculate_ensemble_parameters_after_DA(ensemble, tile, variable_list)
        %     %recalculate ensemble parameters based on DA weights and
        %     %previous ensemble-the new variable centers and widths are set
        %     %by the DA class, plus the variables that are affected by the
        %     %DA
        % 
        %     for i=1:size(variable_list,1)
        %         pos = find(strcmp(variable_list{i,1}, ensemble.TEMP.variable_name));
        %         ensemble.STATVAR.(variable_list{i,1}) = ensemble.TEMP.value_gaussian(pos,1);
        % 
        %     end
        % 
        %     %establish pointers to all classes in the stratigraphy
        %     ensemble = set_pointers2classes(ensemble, tile);
        % 
        %     %rewrite PARA in PROVIDER and all classes in the
        %     %stratigraphy
        %     for i=1:size(ensemble.PARA.modify_class_name,1)
        %         if any(strcmp(variable_list, ensemble.PARA.variable_in_ensemble{i,1}))
        %             tile.RUN_INFO.PPROVIDER.CLASSES.(ensemble.PARA.modify_class_name{i,1}){ensemble.PARA.modify_class_index(i,1),1}.PARA.(ensemble.PARA.variable_in_class{i,1})(ensemble.PARA.index_V(i,1), ensemble.PARA.index_H(i,1)) = ...
        %                 max(ensemble.PARA.lower_bound(i,1), min(ensemble.PARA.upper_bound(i,1), ensemble.STATVAR.(ensemble.PARA.variable_in_ensemble{i,1})));
        %             for j=1:size(ensemble.PARA.modify_class_pointer{i,1},1)
        %                 ensemble.PARA.modify_class_pointer{i,1}(j,1).PARA.(ensemble.PARA.variable_in_class{i,1})(ensemble.PARA.index_V(i,1), ensemble.PARA.index_H(i,1)) = ...
        %                     max(ensemble.PARA.lower_bound(i,1), min(ensemble.PARA.upper_bound(i,1), ensemble.STATVAR.(ensemble.PARA.variable_in_ensemble{i,1})));
        %             end
        %         end
        %     end
        % 
        %     ensemble.PARA.modify_class_pointer = [];
        % end
        
        function ensemble = set_pointers2classes(ensemble, tile)
            for i=1:size(ensemble.PARA.modify_class_name,1)
                ensemble.PARA.modify_class_pointer{i,1} = {};
            end
            
            ensemble = find_classes_in_variable(ensemble, tile, 1);
            
        end
        
        function ensemble = find_classes_in_variable(ensemble, current_class, level)
            class_name = class(current_class);
            if isobject(current_class)
                for i=1:size(ensemble.PARA.modify_class_name, 1)
                    if strcmp(ensemble.PARA.modify_class_name{i,1}, class_name)
                        if  ensemble.PARA.modify_class_index(i,1) == current_class.PARA.class_index
                            dec = 1;
                            for j=1:size(ensemble.PARA.modify_class_pointer{i,1}, 1)
                                if isequal(ensemble.PARA.modify_class_pointer{i,1}(j,1), current_class)
                                    dec = 0;
                                end
                            end
                            if dec
                                ensemble.PARA.modify_class_pointer{i,1} = [ensemble.PARA.modify_class_pointer{i,1}; current_class];
                            end
                        end
                    end
                end
            end
            if isstruct(current_class) || isobject(current_class)
                variables = fieldnames(current_class);
                for i = 1:size(variables,1)
                    if iscell(current_class.(variables{i,1}))
                        if level <= 5
                            for j=1:size(current_class.(variables{i,1}), 1)
                                for k=1:size(current_class.(variables{i,1}), 2)
                                    ensemble = find_classes_in_variable(ensemble, current_class.(variables{i,1}){j,k}, level+1);
                                end
                            end
                        end
                    elseif isstruct(current_class.(variables{i,1}))
                        if level <= 5
                            ensemble = find_classes_in_variable(ensemble, current_class.(variables{i,1}), level+1);
                        end
                    elseif isobject(current_class.(variables{i,1})) && ~strcmp(variables{i,1}, 'RUN_INFO')
                        if level <= 5
                            ensemble = find_classes_in_variable(ensemble, current_class.(variables{i,1}), level+1);
                        end
                    end
                end
            end
        end

    end
end

