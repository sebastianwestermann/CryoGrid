classdef GENERATE_ENSEMBLE < matlab.mixin.Copyable

    properties
        PARA
        CONST
        STATVAR
        TEMP
        PARENT
    end

    methods
        function ensemble = provide_PARA(ensemble)
            ensemble.PARA.no_variable_sets = [];
            ensemble.PARA.variable_set1 = [];
            ensemble.PARA.variable_set2 = [];
            ensemble.PARA.variable_set3 = [];
            ensemble.PARA.variable_set4 = [];
            ensemble.PARA.variable_set5 = [];
            ensemble.PARA.variable_set6 = [];
            ensemble.PARA.variable_set7 = [];
            ensemble.PARA.variable_set8 = [];
            ensemble.PARA.variable_set9 = [];
            ensemble.PARA.variable_set10 = [];

            ensemble.PARA.id_variable = [];
            ensemble.PARA.tag_variable= [];
            ensemble.PARA.tag_value= [];

            ensemble.PARA.add_existing = [];
            ensemble.PARA.mask_class = [];
            ensemble.PARA.mask_class_index = [];
        end

        function ensemble = provide_STATVAR(ensemble)
        end

        function ensemble = provide_CONST(ensemble)
        end

        function ensemble = finalize_init(ensemble)
        end

        function ensemble = generate_ensemble(ensemble)
            % Populate empty proj.STATVAR.(variables) with all possible
            % combinations of variable sets provided in GENERATE_ENSEMBLE

            % 1. add first variable set to STATVAR
            variable_set = eval(['ensemble.PARA.variable_set' num2str(1)]);
            variable_set = rmfield(variable_set,'depth');
            variable_names = fieldnames(variable_set);
            for j = 1:numel(variable_names)
                ensemble.STATVAR.(variable_names{j,1}) = variable_set.(variable_names{j,1});
            end
            %ensemble.STATVAR.key = [1:length(variable_set.(variable_names{1}))]';
            key = [1:length(variable_set.(variable_names{1}))]';
            variables = variable_names;
            
            % 2. Add combinations with remaining variable_sets
            for i = 2:ensemble.PARA.no_variable_sets
                variable_set = eval(['ensemble.PARA.variable_set' num2str(i)]);
                variable_set = rmfield(variable_set,'depth');
                variable_names = fieldnames(variable_set);
                variable_set.key = [1:length(variable_set.(variable_names{1}))]';
                
                % 2.1 find all possible combinations of variables in ensemble.PARENT.STATVAR and variable_set
                combos = combinations(key,variable_set.key);
                
                % 2.2 extend values to accomodate new no. combinantions
                for k = 1:length(variables)
                    temp = ensemble.STATVAR.(variables{k,1}); % 
                    for l = 1:size(combos,1) % Loop through all new combinations and update
                        ensemble.STATVAR.(variables{k,1})(l,1) = temp(key == combos{l,1});
                    end
                end
                
                % 2.3 Add next variable set
                for j = 1:numel(variable_names)
                    ensemble.STATVAR.(variable_names{j,1}) = variable_set.(variable_names{j,1})(combos{:,2});
                end
                key = str2num([num2str(combos{:,1}) num2str(combos{:,2})]);
                variables = [variables; variable_names(j,1)];
            end
            ensemble.TEMP.variables = variables;
            if ~isempty(ensemble.PARA.id_variable)
                ensemble.STATVAR.(ensemble.PARA.id_variable) = [1:size(ensemble.STATVAR.(variable_names{1,1}),1)]';
            end
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
        
    end
end
