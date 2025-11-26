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
            ensemble.STATVAR.key = [1:length(variable_set.(variable_names{1}))]';
            variables = variable_names;
            
            % 2. Add combinations with remaining variable_sets
            for i = 2:ensemble.PARA.no_variable_sets
                variable_set = eval(['ensemble.PARA.variable_set' num2str(i)]);
                variable_set = rmfield(variable_set,'depth');
                variable_names = fieldnames(variable_set);
                variable_set.key = [1:length(variable_set.(variable_names{1}))]';
                
                % 2.1 find all possible combinations of variables in ensemble.PARENT.STATVAR and variable_set
                combos = combinations(ensemble.STATVAR.key,variable_set.key);
                
                % 2.2 extend values to accomodate new no. combinantions
                for k = 1:length(variables)
                    temp = ensemble.STATVAR.(variables{k,1}); % 
                    for l = 1:size(combos,1) % Loop through all new combinations and update
                        ensemble.STATVAR.(variables{k,1})(l,1) = temp(ensemble.STATVAR.key == combos{l,1});
                    end
                end
                
                % 2.3 Add next variable set
                for j = 1:numel(variable_names)
                    ensemble.STATVAR.(variable_names{j,1}) = variable_set.(variable_names{j,1})(combos{:,2});
                end
                
                ensemble.STATVAR.key = str2num([num2str(combos{:,1}) num2str(combos{:,2})]);
                variables = [variables; variable_names(j,1)];
            end
            variables = [variables; 'key'];
            ensemble.TEMP.variables = variables;
        end

        % function ensemble = expand_ensemble(ensemble)
        %     % Expand proj.STATVAR.(variables) with all possible combinations of variable sets provided in GENERATE_ENSEMBLE
        % 
        %     variables = ensemble.PARENT.STATVAR.variables; % Existing variables
        % 
        %     % 2. Add combinations with remaining variable_sets
        %     for i = 1:ensemble.PARA.no_variable_sets
        %         variable_set = eval(['ensemble.PARA.variable_set' num2str(i)]);
        %         variable_set = rmfield(variable_set,'depth');
        %         variable_names = fieldnames(variable_set);
        %         variable_set.key = 1:length(variable_set.(variable_names{1}));
        % 
        %         % 2.1 find all possible combinations of variables in ensemble.PARENT.STATVAR and variable_set
        %         combos = combinations(ensemble.PARENT.STATVAR.key,variable_set.key);
        % 
        %         % 2.2 extend values to accomodate new no. comibantions
        %         for k = 1:length(variables)
        %             temp = ensemble.PARENT.STATVAR.(variables{k,1}); % 
        %             for l = 1:size(combos,1) % Loop through all new combinations and update
        %                 ensemble.PARENT.STATVAR.(variables{k,1})(l) = temp(ensemble.PARENT.STATVAR.key == combos{l,1});
        %             end
        %         end
        % 
        %         % 2.3 Add next variable set
        %         for j = 1:numel(variable_names)
        %             ensemble.PARENT.STATVAR.(variable_names{j,1}) = variable_set.(variable_names{j,1})(combos{:,2});
        %         end
        % 
        %         ensemble.PARENT.STATVAR.key = str2num([num2str(combos{:,1}) num2str(combos{:,2})]);
        %         variables = [variables; variable_names(:,1)];
        %     end
        %     ensemble.PARENT.STATVAR.variables = variables;
        % end

    end
end
