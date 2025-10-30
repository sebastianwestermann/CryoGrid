classdef GENERATE_ENSEMBLE < matlab.mixin.Copyable

    properties
        PARA
        CONST
        STATVAR
        TEMP
        PARENT
    end

    methods
        function dem = provide_PARA(dem)
            dem.PARA.no_variable_sets = [];
            dem.PARA.variable_set1 = [];
            dem.PARA.variable_set2 = [];
            dem.PARA.variable_set3 = [];
            dem.PARA.variable_set4 = [];
            dem.PARA.variable_set5 = [];

        end

        function dem = provide_STATVAR(dem)
        end

        function dem = provide_CONST(dem)
        end

        function dem = finalize_init(dem)
        end

        function dem = generate_ensemble(dem)
            % Populate empty proj.STATVAR.(variables) with all possible
            % combinations of variable sets provided in GENERATE_ENSEMBLE


            % 1. add first variable set to STATVAR
            variable_set = eval(['dem.PARA.variable_set' num2str(1)]);
            variable_set = rmfield(variable_set,'depth');
            variable_names = fieldnames(variable_set);
            for j = 1:numel(variable_names)
                dem.PARENT.STATVAR.(variable_names{j,1}) = variable_set.(variable_names{j,1});
            end
            dem.PARENT.STATVAR.key = 1:length(variable_set.(variable_names{1}));
            variables = variable_names;
            
            % 2. Add combinations with remaining variable_sets
            for i = 2:dem.PARA.no_variable_sets
                variable_set = eval(['dem.PARA.variable_set' num2str(i)]);
                variable_set = rmfield(variable_set,'depth');
                variable_names = fieldnames(variable_set);
                variable_set.key = 1:length(variable_set.(variable_names{1}));
                
                % 2.1 find all possible combinations of variables in dem.PARENT.STATVAR and variable_set
                combos = combinations(dem.PARENT.STATVAR.key,variable_set.key);
                
                % 2.2 extend values to accomodate new no. comibantions
                for k = 1:length(variables)
                    temp = dem.PARENT.STATVAR.(variables{k,1}); % 
                    for l = 1:size(combos,1) % Loop through all new combinations and update
                        dem.PARENT.STATVAR.(variables{k,1})(l) = temp(dem.PARENT.STATVAR.key == combos{l,1});
                    end
                end
                
                % 2.3 Add next variable set
                for j = 1:numel(variable_names)
                    dem.PARENT.STATVAR.(variable_names{j,1}) = variable_set.(variable_names{j,1})(combos{:,2});
                end
                
                dem.PARENT.STATVAR.key = str2num([num2str(combos{:,1}) num2str(combos{:,2})]);
                variables = [variables; variable_names(j,1)];
            end
            dem.PARENT.TEMP.variables = variables;
        end

        function dem = expand_ensemble(dem)
            % Expand proj.STATVAR.(variables) with all possible combinations of variable sets provided in GENERATE_ENSEMBLE
            
            variables = dem.PARENT.STATVAR.variables; % Existing variables
            
            % 2. Add combinations with remaining variable_sets
            for i = 1:dem.PARA.no_variable_sets
                variable_set = eval(['dem.PARA.variable_set' num2str(i)]);
                variable_set = rmfield(variable_set,'depth');
                variable_names = fieldnames(variable_set);
                variable_set.key = 1:length(variable_set.(variable_names{1}));
                
                % 2.1 find all possible combinations of variables in dem.PARENT.STATVAR and variable_set
                combos = combinations(dem.PARENT.STATVAR.key,variable_set.key);
                
                % 2.2 extend values to accomodate new no. comibantions
                for k = 1:length(variables)
                    temp = dem.PARENT.STATVAR.(variables{k,1}); % 
                    for l = 1:size(combos,1) % Loop through all new combinations and update
                        dem.PARENT.STATVAR.(variables{k,1})(l) = temp(dem.PARENT.STATVAR.key == combos{l,1});
                    end
                end
                
                % 2.3 Add next variable set
                for j = 1:numel(variable_names)
                    dem.PARENT.STATVAR.(variable_names{j,1}) = variable_set.(variable_names{j,1})(combos{:,2});
                end
                
                dem.PARENT.STATVAR.key = str2num([num2str(combos{:,1}) num2str(combos{:,2})]);
                variables = [variables; variable_names(:,1)];
            end
            dem.PARENT.STATVAR.variables = variables;
        end

    end
end
