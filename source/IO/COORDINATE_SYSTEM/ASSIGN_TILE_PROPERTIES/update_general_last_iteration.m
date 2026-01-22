%========================================================================
% CryoGrid ASSIGN_TILE_PROPERTIES class update_one2one
% updates variables that have the same name in the CryoGrid class as
% provided by the DATA_PROVODER class
%
% S. Westermann, Dec 2022
%========================================================================

classdef update_general_last_iteration < matlab.mixin.Copyable

    properties
        PARA
        CONST
        PROJ
    end
    
    methods
        function update = provide_PARA(update)

            update.PARA.target_class_name = [];
            update.PARA.target_class_index = [];
            update.PARA.variable = [];
            update.PARA.variable_in_class = [];
            
        end
       
        function update = provide_STATVAR(update)

        end
        
        function update = provide_CONST(update)
            
        end
        
        function update = finalize_init(update)
 
        end
        
        function update = assign_tile_properties(update, run_number)
            
            update = set_pointers2classes(update, update.PROJ.RUN_INFO);

            if isfield(update.PROJ.STATVAR, 'iteration')
                id = max(update.PROJ.STATVAR.iteration);
                valid = find(update.PROJ.STATVAR.iteration == id);
                for i=1:size(update.PARA.target_class_name,1)
                    if ~iscell(update.PROJ.STATVAR.(update.PARA.variable{i,1})(valid(run_number,1),:))
                        update.PROJ.RUN_INFO.PPROVIDER.CLASSES.(update.PARA.target_class_name{i,1}){update.PARA.target_class_index(i,1),1}.PARA.(update.PARA.variable_in_class{i,1}) = ...
                            update.PROJ.STATVAR.(update.PARA.variable{i,1})(valid(run_number,1),:);
                        for j=1:size(update.PARA.modify_class_pointer{i,1},1)
                            update.PARA.modify_class_pointer{i,1}(j,1).PARA.(update.PARA.variable_in_class{i,1}) = update.PROJ.STATVAR.(update.PARA.variable{i,1})(valid(run_number,1),:);
                        end
                    else
                        update.PROJ.RUN_INFO.PPROVIDER.CLASSES.(update.PARA.target_class_name{i,1}){update.PARA.target_class_index(i,1),1}.PARA.(update.PARA.variable_in_class{i,1}) = ...
                            update.PROJ.STATVAR.(update.PARA.variable{i,1}){valid(run_number,1),1};
                        for j=1:size(update.PARA.modify_class_pointer{i,1},1)
                            update.PARA.modify_class_pointer{i,1}(j,1).PARA.(update.PARA.variable_in_class{i,1}) = update.PROJ.STATVAR.(update.PARA.variable{i,1}){valid(run_number,1),1};
                        end
                    end
                end
            else %make compatible for runs without iterations
                for i=1:size(update.PARA.target_class_name,1)
                    if ~iscell(update.PROJ.STATVAR.(update.PARA.variable{i,1})(run_number,:))
                        update.PROJ.RUN_INFO.PPROVIDER.CLASSES.(update.PARA.target_class_name{i,1}){update.PARA.target_class_index(i,1),1}.PARA.(update.PARA.variable_in_class{i,1}) = ...
                            update.PROJ.STATVAR.(update.PARA.variable{i,1})(run_number,:);
                        for j=1:size(update.PARA.modify_class_pointer{i,1},1)
                            update.PARA.modify_class_pointer{i,1}(j,1).PARA.(update.PARA.variable_in_class{i,1}) = update.PROJ.STATVAR.(update.PARA.variable{i,1})(run_number,:);
                        end
                    else
                        update.PROJ.RUN_INFO.PPROVIDER.CLASSES.(update.PARA.target_class_name{i,1}){update.PARA.target_class_index(i,1),1}.PARA.(update.PARA.variable_in_class{i,1}) = ...
                            update.PROJ.STATVAR.(update.PARA.variable{i,1}){run_number,1};
                        for j=1:size(update.PARA.modify_class_pointer{i,1},1)
                            update.PARA.modify_class_pointer{i,1}(j,1).PARA.(update.PARA.variable_in_class{i,1}) = update.PROJ.STATVAR.(update.PARA.variable{i,1}){run_number,1,1};
                        end
                    end
                end

            end
            update.PARA.modify_class_pointer = [];
        end

        
        function update = set_pointers2classes(update, run_info)
            for i=1:size(update.PARA.target_class_name,1)
                update.PARA.modify_class_pointer{i,1} = {};
            end
            
            update = find_classes_in_variable(update, run_info.TILE, 1);
            %run this only from TILE level and below
        end
        
        function update = find_classes_in_variable(update, current_class, level)
            class_name = class(current_class);
            if isobject(current_class)
                for i=1:size(update.PARA.target_class_name, 1)
                    if strcmp(update.PARA.target_class_name{i,1}, class_name)
                        if  update.PARA.target_class_index(i,1) == current_class.PARA.class_index
                            dec = 1;
                            for j=1:size(update.PARA.modify_class_pointer{i,1}, 1)
                                if isequal(update.PARA.modify_class_pointer{i,1}(j,1), current_class)
                                    dec = 0;
                                end
                            end
                            if dec
                                update.PARA.modify_class_pointer{i,1} = [update.PARA.modify_class_pointer{i,1}; current_class];
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
                                    update = find_classes_in_variable(update, current_class.(variables{i,1}){j,k}, level+1);
                                end
                            end
                        end
                    elseif isstruct(current_class.(variables{i,1}))
                        if level <= 5
                            update = find_classes_in_variable(update, current_class.(variables{i,1}), level+1);
                        end
                    elseif isobject(current_class.(variables{i,1})) && ~strcmp(variables{i,1}, 'RUN_INFO')
                        if level <= 5
                            update = find_classes_in_variable(update, current_class.(variables{i,1}), level+1);
                        end
                    end
                end
            end
        end
        
        
%         %-------------param file generation-----
%         function update = param_file_info(update)
%             update = provide_PARA(update);
% 
%             update.PARA.STATVAR = [];
%             update.PARA.class_category = 'ASSIGN_TILE_PROPERTIES';
%             update.PARA.default_value = [];
%             
%             update.PARA.comment.class_name =  {'name of class in which the variable is changed'};
%             update.PARA.options.class_name.name =  'H_LIST';
%             update.PARA.options.class_name.entries_x = {'TILE_1D_standard' 'TILE_1D_standard' 'TILE_1D_standard'};
%             
%             update.PARA.comment.class_index =  {'index of class in which the variable is changed'};
%             update.PARA.options.class_index.name =  'H_LIST';
%             update.PARA.options.class_index.entries_x = {'1' '1' '1'};
% 
%             update.PARA.comment.variable =  {'name of variable which is changed (must be assigned by DATA PROVIDER class)'};
%             update.PARA.options.variable.name =  'H_LIST';
%             update.PARA.options.variable.entries_x = {'latitude' 'longitude' 'altitude'};
% 
%         end
        
    end
end

