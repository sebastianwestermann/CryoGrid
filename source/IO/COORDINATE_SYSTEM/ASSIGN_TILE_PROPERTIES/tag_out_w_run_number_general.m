%========================================================================
% CryoGrid ASSIGN_TILE_PROPERTIES class tag_out_w_run_number
% add the run_number as tag to the OUT class
%
% S. Westermann, Dec 2022
%========================================================================

classdef tag_out_w_run_number_general < matlab.mixin.Copyable

    properties
        PARA
        CONST
        PROJ
    end
    
    methods
        function update = provide_PARA(update)
            
            update.PARA.variables = [];
            update.PARA.add_iteration = 0;
            update.PARA.target_class_name = [];
            update.PARA.target_class_index = [];
            update.PARA.overwrite_existing = 0;
        end
        
        function update = provide_STATVAR(update)

        end
        
        function update = provide_CONST(update)
            
        end
        
        function update = finalize_init(update)
            if ~iscell(update.PARA.variables) && sum(isnan(update.PARA.variables))>0
                update.PARA.variables = [];
            end
        end
        
        function update = assign_tile_properties(update, run_number)

            
            tag2_str = [];
            if update.PARA.add_iteration == 1
                if isfield(update.PROJ.STATVAR, 'iteration')
                    id = max(update.PROJ.STATVAR.iteration);
                    valid = find(update.PROJ.STATVAR.iteration == id);
                    for i=1:size(update.PARA.variables,1)
                        tag2_str = [tag2_str num2str(update.PROJ.STATVAR.(update.PARA.variables{i,1})(valid(run_number(1)),1)) '_'];
                    end
                    tag2_str = [tag2_str num2str(run_number(1)) '_' num2str(id)]; %add run_number and iteration
                else
                    for i=1:size(update.PARA.variables,1)
                        tag2_str = [tag2_str num2str(update.PROJ.STATVAR.(update.PARA.variables{i,1})(run_number(1),1)) '_'];
                    end
                    tag2_str = [tag2_str num2str(run_number(1)) '_1']; %add run_number and iteration                    
                end
            else
                for i=1:size(update.PARA.variables,1)
                    tag2_str = [tag2_str num2str(update.PROJ.STATVAR.(update.PARA.variables{i,1})(run_number(1),1)) '_'];
                end
                tag2_str = [tag2_str num2str(run_number(1))]; %add run_number 
            end


            for i=1:size(update.PARA.target_class_name,1)
                update.PROJ.RUN_INFO.PPROVIDER.CLASSES.(update.PARA.target_class_name{i,1}){update.PARA.target_class_index(i,1),1}.PARA.tag2 = tag2_str;
            end

            if update.PARA.overwrite_existing == 1 %go through existing classes and find the ones to change
                if any(strcmp(fieldnames(update.PROJ.RUN_INFO), 'TILE'))
                    if ~isempty(update.PROJ.RUN_INFO.TILE) && any(strcmp(fieldnames(update.PROJ.RUN_INFO.TILE), 'OUT'))
                        for i=1:size(update.PARA.target_class_name,1)
                            for j=1:size(update.PROJ.RUN_INFO.TILE.OUT, 1)
                                if strcmp(class(update.PROJ.RUN_INFO.TILE.OUT{j,1}), update.PARA.target_class_name{i,1}) && update.PROJ.RUN_INFO.TILE.OUT{j,1}.PARA.class_index == update.PARA.target_class_index(i,1)
                                    update.PROJ.RUN_INFO.TILE.OUT{j,1}.PARA.tag2 = tag2_str;
                                end
                            end
                        end
                    end
                end
            end
       
        end
        
        
        %-------------param file generation-----
        function update = param_file_info(update)
            update = provide_PARA(update);

            update.PARA.STATVAR = [];
            update.PARA.class_category = 'ASSIGN_TILE_PROPERTIES';
            update.PARA.default_value = [];
            
            update.PARA.comment.target_class_name =  {'name of OUT class in which tag is to be added'};
            update.PARA.options.target_class_name.name =  'H_LIST';
            
            update.PARA.comment.target_class_index =  {'index of OUT class in which tag is to be added'};
            update.PARA.options.target_class_index.name =  'H_LIST';

        end
    end
end

