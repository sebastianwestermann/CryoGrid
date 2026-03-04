%========================================================================
% CryoGrid ASSIGN_TILE_PROPERTIES class tag_out_w_run_number
% add the run_number as tag to the OUT class
%
% S. Westermann, Dec 2022
%========================================================================

classdef tag_out_w_run_number < matlab.mixin.Copyable

    properties
        PARA
        CONST
        PROJ
    end
    
    methods
        function update = provide_PARA(update)

            update.PARA.target_class_name = [];
            update.PARA.target_class_index = [];
            update.PARA.append = 0;
            update.PARA.overwrite_existing = 0;
        end
        
        function update = provide_STATVAR(update)

        end
        
        function update = provide_CONST(update)
            
        end
        
        function update = finalize_init(update)

        end
        
        function update = assign_tile_properties(update, run_number)
            for i=1:size(update.PARA.target_class_name,1)
                str = update.PROJ.RUN_INFO.PPROVIDER.CLASSES.(update.PARA.target_class_name{i,1}){update.PARA.target_class_index(i,1),1}.PARA.tag2;
                if update.PARA.append == 0 || isempty(str) || sum(isnan(str))>0
                    update.PROJ.RUN_INFO.PPROVIDER.CLASSES.(update.PARA.target_class_name{i,1}){update.PARA.target_class_index(i,1),1}.PARA.tag2 = num2str(run_number);
                else
                    update.PROJ.RUN_INFO.PPROVIDER.CLASSES.(update.PARA.target_class_name{i,1}){update.PARA.target_class_index(i,1),1}.PARA.tag2 = [str '_' num2str(run_number)];
                end
            end

            if update.PARA.overwrite_existing == 1 %go through existing classes and fid the oes to change
                if any(strcmp(fieldnames(update.PROJ.RUN_INFO), 'TILE'))
                    if ~isempty(update.PROJ.RUN_INFO.TILE) && any(strcmp(fieldnames(update.PROJ.RUN_INFO.TILE), 'OUT'))
                        for i=1:size(update.PARA.target_class_name,1)
                            for j=1:size(update.PROJ.RUN_INFO.TILE.OUT, 1)
                                if strcmp(class(update.PROJ.RUN_INFO.TILE.OUT{j,1}), update.PARA.target_class_name{i,1}) && update.PROJ.RUN_INFO.TILE.OUT{j,1}.PARA.class_index == update.PARA.target_class_index(i,1)
                                    update.PROJ.RUN_INFO.TILE.OUT{j,1}.PARA.tag2 = num2str(run_number);
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

