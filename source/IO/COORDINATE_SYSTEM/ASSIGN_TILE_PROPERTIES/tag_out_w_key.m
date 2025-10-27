%========================================================================
% CryoGrid ASSIGN_TILE_PROPERTIES class tag_out_w_key
% add the run_number as tag to the OUT class
%
% S. Westermann, April 2025
%========================================================================

classdef tag_out_w_key < matlab.mixin.Copyable

    properties
        PARA
        CONST
        PROJ
    end
    
    methods
        function update = provide_PARA(update)

            update.PARA.target_class_name = [];
            update.PARA.target_class_index = [];
        end
        
        function update = provide_STATVAR(update)

        end
        
        function update = provide_CONST(update)
            
        end
        
        function update = finalize_init(update)
 
        end
        
        function update = assign_tile_properties(update, run_number)
            for i=1:size(update.PARA.target_class_name,1)
                update.PROJ.RUN_INFO.PPROVIDER.CLASSES.(update.PARA.target_class_name{i,1}){update.PARA.target_class_index(i,1),1}.PARA.tag = num2str(update.PROJ.RUN_INFO.SPATIAL.STATVAR.key(run_number,1));
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

