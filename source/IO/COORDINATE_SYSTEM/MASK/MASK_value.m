%========================================================================
% CryoGrid DATA_MASK class MASK_value
% selects the region of interest as the target points with defined values
%
% S. Westermann, Dec 2025
%========================================================================

classdef MASK_value < matlab.mixin.Copyable

    properties
        PARENT
        PARA
        CONST
        STATVAR
    end
    
    methods
        
        function mask = provide_PARA(mask)
            mask.PARA.variables = [];
            mask.PARA.value = [];
        end

        function mask = provide_STATVAR(mask)

        end
        
        function mask = provide_CONST(mask)
            
        end
        
        function mask = finalize_init(mask)

        end
        

        function mask = apply_mask(mask)
            for i=1:size(mask.PARA.variables,1)
                if isnan(mask.PARA.value(i,1))
                    mask.PARENT.STATVAR.mask = mask.PARENT.STATVAR.mask & isnan(mask.PARENT.STATVAR.(mask.PARA.variables{i,1}));
                else
                    mask.PARENT.STATVAR.mask = mask.PARENT.STATVAR.mask & mask.PARENT.STATVAR.(mask.PARA.variables{i,1}) == mask.PARA.value(i,1);
                end
            end
        end
        
        
        
        %-------------param file generation-----
        function mask = param_file_info(mask)
            mask = provide_PARA(mask);
            
            mask.PARA.STATVAR = [];
            mask.PARA.class_category = 'DATA_MASK';
            mask.PARA.default_value = [];
            mask.PARA.options = [];
            
            mask.PARA.comment.min_altitude = {'minimum altitude of target points'};
            
            mask.PARA.comment.max_altitude = {'maximum altitude of target points'};
        end
     
            
    end
end

