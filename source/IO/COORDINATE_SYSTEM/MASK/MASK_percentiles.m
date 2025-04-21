%========================================================================
% CryoGrid DATA_MASK class MASK_percentiles
% selects the region of interest as the target points in an altitudinal
% range
%
% S. Westermann, April 2025
%========================================================================

classdef MASK_percentiles < matlab.mixin.Copyable

    properties
        PARENT
        PARA
        CONST
        STATVAR
    end
    
    methods
        
        function mask = provide_PARA(mask)
            mask.PARA.variables = [];
            mask.PARA.min_percentile = []; % 0 = minimum of value range, 1 = maximum of value range 
            mask.PARA.max_percentile = [];
            mask.PARA.additive = [];
        end

        function mask = provide_STATVAR(mask)

        end
        
        function mask = provide_CONST(mask)
            
        end
        
        function mask = finalize_init(mask)

        end
        

        function mask = apply_mask(mask)
            valid = double(mask.PARENT.STATVAR.mask .*0);
            for i=1:size(mask.PARA.variables,1)
                [~, ind] = sort(mask.PARENT.STATVAR.(mask.PARA.variables{i,1}));
                lb = max(1, round(mask.PARA.min_percentile(i,1) .* size(mask.PARENT.STATVAR.(mask.PARA.variables{i,1}),1)));
                ub = max(1, round(mask.PARA.max_percentile(i,1) .* size(mask.PARENT.STATVAR.(mask.PARA.variables{i,1}),1)));  
                valid(ind(lb:ub,1),1) = valid(ind(lb:ub,1),1) + 1;
            end
            if mask.PARA.additive
                valid = valid >0;
            else
                valid = (valid==size(mask.PARA.variables,1))
            end

            mask.PARENT.STATVAR.mask = mask.PARENT.STATVAR.mask & valid;

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

