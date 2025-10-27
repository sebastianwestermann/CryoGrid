%TRANSFORM class doy_linear_fit

% S. Westermann Dec 2023


classdef doy_linear_fit < matlab.mixin.Copyable
    
    properties
        PARA
        STATVAR
        TEMP
        CONST
    end
    
    methods
        
        function transform = provide_PARA(transform)
           transform.PARA.window = []; %in days
           transform.PARA.invalid_threshold = [];
        end
        
        function transform = provide_CONST(transform)
            
        end
        
        function transform = provide_STATVAR(transform)

        end
        
        function transform = finalize_init(transform, tile)
            if isempty(transform.PARA.invalid_threshold) || isnan(transform.PARA.invalid_threshold)
                 transform.PARA.invalid_threshold = Inf;
            end
        end
        
        function transform = fit_transform(transform, forcing, tile)
            overlap_pairs_time = transform.TEMP.overlap_pairs_time;
            overlap_pairs = transform.TEMP.overlap_pairs;
            
            del_range = isnan(overlap_pairs(:,1)) | isnan(overlap_pairs(:,2)) | abs(overlap_pairs(:,1)-overlap_pairs(:,2)) > transform.PARA.invalid_threshold;
            overlap_pairs(del_range, :) = [];
            overlap_pairs_time(del_range, :) = [];
            
            doy = floor(overlap_pairs_time - datenum(year(overlap_pairs_time), 1, 1))+1;

            slope = [];
            intercept = [];
            
%             overlap_pairs=[overlap_pairs overlap_pairs(:,1).*NaN];
            
            transform.STATVAR.slope= zeros(366,1).*NaN;
            transform.STATVAR.intercept = zeros(366,1).*NaN;
            for i=1:365
                doy_corrected = doy;
                doy_corrected(doy_corrected > i+transform.PARA.window+1) = doy_corrected(doy_corrected > i + transform.PARA.window+1) - 365;
                doy_corrected(doy_corrected < i-transform.PARA.window-1) = doy_corrected(doy_corrected < i - transform.PARA.window-1) + 365;
                range = find(doy_corrected >= i-transform.PARA.window & doy_corrected <= i+transform.PARA.window);
                
                fit_param = polyfit(overlap_pairs(range,1), overlap_pairs(range,2), 1);
                transform.STATVAR.slope(i) = fit_param(1);
                transform.STATVAR.intercept(i) = fit_param(2);
     
            end
            transform.STATVAR.slope(366) = transform.STATVAR.slope(365);
            transform.STATVAR.intercept(366) = transform.STATVAR.intercept(365);
        end
        
        
        function forcing_corrected = apply_transform(transform, forcing, tile)
            
            forcing_corrected = forcing.DATA.(transform.PARA.variable);
            
            for i=1:366
                range = find(floor(forcing.DATA.timeForcing - datenum(year(forcing.DATA.timeForcing), 1, 1))+1 == i);
                forcing_corrected(range,1) = transform.STATVAR.intercept(i) + transform.STATVAR.slope(i) .* forcing.DATA.(transform.PARA.variable)(range,1);
            end
        end
        
        
    end
end