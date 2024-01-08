classdef doy_linear_fit_zero_intercept < matlab.mixin.Copyable
    
    properties
        PARA
        STATVAR
        CONST
    end
    
    methods
        
        function bias_correct = provide_PARA(bias_correct)
           bias_correct.PARA.window = []; %in days
           bias_correct.PARA.invalid_threshold = [];
        end
        
        function bias_correct = provide_CONST(bias_correct)

        end
        
        function bias_correct = provide_STATVAR(bias_correct)
            
        end
        
        function bias_correct = finalize_init(bias_correct, tile)
            if isempty(bias_correct.PARA.invalid_threshold) || isnan(bias_correct.PARA.invalid_threshold)
                 bias_correct.PARA.invalid_threshold = Inf;
            end
        end
        
        function [forcing_corrected, overlap_pairs] = bias_correct(bias_correct, forcing, carrier_class, tile)
            overlap_pairs_time = bias_correct.STATVAR.overlap_pairs_time;
            overlap_pairs = bias_correct.STATVAR.overlap_pairs;
            
            del_range = isnan(overlap_pairs(:,1)) | isnan(overlap_pairs(:,2)) | abs(overlap_pairs(:,1)-overlap_pairs(:,2)) > bias_correct.PARA.invalid_threshold;
            overlap_pairs(del_range, :) = [];
            overlap_pairs_time(del_range, :) = [];
            doy = floor(overlap_pairs_time - datenum(year(overlap_pairs_time), 1, 1))+1;


            forcing_corrected = carrier_class.DATA.(bias_correct.PARA.variable);
            
           overlap_pairs=[overlap_pairs overlap_pairs(:,1).*NaN];
            
            warning off all
            
            for i=1:365
                doy_corrected = doy;
                doy_corrected(doy_corrected > i+bias_correct.PARA.window+1) = doy_corrected(doy_corrected > i + bias_correct.PARA.window+1) - 365;
                doy_corrected(doy_corrected < i-bias_correct.PARA.window-1) = doy_corrected(doy_corrected < i - bias_correct.PARA.window-1) + 365;
                range = find(doy_corrected >= i-bias_correct.PARA.window & doy_corrected <= i+bias_correct.PARA.window);
                
                P = nlinfit(overlap_pairs(range,1), overlap_pairs(range,2), @(a,x) a.*x, 1);
                P(isnan(P))=1;
                P = min(max(P, 0.5), 2);

               overlap_pairs(range,3) = P.*overlap_pairs(range,1);

                
                range = find(floor(carrier_class.DATA.timeForcing - datenum(year(carrier_class.DATA.timeForcing), 1, 1))+1 == i);
                forcing_corrected(range,1) = P .* forcing_corrected(range,1);
            end
            i = 366;
            range = find(floor(carrier_class.DATA.timeForcing - datenum(year(carrier_class.DATA.timeForcing), 1, 1))+1 == i); %same as i=365
            forcing_corrected(range,1) =  P .* forcing_corrected(range,1);
        end
        
        
    end
end