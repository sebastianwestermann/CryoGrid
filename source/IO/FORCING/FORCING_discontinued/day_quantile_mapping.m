classdef day_quantile_mapping < matlab.mixin.Copyable
    
    properties
        PARA
        STATVAR
        CONST
    end
    
    methods
        
        function bias_correct = provide_PARA(bias_correct)
           bias_correct.PARA.window = []; %in days
           bias_correct.PARA.invalid_threshold = [];
           bias_correct.PARA.number_of_bins = [];
        end
        
        function bias_correct = provide_CONST(bias_correct)
            
        end
        
        function bias_correct = provide_STATVAR(bias_correct)
            
        end
        
        function bias_correct = finalize_init(bias_correct, tile)

        end
        
        function forcing_corrected = bias_correct(bias_correct, reference_class, tile)
            overlap_pairs_time = bias_correct.STATVAR.overlap_pairs_time;
            overlap_pairs = bias_correct.STATVAR.overlap_pairs;
            
            overlap_pairs(isnan(overlap_pairs(:,1)) | isnan(overlap_pairs(:,2)) | abs(overlap_pairs(:,1)-overlap_pairs(:,2)) > bias_correct.PARA.invalid_threshold, :) = [];
            doy = floor(bias_correct.STATVAR.overlap_pairs_time - datenum(year(bias_correct.STATVAR.overlap_pairs_time), 1, 1))+1;

            forcing_corrected = reference_class.DATA.(bias_correct.PARA.variable);
            
            for i=1:365
                doy_corrected = doy;
                doy_corrected(doy_corrected > i+bias_correct.PARA.window+1) = doy_corrected(doy_corrected > i + bias_correct.PARA.window+1) - 365;
                doy_corrected(doy_corrected < i-bias_correct.PARA.window-1) = doy_corrected(doy_corrected < i - bias_correct.PARA.window-1) + 365;
                range = find(doy_corrected >= i-bias_correct.PARA.window & doy_corrected <= i+bias_correct.PARA.window);
                
                quantile_factor = [];
                carrier_sorted = sort(overlap_pairs(range,1);
                reference_sorted = sort(overlap_pairs(range,2);
                for j=1:bias_correct.PARA.number_of_bins
                    range = round((j-1) .* size(carrier_sorted,1)./bias_correct.PARA.number_of_bins +1):round(j.* size(carrier_sorted,1)./bias_correct.PARA.number_of_bins);
                    quantile_factor = [quantile_factor; [mean(carrier_sorted(range,1))./mean(reference_sorted(range,1)) reference_sorted(range(1),1)]];
                    if isinf(quantile_factor(end,1)) || isnan(quantile_factor(end,1))
                        quantile_factor(end,1)=1;
                    end
                end
                
                
                range = find(floor(reference_class.DATA.timeForcing - datenum(year(reference_class.DATA.timeForcing), 1, 1))+1 == i);
                for j=1:length(range)
                   factor = quantile_factor(:,2)(find(quantile_factor(:,2) < forcing_corrected(range(j),1), 1, 'last'),1)
                   forcing_corrected(range(j),1) = forcing_corrected(range(j),1) .* factor;
                end                
            end
            i = 366;
            range = find(floor(reference_class.DATA.timeForcing - datenum(year(reference_class.DATA.timeForcing), 1, 1))+1 == i); %same as i=365
            for j=1:length(range)
                factor = quantile_factor(:,2)(find(quantile_factor(:,2) < forcing_corrected(range(j),1), 1, 'last'),1)
                forcing_corrected(range(j),1) = forcing_corrected(range(j),1) .* factor;
            end
        end
        
        
    end
end