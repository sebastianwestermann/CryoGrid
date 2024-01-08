%TRANSFORM class doy_linear_fit

% S. Westermann Dec 2023


classdef doy_quantile_mapping < matlab.mixin.Copyable
    
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
           transform.PARA.number_of_quantiles = [];
           transform.PARA.recompute_quantiles = []; %1:quantiles are computed for the target forcing time series; 0: quantiles from overlap_paris are used
           transform.PARA.relative_correction = 0;
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

            transform.STATVAR.quantile_correction = zeros(366, transform.PARA.number_of_quantiles).*NaN;
            transform.STATVAR.quantiles = zeros(366, transform.PARA.number_of_quantiles-1).*NaN;

            for i=1:365
                doy_corrected = doy;
                doy_corrected(doy_corrected > i+transform.PARA.window+1) = doy_corrected(doy_corrected > i + transform.PARA.window+1) - 365;
                doy_corrected(doy_corrected < i-transform.PARA.window-1) = doy_corrected(doy_corrected < i - transform.PARA.window-1) + 365;
                range = find(doy_corrected >= i-transform.PARA.window & doy_corrected <= i+transform.PARA.window);
                
                %1: carrier
                overlap_pairs_doy = overlap_pairs(range,:);
                quantiles_carrier = [-Inf quantile(overlap_pairs_doy(:,1), transform.PARA.number_of_quantiles-1) Inf];
                quantiles_reference = [-Inf quantile(overlap_pairs_doy(:,2), transform.PARA.number_of_quantiles-1) Inf];
                
                transform.STATVAR.quantiles(i,:) = quantiles_carrier(1,2:end-1);
                
                quantiles_mean_carrier = zeros(1, transform.PARA.number_of_quantiles).*NaN;
                quantiles_mean_reference = zeros(1, transform.PARA.number_of_quantiles).*NaN;
                
                for j=1:transform.PARA.number_of_quantiles
                    quantiles_mean_carrier(1,j) = nanmean(overlap_pairs_doy(overlap_pairs_doy(:,1)>=quantiles_carrier(1,j) & overlap_pairs_doy(:,1)<quantiles_carrier(1,j+1),1));
                    quantiles_mean_reference(1,j) = nanmean(overlap_pairs_doy(overlap_pairs_doy(:,2)>=quantiles_reference(1,j) & overlap_pairs_doy(:,2)<quantiles_reference(1,j+1),2));                    
                end
                if transform.PARA.relative_correction == 1
                    correction = quantiles_mean_reference ./ quantiles_mean_carrier;
                    correction(isnan(correction)) = 1;
                    transform.STATVAR.quantile_correction(i,:) = correction;
                else
                    correction = quantiles_mean_reference - quantiles_mean_carrier;                     
                    correction(isnan(correction)) = 0;
                    transform.STATVAR.quantile_correction(i,:) = correction;
                end
            end
            transform.STATVAR.quantile_correction(366,:) = transform.STATVAR.quantile_correction(365,:);
            transform.STATVAR.quantiles(366,:) = transform.STATVAR.quantiles(365,:);
        end
        
        
        function forcing_corrected = apply_transform(transform, forcing, tile)
            
            forcing_corrected = forcing.DATA.(transform.PARA.variable);
            doy = floor(forcing.DATA.timeForcing - datenum(year(forcing.DATA.timeForcing), 1, 1))+1;
            for i=1:366
                if transform.PARA.recompute_quantiles == 1
                    doy_corrected = doy;
                    doy_corrected(doy_corrected > i+transform.PARA.window+1) = doy_corrected(doy_corrected > i + transform.PARA.window+1) - 365;
                    doy_corrected(doy_corrected < i-transform.PARA.window-1) = doy_corrected(doy_corrected < i - transform.PARA.window-1) + 365;
                    range = find(doy_corrected >= i-transform.PARA.window & doy_corrected <= i+transform.PARA.window);
                    quantiles_carrier = [quantile(forcing.DATA.(transform.PARA.variable)(range,1), transform.PARA.number_of_quantiles-1) Inf];
                else
                    quantiles_carrier = [transform.STATVAR.quantiles(i,:) Inf];
                end
                range = find(floor(forcing.DATA.timeForcing - datenum(year(forcing.DATA.timeForcing), 1, 1))+1 == i);
                
                quantile_score = 1;
                for j = 1:transform.PARA.number_of_quantiles
                    quantile_score = quantile_score + double(forcing_corrected(range,1)>=quantiles_carrier(1,j));
                end
                correction = quantile_score.*0;
                for j = 1:transform.PARA.number_of_quantiles
                    correction(quantile_score == j,1) = transform.STATVAR.quantile_correction(i,j);
                end
                
                if transform.PARA.relative_correction == 1
                    forcing_corrected(range,1) = forcing_corrected(range,1) .* correction;
                else
                    forcing_corrected(range,1) = forcing_corrected(range,1) + correction;                    
                end
                
            end
        end
        
        
    end
end