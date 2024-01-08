%TRANSFORM class month_offset -DISCONTINUED

% S. Westermann Dec 2023


classdef month_offset_Sin < matlab.mixin.Copyable
    
    properties
        PARA
        STATVAR
        TEMP
        CONST
    end
    
    methods
        
        function transform = provide_PARA(transform)
           transform.PARA.trend_reference_to_carrier = 0;
           transform.PARA.number_of_quantiles = 10;
        end
        
        function transform = provide_CONST(transform)
            
        end
        
        function transform = provide_STATVAR(transform)

        end
        
        function transform = finalize_init(transform, tile)

        end
        
        function transform = fit_transform(transform, forcing, tile)
            overlap_pairs_time = transform.TEMP.overlap_pairs_time;
            overlap_pairs = transform.TEMP.overlap_pairs;
            
            %add S_TOA overlap paris 
            S_TOA_carrier = [];
            S_TOA_reference = [];
            for i=1:size(overlap_pairs_time,1)
                S_TOA_carrier = [S_TOA_carrier; mean(forcing.CARRIER.DATA.S_TOA(forcing.CARRIER.DATA.timeForcing>=datenum(year(overlap_pairs_time(i,1)), month(overlap_pairs_time(i,1)),1) & ...
                  forcing.CARRIER.DATA.timeForcing < datenum(year(overlap_pairs_time(i,1)), month(overlap_pairs_time(i,1))+1, 1) ,1))];
                S_TOA_reference = [S_TOA_reference; mean(forcing.REFERENCE.DATA.S_TOA(forcing.REFERENCE.DATA.timeForcing>=datenum(year(overlap_pairs_time(i,1)), month(overlap_pairs_time(i,1)),1) & ...
                  forcing.REFERENCE.DATA.timeForcing < datenum(year(overlap_pairs_time(i,1)), month(overlap_pairs_time(i,1))+1, 1) ,1))];
            end
            
            transform.TEMP.overlap_pairs_S_TOA = [S_TOA_carrier S_TOA_reference];
            transform.TEMP.overlap_pairs_kd = transform.TEMP.overlap_pairs ./ transform.TEMP.overlap_pairs_S_TOA;

            m = month(overlap_pairs_time);

            transform.STATVAR.average_reference = zeros(12, 1).*NaN;
            transform.STATVAR.average_carrier = zeros(12, 1).*NaN;

            for i=1:12
                transform.STATVAR.average_reference_Sin(i,1) = mean(overlap_pairs(find(m==i),2));
                transform.STATVAR.average_carrier_Sin(i,1) = mean(overlap_pairs(find(m==i),1));
                transform.STATVAR.average_reference_S_TOA(i,1) = mean(transform.TEMP.overlap_pairs_S_TOA(find(m==i),2));
                transform.STATVAR.average_carrier_S_TOA(i,1) = mean(transform.TEMP.overlap_pairs_S_TOA(find(m==i),1));
                transform.STATVAR.average_reference_kd(i,1) = nanmean(transform.TEMP.overlap_pairs_kd(find(m==i),2));
                transform.STATVAR.average_carrier_kd(i,1) = nanmean(transform.TEMP.overlap_pairs_kd(find(m==i),1));
            end
            
            
            

        end
        
        
        function forcing_corrected = apply_transform(transform, forcing, tile)
            
            forcing_corrected = forcing.DATA.(transform.PARA.variable);
            mo = month(forcing.CARRIER.DATA.timeForcing);
            ye = year(forcing.CARRIER.DATA.timeForcing);
            
            m= mo(1,1);
            y = ye(1,1);
            
            while datenum(y,m,1)<=datenum(ye(end), mo(end),1)
                range = find(mo == m & ye == y);
                if  transform.PARA.trend_reference_to_carrier
                    average_carrier_Sin = mean(forcing.CARRIER.DATA.(transform.PARA.variable)(range,1));
                    average_reference_Sin = mean(forcing.REFERENCE.DATA.(transform.PARA.variable)(forcing.REFERENCE.DATA.timeForcing>=datenum(y,m,1) & forcing.REFERENCE.DATA.timeForcing < datenum(y,m+1,1),1));
                    average_carrier_S_TOA = mean(forcing.CARRIER.DATA.S_TOA(range,1));
                    average_reference_S_TOA = mean(forcing.REFERENCE.DATA.S_TOA(forcing.REFERENCE.DATA.timeForcing>=datenum(y,m,1) & forcing.REFERENCE.DATA.timeForcing < datenum(y,m+1,1),1));
                    average_carrier_kd = average_carrier_Sin ./ average_carrier_S_TOA;
                    average_reference_kd = average_reference_Sin ./ average_reference_S_TOA;
                    
                    relative_change_S_TOA = 1 ./ (transform.STATVAR.average_reference_S_TOA(m,1) ./ average_reference_S_TOA) .* (transform.STATVAR.average_carrier_S_TOA(m,1) ./ average_carrier_S_TOA);
                    relative_change_kd = 1 ./ (transform.STATVAR.average_reference_kd(m,1) ./ average_reference_kd) .* (transform.STATVAR.average_carrier_kd(m,1) ./ average_carrier_kd);
                    
                    forcing_corrected(range,1) = forcing_corrected(range,1) ./ (transform.STATVAR.average_reference(m,1) ./ average_reference) .* (transform.STATVAR.average_carrier(m,1) ./ average_carrier);
                else
                    
                    forcing_corrected(range,1) = forcing_corrected(range,1) .* transform.STATVAR.average_reference(m,1) ./ transform.STATVAR.average_carrier(m,1);
                    
                end

                m=m+1;
                if m==13
                    y=y+1;
                    m=1;
                end

            end
        end
        
        
    end
end