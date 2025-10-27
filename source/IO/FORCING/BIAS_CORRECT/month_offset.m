%TRANSFORM class month_offset

% S. Westermann Dec 2023


classdef month_offset < matlab.mixin.Copyable
    
    properties
        PARA
        STATVAR
        TEMP
        CONST
    end
    
    methods
        
        function transform = provide_PARA(transform)

           transform.PARA.relative_correction = 0;
           transform.PARA.trend_reference_to_carrier = 0;
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
            
            
            m = month(overlap_pairs_time);

            transform.STATVAR.average_reference = zeros(12, 1).*NaN;
            transform.STATVAR.average_carrier = zeros(12, 1).*NaN;

            for i=1:12
                transform.STATVAR.average_reference(i,1) = mean(overlap_pairs(find(m==i),2));
                transform.STATVAR.average_carrier(i,1) = mean(overlap_pairs(find(m==i),1));
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
                    average_carrier = mean(forcing.CARRIER.DATA.(transform.PARA.variable)(range,1));
                    average_reference = mean(forcing.REFERENCE.DATA.(transform.PARA.variable)(forcing.REFERENCE.DATA.timeForcing>=datenum(y,m,1) & forcing.REFERENCE.DATA.timeForcing < datenum(y,m+1,1),1));
                    
                    if transform.PARA.relative_correction == 1
                        forcing_corrected(range,1) = forcing_corrected(range,1) ./ (transform.STATVAR.average_reference(m,1) ./ average_reference) .* (transform.STATVAR.average_carrier(m,1) ./ average_carrier);
                    else
                        forcing_corrected(range,1) = forcing_corrected(range,1) - (transform.STATVAR.average_reference(m,1) - average_reference) + (transform.STATVAR.average_carrier(m,1) - average_carrier);
                    end
                else
                    if transform.PARA.relative_correction == 1
                        forcing_corrected(range,1) = forcing_corrected(range,1) .* transform.STATVAR.average_reference(m,1) ./ transform.STATVAR.average_carrier(m,1);
                    else
                        forcing_corrected(range,1) = forcing_corrected(range,1) + transform.STATVAR.average_reference(m,1) - transform.STATVAR.average_carrier(m,1);
                    end
                    
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