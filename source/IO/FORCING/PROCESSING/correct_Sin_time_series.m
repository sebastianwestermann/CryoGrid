%========================================================================
% CryoGrid FORCING processing class 
%
% Authors:
% S. Westermann, December 2023
%
%========================================================================

classdef correct_Sin_time_series < process_BASE
    

    methods
        function proc = provide_PARA(proc)
            proc.PARA.number_of_quantiles = 20;
        end
        
        
        function proc = provide_CONST(proc)

        end
        
        
        function proc = provide_STATVAR(proc)
            
        end
        
        
        function proc = finalize_init(proc, tile)
            
        end
        
        
        function forcing = process(proc, forcing, tile)
            
            mo = month(forcing.DATA.timeForcing);
            ye = year(forcing.DATA.timeForcing);
            
            m= mo(1,1);
            y = ye(1,1);
            
            while datenum(y,m,1)<=datenum(ye(end), mo(end),1)
                S_in_uncorrected = forcing.CARRIER.DATA.Sin(forcing.CARRIER.DATA.timeForcing >= datenum(y, m,1) & forcing.CARRIER.DATA.timeForcing < datenum(y, m+1, 1), 1);
                S_in_corrected = forcing.DATA.Sin(mo == m & ye == y, 1);
                mean_Sin_corrected = mean(S_in_corrected);
                S_TOA_uncorrected = forcing.CARRIER.DATA.S_TOA(forcing.CARRIER.DATA.timeForcing >= datenum(y, m,1) & forcing.CARRIER.DATA.timeForcing < datenum(y, m+1, 1), 1);
                S_TOA_corrected = forcing.DATA.S_TOA(mo == m & ye == y, 1);
                
                kd_target_change = (mean(S_in_corrected) ./ mean(S_TOA_corrected)) ./ (mean(S_in_uncorrected) ./ mean(S_TOA_uncorrected));
                
                %compute daily averages
                %                 S_TOA_av = [];
                %                 S_in_av = [];
                %                 for i=datenum(y,m,1):datenum(y,m+1,1)-1
                %                     S_TOA_av = [S_TOA_av; mean(forcing.CARRIER.DATA.S_TOA(floor(forcing.CARRIER.DATA.timeForcing)==i))];
                %                     S_in_av = [S_in_av; mean(forcing.CARRIER.DATA.Sin(floor(forcing.CARRIER.DATA.timeForcing)==i))];
                %                 end
                %                 kd = S_in_av./S_TOA_av;
                kd = S_in_uncorrected ./ S_TOA_uncorrected;
                kd_NaN = isnan(kd) | isinf(kd);
                kd(kd_NaN) = [];
                
                if size(kd,1)>proc.PARA.number_of_quantiles
                    kd_quantiles = quantile(kd,proc.PARA.number_of_quantiles-1);
                    mean_quantiles = [];
                    kd_quantiles = [-Inf kd_quantiles Inf];
                    
                    quantile_number = kd .*0;
                    for i=1:proc.PARA.number_of_quantiles
                        mean_quantiles = [mean_quantiles; mean(kd(kd>kd_quantiles(1,i) & kd<=kd_quantiles(1,i+1)))];
                        quantile_number(kd>kd_quantiles(1,i) & kd<=kd_quantiles(1,i+1)) = i;
                    end
                    
                    quantile_correction = mean_quantiles.*0+1;
                    N_fixed = 1;
                    dec = 0;
                    if kd_target_change > 1
                        while ~dec
                            quantile_correction(1:end-N_fixed) = (kd_target_change .* proc.PARA.number_of_quantiles - sum(quantile_correction(end-N_fixed+1:end))) ./ (proc.PARA.number_of_quantiles - N_fixed);
                            mean_quantiles_corrected = quantile_correction .* mean_quantiles;
                            if  N_fixed < size(quantile_correction,1) && mean_quantiles_corrected(end-N_fixed) > mean_quantiles_corrected(end)
                                quantile_correction(end-N_fixed) = mean_quantiles(end)./mean_quantiles(end-N_fixed);
                                N_fixed = N_fixed +1;
                            else
                                dec=1;
                            end
                        end
                    elseif kd_target_change < 1
                        while ~dec
                            quantile_correction(N_fixed+1:end) = (kd_target_change .* proc.PARA.number_of_quantiles - sum(quantile_correction(1:N_fixed))) ./ (proc.PARA.number_of_quantiles - N_fixed);
                            mean_quantiles_corrected = quantile_correction .* mean_quantiles;
                            if  N_fixed < size(quantile_correction,1) && mean_quantiles_corrected(N_fixed+1) < mean_quantiles_corrected(1)
                                quantile_correction(N_fixed+1) = mean_quantiles(1)./mean_quantiles(N_fixed+1);
                                N_fixed = N_fixed +1;
                            else
                                dec=1;
                            end
                        end
                    end
                    %                 quantile_number = repmat(quantile_number, 1, round(1./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1))));
                    %                 quantile_number = quantile_number';
                    %                 quantile_number = quantile_number(:);
                    
                    quantile_number_all = S_in_corrected .*0;
                    quantile_number_all(~kd_NaN) = quantile_number;
                    
                    
                    for i=1:proc.PARA.number_of_quantiles
                        S_in_corrected(quantile_number_all ==i) = S_in_uncorrected(quantile_number_all ==i) ./ S_TOA_uncorrected(quantile_number_all ==i) .* quantile_correction(i) .* S_TOA_corrected(quantile_number_all ==i);
                    end
                    S_in_corrected(isnan(S_in_corrected)) = 0;
                    
                    
                    %correct residual errors
                    if mean(S_in_corrected)~=0
                        forcing.DATA.Sin(mo == m & ye == y, 1) = S_in_corrected ./ mean(S_in_corrected) .* mean_Sin_corrected;
                    else
                        forcing.DATA.Sin(mo == m & ye == y, 1) = 0;
                    end
                end
                m=m+1;
                if m==13
                    y=y+1;
                    m=1;
                end
            end

        end
        
        

        
                %-------------param file generation-----
%         function post_proc = param_file_info(post_proc)
%             post_proc = provide_PARA(post_proc);
% 
%             post_proc.PARA.STATVAR = [];
%             post_proc.PARA.class_category = 'FORCING POST_PROCESSING';
%             post_proc.PARA.options = [];
%             
%             post_proc.PARA.eliminate_fraction = [];
%             post_proc.PARA.survive_fraction = [];
%                         
%             post_proc.PARA.default_value.window_size = {7};
%             post_proc.PARA.comment.window_size = {'window size in days within which precipitation is reallocated'};
%             
%             post_proc.PARA.default_value.eliminate_fraction = {0.5};
%             post_proc.PARA.comment.eliminate_fraction = {'fraction of smallest precipitation events (= timestamps with precipitation) that is reallocated to larger events'};
%             
%             post_proc.PARA.default_value.survive_fraction = {0.5};  
%             post_proc.PARA.comment.survive_fraction = {'fraction of largest precipitation events (= timestamps with precipitation) that the small events are reallocated to'};
%             
%         end
        
    end
    
end

