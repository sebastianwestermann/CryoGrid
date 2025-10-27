%========================================================================
% CryoGrid FORCING processing class detrend_reference_and_carrier
%
% Authors:
% S. Westermann, December 2023
%
%========================================================================

classdef detrend_reference_and_carrier_monthly < process_BASE
    

    methods
        function proc = provide_PARA(proc)
%            proc.PARA.number_of_segments_carrier = []; 
            proc.PARA.detrend_interval = [];
            proc.PARA.variables = [];
            proc.PARA.relative_correction = [];
        end
        
        
        function proc = provide_CONST(proc)

        end
        
        
        function proc = provide_STATVAR(proc)
            
        end
        
        
        function proc = finalize_init(proc, tile)

        end
        
        
        function forcing = process(proc, forcing, tile)
            
            %CARRIER
            for i=1:size(proc.PARA.variables,1)
                forcing.CARRIER.STATVAR.(proc.PARA.variables{i,1}).year_points = [];
                forcing.CARRIER.STATVAR.(proc.PARA.variables{i,1}).variable_points = [];
                forcing.CARRIER.STATVAR.(proc.PARA.variables{i,1}).variable_trend = [];
%                 monthly_averages = zeros(year(forcing.CARRIER.DATA.timeForcing(end,1))-year(forcing.CARRIER.DATA.timeForcing(1,1))+1,12).*NaN;
%                 monthly_averages_detrended = zeros(year(forcing.CARRIER.DATA.timeForcing(end,1))-year(forcing.CARRIER.DATA.timeForcing(1,1))+1,12).*NaN;
%                 year_vector = [year(forcing.CARRIER.DATA.timeForcing(1,1)):year(forcing.CARRIER.DATA.timeForcing(end,1))]';
%                 y = year(forcing.CARRIER.DATA.timeForcing(1,1));
%                 m = month(forcing.CARRIER.DATA.timeForcing(1,1));
%                 year_count = 1;
%                 while datenum(y,m,1) <= datenum(year(forcing.CARRIER.DATA.timeForcing(end,1)), month(forcing.CARRIER.DATA.timeForcing(end,1)),1)
%                     monthly_averages(year_count,m) = mean(forcing.CARRIER.DATA.(proc.PARA.variables{i,1})(forcing.CARRIER.DATA.timeForcing>=datenum(y,m,1) & forcing.CARRIER.DATA.timeForcing<datenum(y,m+1,1),1));
%                     
%                     m=m+1;
%                     if m==13
%                         y=y+1;
%                         year_count = year_count+1;
%                         m=1;
%                     end
%                 end
                
                %changed Nov 
                monthly_averages = zeros(proc.PARA.detrend_interval(2)-proc.PARA.detrend_interval(1)+1,12).*NaN;
                monthly_averages_detrended = zeros(proc.PARA.detrend_interval(2)-proc.PARA.detrend_interval(1)+1,12).*NaN;
                year_vector = [proc.PARA.detrend_interval(1):proc.PARA.detrend_interval(2)]';
                y = proc.PARA.detrend_interval(1);
                m = 1;
                year_count = 1;
                while datenum(y,m,1) <= datenum(proc.PARA.detrend_interval(2), 12, 31)
                    monthly_averages(year_count,m) = mean(forcing.CARRIER.DATA.(proc.PARA.variables{i,1})(forcing.CARRIER.DATA.timeForcing>=datenum(y,m,1) & forcing.CARRIER.DATA.timeForcing<datenum(y,m+1,1),1));
                    
                    m=m+1;
                    if m==13
                        y=y+1;
                        year_count = year_count+1;
                        m=1;
                    end
                end
                
                forcing.CARRIER.STATVAR.(proc.PARA.variables{i,1}).mean_monthly_averages = mean(monthly_averages);
                %normalize to mean
                if proc.PARA.relative_correction(i,1)
                    monthly_averages = monthly_averages ./ repmat(mean(monthly_averages), size(monthly_averages,1),1);
                    monthly_averages(isnan(monthly_averages) | isinf(monthly_averages)) = 1;
                end
                
                
                year_points = [year_vector(1,1)-0.5; year_vector(end,1)+0.5];
                for m=1:12
                    fit_params = nlinfit(year_vector ,monthly_averages(:,m), @(fit_params, y_vec)piece_wise_linear(fit_params, y_vec, year_points), [mean(monthly_averages(:,m)); mean(monthly_averages(:,m))]);
                    
                    unconstrained_years = find(isnan(year_points));
                    year_points(unconstrained_years) = fit_params(1:size(unconstrained_years,1));
                    variable_points = fit_params(size(unconstrained_years,1)+1:end, 1);
                    if ~proc.PARA.relative_correction(i,1)
                        variable_trend = piece_wise_linear(fit_params, year_vector, year_points) - variable_points(1,1);
                    else
                        variable_trend = piece_wise_linear(fit_params, year_vector, year_points);
                    end
                        
                    forcing.CARRIER.STATVAR.(proc.PARA.variables{i,1}).year_points = [forcing.CARRIER.STATVAR.(proc.PARA.variables{i,1}).year_points year_points];
                    forcing.CARRIER.STATVAR.(proc.PARA.variables{i,1}).variable_points = [forcing.CARRIER.STATVAR.(proc.PARA.variables{i,1}).variable_points variable_points];
                    forcing.CARRIER.STATVAR.(proc.PARA.variables{i,1}).variable_trend = [forcing.CARRIER.STATVAR.(proc.PARA.variables{i,1}).variable_trend variable_trend];
                    if ~proc.PARA.relative_correction(i,1)
                        monthly_averages_detrended(:,m) = monthly_averages(:,m) - variable_trend;
                    else
                        monthly_averages_detrended(:,m) = monthly_averages(:,m) ./ variable_trend;
                    end
                end
                forcing.CARRIER.STATVAR.(proc.PARA.variables{i,1}).monthly_std = nanstd(monthly_averages_detrended);
                forcing.CARRIER.STATVAR.(proc.PARA.variables{i,1}).relative = proc.PARA.relative_correction(i,1);
                
                y = year(forcing.CARRIER.DATA.timeForcing(1,1));
                m = month(forcing.CARRIER.DATA.timeForcing(1,1));
                year_count = 1;
                while datenum(y,m,1) <= datenum(proc.PARA.detrend_interval(2), 12, 31)
%                                     while datenum(y,m,1) <= datenum(year(forcing.CARRIER.DATA.timeForcing(end,1)), month(forcing.CARRIER.DATA.timeForcing(end,1)),1)
                    if ~proc.PARA.relative_correction(i,1)
                    forcing.CARRIER.DATA.(proc.PARA.variables{i,1})(forcing.CARRIER.DATA.timeForcing>=datenum(y,m,1) & forcing.CARRIER.DATA.timeForcing<datenum(y,m+1,1),1) = ...
                        forcing.CARRIER.DATA.(proc.PARA.variables{i,1})(forcing.CARRIER.DATA.timeForcing>=datenum(y,m,1) & forcing.CARRIER.DATA.timeForcing<datenum(y,m+1,1),1) - forcing.CARRIER.STATVAR.(proc.PARA.variables{i,1}).variable_trend(year_count,m);
                    else
                       forcing.CARRIER.DATA.(proc.PARA.variables{i,1})(forcing.CARRIER.DATA.timeForcing>=datenum(y,m,1) & forcing.CARRIER.DATA.timeForcing<datenum(y,m+1,1),1) = ...
                        forcing.CARRIER.DATA.(proc.PARA.variables{i,1})(forcing.CARRIER.DATA.timeForcing>=datenum(y,m,1) & forcing.CARRIER.DATA.timeForcing<datenum(y,m+1,1),1) ./ forcing.CARRIER.STATVAR.(proc.PARA.variables{i,1}).variable_trend(year_count,m);
                    end
                    m=m+1;
                    if m==13
                        y=y+1;
                        year_count = year_count+1;
                        m=1;
                    end
                end
            end
            
            %REFERENCE
            for i=1:size(proc.PARA.variables,1)
                forcing.REFERENCE.STATVAR.(proc.PARA.variables{i,1}).year_points = [];
                forcing.REFERENCE.STATVAR.(proc.PARA.variables{i,1}).variable_points = [];
                forcing.REFERENCE.STATVAR.(proc.PARA.variables{i,1}).variable_trend = [];
                monthly_averages = zeros(year(forcing.REFERENCE.DATA.timeForcing(end,1))-year(forcing.REFERENCE.DATA.timeForcing(1,1))+1,12).*NaN;
                monthly_averages_detrended = zeros(year(forcing.REFERENCE.DATA.timeForcing(end,1))-year(forcing.REFERENCE.DATA.timeForcing(1,1))+1,12).*NaN;
                year_vector = [year(forcing.REFERENCE.DATA.timeForcing(1,1)):year(forcing.REFERENCE.DATA.timeForcing(end,1))]';
                y = year(forcing.REFERENCE.DATA.timeForcing(1,1));
                m = month(forcing.REFERENCE.DATA.timeForcing(1,1));
                year_count = 1;
                while datenum(y,m,1) <= datenum(year(forcing.REFERENCE.DATA.timeForcing(end,1)), month(forcing.REFERENCE.DATA.timeForcing(end,1)),1)
                    monthly_averages(year_count,m) = mean(forcing.REFERENCE.DATA.(proc.PARA.variables{i,1})(forcing.REFERENCE.DATA.timeForcing>=datenum(y,m,1) & forcing.REFERENCE.DATA.timeForcing<datenum(y,m+1,1),1));
                    
                    m=m+1;
                    if m==13
                        y=y+1;
                        year_count = year_count+1;
                        m=1;
                    end
                end
                
                range_CARRIER = year_vector(:,1) >= year(forcing.CARRIER.DATA.timeForcing(1,1)) & year_vector(:,end) <= year(forcing.CARRIER.DATA.timeForcing(end,1));

                forcing.REFERENCE.STATVAR.(proc.PARA.variables{i,1}).mean_monthly_averages = mean(monthly_averages(range_CARRIER,:));
                %normalize to mean
                if proc.PARA.relative_correction(i,1)
                    monthly_averages = monthly_averages ./ repmat(mean(monthly_averages(range_CARRIER,:)), size(monthly_averages,1),1);
                    monthly_averages(isnan(monthly_averages) | isinf(monthly_averages)) = 0;
                end
                
                
                for m=1:12
                    year_points = [year_vector(1,1)-0.5; forcing.CARRIER.STATVAR.(proc.PARA.variables{i,1}).year_points(:,m); year_vector(end,1)+0.5];

                    fit_params = nlinfit(year_vector ,monthly_averages(:,m), @(fit_params, y_vec)piece_wise_linear(fit_params, y_vec, year_points), repmat(mean(monthly_averages(:,m)), size(year_points,1), 1));
                    
                    unconstrained_years = find(isnan(year_points));
                    year_points(unconstrained_years) = fit_params(1:size(unconstrained_years,1));
                    variable_points = fit_params(size(unconstrained_years,1)+1:end, 1);
                    
                    if ~proc.PARA.relative_correction(i,1)
                        variable_trend = piece_wise_linear(fit_params, year_vector, year_points) - variable_points(1,1);
                    else
                        variable_trend = piece_wise_linear(fit_params, year_vector, year_points);
                    end
                    forcing.REFERENCE.STATVAR.(proc.PARA.variables{i,1}).year_points = [forcing.REFERENCE.STATVAR.(proc.PARA.variables{i,1}).year_points year_points];
                    forcing.REFERENCE.STATVAR.(proc.PARA.variables{i,1}).variable_points = [forcing.REFERENCE.STATVAR.(proc.PARA.variables{i,1}).variable_points variable_points];
                    forcing.REFERENCE.STATVAR.(proc.PARA.variables{i,1}).variable_trend = [forcing.REFERENCE.STATVAR.(proc.PARA.variables{i,1}).variable_trend variable_trend];
                    if ~proc.PARA.relative_correction(i,1)
                        monthly_averages_detrended(:,m) = monthly_averages(:,m) - variable_trend;
                    else
                        monthly_averages_detrended(:,m) = monthly_averages(:,m) ./ variable_trend;
                        monthly_averages_detrended(monthly_averages(:,m)==0,m) = 0;
                    end
                end
                forcing.REFERENCE.STATVAR.(proc.PARA.variables{i,1}).monthly_std = nanstd(monthly_averages_detrended(range_CARRIER,:));
                       
                
                monthly_difference = monthly_averages_detrended - repmat(mean(monthly_averages_detrended,1), size(monthly_averages_detrended,1),1);
                monthly_difference = monthly_difference .* (repmat(forcing.CARRIER.STATVAR.(proc.PARA.variables{i,1}).monthly_std./forcing.REFERENCE.STATVAR.(proc.PARA.variables{i,1}).monthly_std, size(monthly_averages_detrended,1),1)-1);
                
                y = year(forcing.REFERENCE.DATA.timeForcing(1,1));
                m = month(forcing.REFERENCE.DATA.timeForcing(1,1));
                year_count = 1;
                while datenum(y,m,1) <= datenum(year(forcing.REFERENCE.DATA.timeForcing(end,1)), month(forcing.REFERENCE.DATA.timeForcing(end,1)),1)
                    if ~proc.PARA.relative_correction(i,1)
                        forcing.REFERENCE.DATA.(proc.PARA.variables{i,1})(forcing.REFERENCE.DATA.timeForcing>=datenum(y,m,1) & forcing.REFERENCE.DATA.timeForcing<datenum(y,m+1,1),1) = ...
                            forcing.REFERENCE.DATA.(proc.PARA.variables{i,1})(forcing.REFERENCE.DATA.timeForcing>=datenum(y,m,1) & forcing.REFERENCE.DATA.timeForcing<datenum(y,m+1,1),1) ...
                            - forcing.REFERENCE.STATVAR.(proc.PARA.variables{i,1}).variable_trend(year_count,m) + monthly_difference(year_count,m);
                    else
                        forcing.REFERENCE.DATA.(proc.PARA.variables{i,1})(forcing.REFERENCE.DATA.timeForcing>=datenum(y,m,1) & forcing.REFERENCE.DATA.timeForcing<datenum(y,m+1,1),1) = ...
                            forcing.REFERENCE.DATA.(proc.PARA.variables{i,1})(forcing.REFERENCE.DATA.timeForcing>=datenum(y,m,1) & forcing.REFERENCE.DATA.timeForcing<datenum(y,m+1,1),1) ./ ...
                            forcing.REFERENCE.STATVAR.(proc.PARA.variables{i,1}).variable_trend(year_count,m);
                    end
                    m=m+1;
                    if m==13
                        y=y+1;
                        year_count = year_count+1;
                        m=1;
                    end
                end
                %adapt std for  relative correction
                if proc.PARA.relative_correction(i,1)
                    month_ref = month(forcing.REFERENCE.DATA.timeForcing);
                    for m=1:12
                        mean_ref = mean(forcing.REFERENCE.DATA.(proc.PARA.variables{i,1})(month_ref==m,1));
                        data =forcing.REFERENCE.DATA.(proc.PARA.variables{i,1})(month_ref==m,1);
                        deviation_from_mean = (data - mean_ref) ./mean_ref;
                        deviation_from_mean(isnan(deviation_from_mean) | isinf(deviation_from_mean)) = 1;
                        smaller_than_zero = deviation_from_mean <0;
                        deviation_from_mean(smaller_than_zero) = mean_ref ./ data(smaller_than_zero); 
                        deviation_from_mean = deviation_from_mean .* forcing.CARRIER.STATVAR.(proc.PARA.variables{i,1}).monthly_std(1,m)./forcing.REFERENCE.STATVAR.(proc.PARA.variables{i,1}).monthly_std(1,m);
                        deviation_from_mean(isnan(deviation_from_mean) | isinf(deviation_from_mean)) = 1;
                        
                        data_revised  =  deviation_from_mean .* mean_ref + mean_ref;
                        data_revised(smaller_than_zero) = mean_ref ./ deviation_from_mean(smaller_than_zero);
                        
                        forcing.REFERENCE.DATA.(proc.PARA.variables{i,1})(month_ref==m,1) = data_revised; %deviation_from_mean .* mean_ref + mean_ref;
                    end
                     
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

