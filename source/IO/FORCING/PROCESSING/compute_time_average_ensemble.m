%========================================================================
% CryoGrid FORCING post-processing class compute_time_average_ensemble
%
% averaging of forcing data for SEB calculations,  including ensemble generation 
% all variables stay the same

% Authors:
% S. Westermann, December 2022
%
%========================================================================

%ATTENTION: CLASS; changed, must be adapted in the ENSEMBLE  generation/DA classes for ESA CCI 

classdef compute_time_average_ensemble < FORCING_base

    properties
        
    end
    
    methods
        function proc = provide_PARA(proc)
            
            proc.PARA.averaging_period = [];  %in days
            proc.PARA.all_snow_T = [];
            proc.PARA.all_rain_T = [];
            
            %these can be written by ensemble class
            proc.PARA.ensemble_size = 1;
            proc.PARA.absolute_change_Tair = 0;
            proc.PARA.snow_fraction = 1;
            proc.PARA.rain_fraction = 1;
            proc.PARA.relative_change_Sin = 1;                 
        end
        
        
        function proc = provide_CONST(proc)
            proc.CONST.L_f = []; 
            proc.CONST.sigma = [];
            proc.CONST.day_sec = [];
            proc.CONST.Tmfw = [];
        end
        
        
        function proc = provide_STATVAR(proc)
            
        end
        
        
        function proc = finalize_init(proc, tile)
            
            if size(proc.PARA.absolute_change_Tair,2)==1 
                proc.PARA.absolute_change_Tair = repmat(proc.PARA.absolute_change_Tair,1, proc.PARA.ensemble_size);
            end   
            if size(proc.PARA.snow_fraction,2)==1 
                proc.PARA.snow_fraction = repmat(proc.PARA.snow_fraction,1, proc.PARA.ensemble_size);
            end
            if size(proc.PARA.rain_fraction,2)==1
                proc.PARA.rain_fraction = repmat(proc.PARA.rain_fraction,1, proc.PARA.ensemble_size);
            end
            if size(proc.PARA.relative_change_Sin,2)==1
                proc.PARA.relative_change_Sin = repmat(proc.PARA.relative_change_Sin,1, proc.PARA.ensemble_size);
            end
        end
        
        
        function forcing = process(proc, forcing, tile)

            if ~isfield(forcing.DATA, 'full_res_data') %CHANGED SW Jan 2024, should make the function "post_process_adjust_forcing" redundant
                data_full = forcing.DATA;
            else
                data_full = forcing.DATA.full_res_data;
            end
            
            forcing.DATA = [];
            forcing.DATA.snowfall = [];
            forcing.DATA.rainfall = [];
            forcing.DATA.Tair = [];
            forcing.DATA.wind = [];
            forcing.DATA.q = [];
            forcing.DATA.p = [];
            forcing.DATA.Lin = [];
            forcing.DATA.Sin = [];
            forcing.DATA.timeForcing = [];
            
            for i = data_full.timeForcing(1,1):proc.PARA.averaging_period:data_full.timeForcing(end,1)-proc.PARA.averaging_period
                range = find(data_full.timeForcing>=i & data_full.timeForcing < min(data_full.timeForcing(end,1), i + proc.PARA.averaging_period));
                forcing.DATA.timeForcing = [forcing.DATA.timeForcing; mean(data_full.timeForcing(range,1))];
                forcing.DATA.Tair = [forcing.DATA.Tair; mean(data_full.Tair(range,1)) + proc.PARA.absolute_change_Tair];
                forcing.DATA.wind = [forcing.DATA.wind; mean(data_full.wind(range,1))];
                forcing.DATA.q = [forcing.DATA.q; mean(data_full.q(range,1))];
                forcing.DATA.p = [forcing.DATA.p; mean(data_full.p(range,1))];

                sky_emissivity = data_full.Lin(range,1) ./ (data_full.Tair(range,1)+273.15).^4 ./ proc.CONST.sigma;
                forcing.DATA.Lin = [forcing.DATA.Lin; mean(sky_emissivity .* proc.CONST.sigma .* (data_full.Tair(range,1) + 273.15 + proc.PARA.absolute_change_Tair).^4)];
                forcing.DATA.Sin = [forcing.DATA.Sin; mean(data_full.Sin(range,1) .*  (1+proc.PARA.relative_change_Sin))];
                
                sf = 0;
                rf = 0;
                for j=1:size(range,1)
                    precip = data_full.snowfall(range(j),1) + data_full.rainfall(range(j),1);
                    factor = max(0, min(1, (data_full.Tair(range(j),1) + proc.PARA.absolute_change_Tair - proc.PARA.all_snow_T) ./ max(1e-12, (proc.PARA.all_rain_T - proc.PARA.all_snow_T))));
                    sf = sf + precip.*(1 - factor);
                    rf = rf + precip.*factor;
                end
                forcing.DATA.snowfall = [forcing.DATA.snowfall; sf./size(range,1) .* proc.PARA.snow_fraction];
                forcing.DATA.rainfall = [forcing.DATA.rainfall; rf./size(range,1) .* proc.PARA.rain_fraction];
                                
            end
            
            if size(forcing.PARA.heatFlux_lb,2) == 1
                tile.PARA.geothermal = repmat(forcing.PARA.heatFlux_lb, 1, proc.PARA.ensemble_size);
            end
            
            %overwrite target variables in TEMP in FORCING
            forcing.TEMP = [];
            forcing.TEMP.Lin = 0;
            forcing.TEMP.Sin = 0;
            forcing.TEMP.q = 0;
            forcing.TEMP.p = 0;
            forcing.TEMP.Tair = 0;
            forcing.TEMP.wind = 0;
            forcing.TEMP.rainfall = 0;
            forcing.TEMP.snowfall = 0;
            
            forcing.DATA.full_res_data = data_full;
        end
        
        
%         function forcing = post_process_adjust_forcing(proc, forcing, tile)
% 
%             data_full = forcing.DATA.full_res_data;
%             forcing.DATA = [];
%             forcing.DATA.snowfall = [];
%             forcing.DATA.rainfall = [];
%             forcing.DATA.Tair = [];
%             forcing.DATA.wind = [];
%             forcing.DATA.q = [];
%             forcing.DATA.p = [];
%             forcing.DATA.Lin = [];
%             forcing.DATA.Sin = [];
%             forcing.DATA.timeForcing = [];
%             
%             for i = data_full.timeForcing(1,1):proc.PARA.averaging_period:data_full.timeForcing(end,1)-proc.PARA.averaging_period
%                 range = find(data_full.timeForcing>=i & data_full.timeForcing < min(data_full.timeForcing(end,1), i + proc.PARA.averaging_period));
%                 forcing.DATA.timeForcing = [forcing.DATA.timeForcing; mean(data_full.timeForcing(range,1))];
%                 forcing.DATA.Tair = [forcing.DATA.Tair; mean(data_full.Tair(range,1)) + proc.PARA.absolute_change_Tair];
%                 forcing.DATA.wind = [forcing.DATA.wind; mean(data_full.wind(range,1))];
%                 forcing.DATA.q = [forcing.DATA.q; mean(data_full.q(range,1))];
%                 forcing.DATA.p = [forcing.DATA.p; mean(data_full.p(range,1))];
% 
%                 sky_emissivity = data_full.Lin(range,1) ./ (data_full.Tair(range,1)+273.15).^4 ./ proc.CONST.sigma;
%                 forcing.DATA.Lin = [forcing.DATA.Lin; mean(sky_emissivity .* proc.CONST.sigma .* (data_full.Tair(range,1) + 273.15 + proc.PARA.absolute_change_Tair).^4)];
%                 forcing.DATA.Sin = [forcing.DATA.Sin; mean(data_full.Sin(range,1) .*  (1+proc.PARA.relative_change_Sin))];
%                 
%                 sf = 0;
%                 rf = 0;
%                 for j=1:size(range,1)
%                     precip = data_full.snowfall(range(j),1) + data_full.rainfall(range(j),1);
%                     factor = max(0, min(1, (data_full.Tair(range(j),1) + proc.PARA.absolute_change_Tair - proc.PARA.all_snow_T) ./ max(1e-12, (proc.PARA.all_rain_T - proc.PARA.all_snow_T))));
%                     sf = sf + precip.*(1 - factor);
%                     rf = rf + precip.*factor;
%                 end
%                 forcing.DATA.snowfall = [forcing.DATA.snowfall; sf./size(range,1) .* proc.PARA.snow_fraction];
%                 forcing.DATA.rainfall = [forcing.DATA.rainfall; rf./size(range,1) .* proc.PARA.rain_fraction];
%                                 
%             end
%             
%             %overwrite target variables in TEMP in FORCING
%             forcing.TEMP = [];
%             forcing.TEMP.Lin = 0;
%             forcing.TEMP.Sin = 0;
%             forcing.TEMP.q = 0;
%             forcing.TEMP.p = 0;
%             forcing.TEMP.Tair = 0;
%             forcing.TEMP.wind = 0;
%             forcing.TEMP.rainfall = 0;
%             forcing.TEMP.snowfall = 0;
%             
%             forcing.DATA.full_res_data = data_full;
%         end
        
%                 %-------------param file generation-----
%         function proc = param_file_info(proc)
%             proc = provide_PARA(proc);
% 
%             proc.PARA.STATVAR = [];
%             proc.PARA.class_category = 'FORCING POST_PROCESSING';
%             proc.PARA.options = [];
%             
%             proc.PARA.eliminate_fraction = [];
%             proc.PARA.survive_fraction = [];
%                         
%             proc.PARA.default_value.window_size = {7};
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