%========================================================================
% CryoGrid FORCING processing class convert2ESA_CCI_ensemble2
%
% compute ESA CCI forcing variables from "normal" forcing data, with ensemble
% generation and calculation of sublimation

% Authors:
% S. Westermann, December 2022
%
%========================================================================

classdef convert2ESA_CCI_ensemble2 < FORCING_base
    
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
            proc.PARA.wind_speed_class = 5;
            proc.PARA.rain_fraction = 1;
            proc.PARA.relative_change_Sin = 1;     
            proc.PARA.relative_change_degree_day_factor = 1;
            
            proc.PARA.emissivity_snow = 0.99; % Snow emissivity (assumed known).
            
            proc.PARA.taus=0.0025; % Threshold snowfall for resetting to maximum [m w.e.].
            proc.PARA.taua=0.008; % Time constant for snow albedo change in non-melting conditions [/day].
            proc.PARA.tauf=0.24; % Time constant for snow albedo change in melting conditions [/day].
            
            proc.PARA.albsmax=0.85; % Maximum snow albedo.
            proc.PARA.albsmin=0.5; % Minimum snow albedo.
            
            proc.PARA.degree_day_factor=0.2/1e2 .* 3.34e8 ./ (24.*3600); % Restricted degree day factor (m*degC/day , value from Burbaker et al. 1996)  ./day_sec .* L_i
            proc.PARA.degree_day_factor=0.5/1e2 .* 3.34e8 ./ (24.*3600); 
        end
        
        
        function proc = provide_CONST(proc)
            proc.CONST.L_f = []; 
            proc.CONST.sigma = [];
            proc.CONST.day_sec = [];
            proc.CONST.Tmfw = [];
            proc.CONST.day_sec = [];
        end
        
        
        function proc = provide_STATVAR(proc)
            
        end
        
        
        function proc = finalize_init(proc, tile)

            proc.PARA.ensemble_size = tile.PARA.tile_size;
            proc.STATVAR.Lupwelling = proc.PARA.emissivity_snow.*proc.CONST.sigma.*proc.CONST.Tmfw.^4; % upwelling longwave radiation for melting snow, T=273.15K
            proc.STATVAR.albedo = repmat(proc.PARA.albsmax, 1, proc.PARA.ensemble_size); %initialize 1st albedo values 
            
            if size(proc.PARA.absolute_change_Tair,2)==1 && size(proc.PARA.absolute_change_Tair,1)==1
                proc.PARA.absolute_change_Tair = repmat(proc.PARA.absolute_change_Tair,1, proc.PARA.ensemble_size);
            elseif size(proc.PARA.absolute_change_Tair,2)==1 && size(proc.PARA.absolute_change_Tair,1) > 1
                proc.PARA.absolute_change_Tair = proc.PARA.absolute_change_Tair';
            end   
            if size(proc.PARA.snow_fraction,2)==1 && size(proc.PARA.snow_fraction,1)==1 
                proc.PARA.snow_fraction = repmat(proc.PARA.snow_fraction,1, proc.PARA.ensemble_size);
            elseif size(proc.PARA.snow_fraction,2)==1 && size(proc.PARA.snow_fraction,1) > 1 
                proc.PARA.snow_fraction = proc.PARA.snow_fraction';
            end
            if size(proc.PARA.rain_fraction,2)==1 && size(proc.PARA.rain_fraction,1)==1 
                proc.PARA.rain_fraction = repmat(proc.PARA.rain_fraction,1, proc.PARA.ensemble_size);
            elseif size(proc.PARA.rain_fraction,2)==1 && size(proc.PARA.rain_fraction,1) > 1 
                proc.PARA.rain_fraction = proc.PARA.rain_fraction';
            end
            if size(proc.PARA.relative_change_Sin,2)==1 && size(proc.PARA.relative_change_Sin,1)==1
                proc.PARA.relative_change_Sin = repmat(proc.PARA.relative_change_Sin,1, proc.PARA.ensemble_size);
            elseif size(proc.PARA.relative_change_Sin,2)==1 && size(proc.PARA.relative_change_Sin,1) > 1
                proc.PARA.relative_change_Sin = proc.PARA.relative_change_Sin';
            end
            if size(proc.PARA.wind_speed_class,2)==1 && size(proc.PARA.wind_speed_class,1)==1
                proc.PARA.wind_speed_class = repmat(proc.PARA.wind_speed_class,1, proc.PARA.ensemble_size);
            elseif size(proc.PARA.wind_speed_class,2)==1 && size(proc.PARA.wind_speed_class,1) > 1
                proc.PARA.wind_speed_class = proc.PARA.wind_speed_class';
            end
            if size(proc.PARA.relative_change_degree_day_factor,2)==1 && size(proc.PARA.relative_change_degree_day_factor,1)==1
                proc.PARA.relative_change_degree_day_factor = repmat(proc.PARA.relative_change_degree_day_factor,1, proc.PARA.ensemble_size);
            elseif size(proc.PARA.relative_change_degree_day_factor,2)==1 && size(proc.PARA.relative_change_degree_day_factor,1) > 1
                proc.PARA.relative_change_degree_day_factor = proc.PARA.relative_change_degree_day_factor';
            end
        end
        
        
        function forcing = process(proc, forcing, tile)
            
            data_full = forcing.DATA;
            forcing.DATA = [];
            forcing.DATA.snowfall = [];
            forcing.DATA.rainfall = [];
            forcing.DATA.sublimation = [];
            forcing.DATA.melt = [];
            forcing.DATA.surfT = [];
            forcing.DATA.timeForcing = [];
            forcing.DATA.albedo = [];
            albedo_reset = 0;

            
            for i = data_full.timeForcing(1,1):proc.PARA.averaging_period:data_full.timeForcing(end,1)-proc.PARA.averaging_period
                range = find(data_full.timeForcing>=i & data_full.timeForcing < min(data_full.timeForcing(end,1), i + proc.PARA.averaging_period));
                forcing.DATA.timeForcing = [forcing.DATA.timeForcing; mean(data_full.timeForcing(range,1))];
                forcing.DATA.surfT = [forcing.DATA.surfT; mean(data_full.Tair(range,1)) + proc.PARA.absolute_change_Tair];
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
                
                %sublimation
                rho_a = rho_air(proc, data_full.p(range,1), min(0, data_full.Tair(range,1) + proc.PARA.absolute_change_Tair) + 273.15);
                rho_w = 1000;
                kappa = 0.4;
                q_surf = 0.622.*satPresIce(proc, min(0, data_full.Tair(range,1)+proc.PARA.absolute_change_Tair) + 273.15) ./ data_full.p(range,1);
                sublimation = rho_a ./ rho_w.*kappa.^2.*proc.PARA.wind_speed_class ./(log(2./1e-3)).^2 .* max(0, (-data_full.q(range,1) + q_surf));
                forcing.DATA.sublimation = [forcing.DATA.sublimation; mean(sublimation) .*1000 .* proc.CONST.day_sec];
                
                melt_depth = 0;
              %  sublimation2 = 0;
                for j = 0:proc.PARA.averaging_period-1 %loop over individual days
                    range = find(data_full.timeForcing>=i+j & data_full.timeForcing < min(data_full.timeForcing(end,1), i+j+1));
                    
                    % Ablation term
                    Lin = 0;
                    Sin = 0;
                    sf = 0;
                    for k=1:size(range,1)
                        sky_emissivity = data_full.Lin(range(k),1) ./ (data_full.Tair(range(k),1)+273.15).^4 ./ proc.CONST.sigma;
                        Lin = Lin + sky_emissivity .* proc.CONST.sigma .* (data_full.Tair(range(k),1) + 273.15 + proc.PARA.absolute_change_Tair).^4;
                        Sin = Sin + data_full.Sin(range(k),1) .*  proc.PARA.relative_change_Sin;
                        precip = data_full.snowfall(range(k),1) + data_full.rainfall(range(k),1);
                        factor = max(0, min(1, (data_full.Tair(range(k),1) + proc.PARA.absolute_change_Tair - proc.PARA.all_snow_T) ./ max(1e-12, (proc.PARA.all_rain_T - proc.PARA.all_snow_T))));
                        sf = sf + precip.*(1 - factor);
                        
%                         %sublimation
%                         surf_T = get_surf_T(proc, forcing, Sin, Lin, data_full.Tair(range(k),1) + proc.PARA.absolute_change_Tair, proc.STATVAR.albedo);
%                         rho_a = rho_air(proc, data_full.p(range(k),1), min(0, data_full.Tair(range(k),1) + proc.PARA.absolute_change_Tair) + 273.15);
%                         q_surf = 0.622.*satPresIce(proc, min(0, surf_T) + 273.15) ./ data_full.p(range(k),1);
%                         sublimation2 = sublimation2 +  max(0, rho_a ./ rho_w.*kappa.^2.*proc.PARA.wind_speed_class ./(log(2./1e-3)).^2 .* max(0, (-data_full.q(range(k),1) + q_surf)));
                    end
                    
                    LW_net = proc.PARA.emissivity_snow .* Lin ./ size(range,1) - proc.STATVAR.Lupwelling; % Net  longwave in W/m2
                    SW_net = (1-proc.STATVAR.albedo) .* Sin ./ size(range,1); % Net shortwave
                    %SH_net = proc.PARA.relative_change_degree_day_factor .* proc.PARA.degree_day_factor .* mean(data_full.Tair(range,1)); % Warming through turbulent heat fluxes, parametrized using a restricted degree day approach.
                    SH_net = 8e-3.*1005 .* (mean(data_full.Tair(range,1)) + proc.PARA.absolute_change_Tair); % Warming through turbulent heat fluxes, parametrized using a restricted degree day approach.
                    
                    daily_melt_depth = (LW_net + SW_net + SH_net) .* proc.CONST.day_sec ./ proc.CONST.L_f .* 1000; %in mm/day
                    melt_depth = melt_depth + daily_melt_depth; % Melt depth over the time step.

                    
                    % Update snow albedo for next step.
                    % Latest ECMWF "continuous reset" snow albedo scheme (Dutra et al. 2010)
                    new_snow = sf./size(range,1) .* proc.PARA.snow_fraction; % mean(data_full.snowfall(range,1)); %in mm/day

                    net_acc = new_snow - max(0,daily_melt_depth); % Net accumulation for one day time-step.
                    constr = net_acc>0;
                    proc.STATVAR.albedo(1, constr) = proc.STATVAR.albedo(1, constr) + min(1,net_acc(1, constr)./(proc.PARA.taus .* 1000)) .* (proc.PARA.albsmax - proc.STATVAR.albedo(1, constr));
                    constr = net_acc==0; %"Steady" case (linear decay)
                    proc.STATVAR.albedo(1, constr) = proc.STATVAR.albedo(1, constr) - proc.PARA.taua;
                    constr = net_acc<0;
                    proc.STATVAR.albedo(1, constr) = (proc.STATVAR.albedo(1, constr) - proc.PARA.albsmin) .* exp(-proc.PARA.tauf) + proc.PARA.albsmin;
                    proc.STATVAR.albedo(proc.STATVAR.albedo < proc.PARA.albsmin) = proc.PARA.albsmin;
                    if albedo_reset == 1 && month(i)>=9
                        albedo_reset=0;
                        proc.STATVAR.albedo = proc.STATVAR.albedo .*0 + proc.PARA.albsmax;
                    elseif albedo_reset == 0 && month(i)>=12 
                        albedo_reset=1;
                    end

                end
               % melt_depth(melt_depth <0) = 0;
                forcing.DATA.melt = [forcing.DATA.melt; melt_depth ./ proc.PARA.averaging_period];  %in mm/day
                forcing.DATA.albedo = [forcing.DATA.albedo; proc.STATVAR.albedo];
           %     forcing.DATA.sublimation2 = [forcing.DATA.sublimation2; sublimation2 ./size(range,1) ./ proc.PARA.averaging_period .* 1000 .* proc.CONST.day_sec];
            end
            
            if size(forcing.PARA.heatFlux_lb,2) == 1
                tile.PARA.geothermal = repmat(forcing.PARA.heatFlux_lb, 1, proc.PARA.ensemble_size);
            end
            
            %overwrite target variables in TEMP in FORCING
            forcing.TEMP = [];
            forcing.TEMP.snowfall=0;
            forcing.TEMP.rainfall=0;
            forcing.TEMP.melt = 0;
            forcing.TEMP.sublimation = [];
            forcing.TEMP.surfT = 0;
        end
        
        function p = satPresIce(proc, T) %saturation pressure ice, Magnus formula
            p = 6.112.* 100.* exp(22.46.*(T-273.15)./(272.61-273.15+T));
        end
        
        function rho = rho_air(proc, p, T) %air density [kg m^(-3)]
            rho = p./(287.058.*T); 
        end
        
%         function surf_T = get_surf_T(proc, forcing, Sin, Lin, Tair,
%         albedo) %calculates analytic surface T taking only radiation and
%         sensible heat flux into account -> sublimation calulated hereof
%         looks too big, so it seems like the normal way of calculating
%         sublimation is more successful%             
%             
%             a = -proc.PARA.emissivity_snow .* forcing.CONST.sigma;
%             d = - proc.PARA.relative_change_degree_day_factor .* proc.PARA.degree_day_factor;
%             e = proc.PARA.emissivity_snow .* Lin + (1-albedo) .* Sin + proc.PARA.relative_change_degree_day_factor .* proc.PARA.degree_day_factor .* (Tair + 273.15);
%             
%             q = d./a;
%             delta_0 = 12.*a.*e;
%             delta_1 = 27 .* a.* d.^2;
%             delta = 256 .* a.^3 .* e.^3 - 27 .* a.^2 .* d.^4;
%             Q = (0.5.*(delta_1 + sqrt(delta_1.^2-4.*delta_0.^3))).^(1./3);
%             S = 0.5 .* sqrt(1./(3.*a).*(Q + delta_0./Q));
%             
%             surf_T = -S + 0.5.*sqrt(-4.* S.^2 + q ./ S)-273.15;
%             
%         end
        
%                 %-------------param file generation-----
%         function proc = param_file_info(proc)
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