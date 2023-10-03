%========================================================================
% CryoGrid FORCING post-processing class 

% Authors:
% S. Westermann, December 2022
%
%========================================================================

classdef convert2ESA_CCI_ensemble2 < FORCING_base 
    
    properties
        
    end
    
    methods
        function post_proc = provide_PARA(post_proc)
            
            post_proc.PARA.averaging_period = [];  %in days
            post_proc.PARA.all_snow_T = [];
            post_proc.PARA.all_rain_T = [];
            
            %these can be written by ensemble class
            post_proc.PARA.ensemble_size = 1;
            post_proc.PARA.absolute_change_Tair = 0;
            post_proc.PARA.snow_fraction = 1;
            post_proc.PARA.wind_speed_class = 5;
            post_proc.PARA.rain_fraction = 1;
            post_proc.PARA.relative_change_Sin = 1;     
            post_proc.PARA.relative_change_degree_day_factor = 1;
            
            post_proc.PARA.emissivity_snow = 0.99; % Snow emissivity (assumed known).
            
            post_proc.PARA.taus=0.0025; % Threshold snowfall for resetting to maximum [m w.e.].
            post_proc.PARA.taua=0.008; % Time constant for snow albedo change in non-melting conditions [/day].
            post_proc.PARA.tauf=0.24; % Time constant for snow albedo change in melting conditions [/day].
            
            post_proc.PARA.albsmax=0.85; % Maximum snow albedo.
            post_proc.PARA.albsmin=0.5; % Minimum snow albedo.
            
            post_proc.PARA.degree_day_factor=0.2/1e2 .* 3.34e8 ./ (24.*3600); % Restricted degree day factor (m*degC/day , value from Burbaker et al. 1996)  ./day_sec .* L_i
        end
        
        
        function post_proc = provide_CONST(post_proc)
            post_proc.CONST.L_f = []; 
            post_proc.CONST.sigma = [];
            post_proc.CONST.day_sec = [];
            post_proc.CONST.Tmfw = [];
            post_proc.CONST.day_sec = [];
        end
        
        
        function post_proc = provide_STATVAR(post_proc)
            
        end
        
        
        function post_proc = finalize_init(post_proc, tile)
            post_proc.STATVAR.Lupwelling = post_proc.PARA.emissivity_snow.*post_proc.CONST.sigma.*post_proc.CONST.Tmfw.^4; % upwelling longwave radiation for melting snow, T=273.15K
            post_proc.STATVAR.albedo = repmat(post_proc.PARA.albsmax, 1, post_proc.PARA.ensemble_size); %initialize 1st albedo values 
            
            if size(post_proc.PARA.absolute_change_Tair,2)==1 
                post_proc.PARA.absolute_change_Tair = repmat(post_proc.PARA.absolute_change_Tair,1, post_proc.PARA.ensemble_size);
            end   
            if size(post_proc.PARA.wind_speed_class,2)==1
                post_proc.PARA.wind_speed_class = repmat(post_proc.PARA.wind_speed_class,1, post_proc.PARA.ensemble_size);
            end
            if size(post_proc.PARA.snow_fraction,2)==1 
                post_proc.PARA.snow_fraction = repmat(post_proc.PARA.snow_fraction,1, post_proc.PARA.ensemble_size);
            end
            if size(post_proc.PARA.rain_fraction,2)==1
                post_proc.PARA.rain_fraction = repmat(post_proc.PARA.rain_fraction,1, post_proc.PARA.ensemble_size);
            end
            if size(post_proc.PARA.relative_change_Sin,2)==1
                post_proc.PARA.relative_change_Sin = repmat(post_proc.PARA.relative_change_Sin,1, post_proc.PARA.ensemble_size);
            end
            if size(post_proc.PARA.relative_change_degree_day_factor,2)==1
                post_proc.PARA.relative_change_degree_day_factor = repmat(post_proc.PARA.relative_change_degree_day_factor,1, post_proc.PARA.ensemble_size);
            end
        end
        
        
        function forcing = post_process(post_proc, forcing, tile)
            
            data_full = forcing.DATA;
            forcing.DATA = [];
            forcing.DATA.snowfall = [];
            forcing.DATA.rainfall = [];
            forcing.DATA.sublimation = [];
            forcing.DATA.sublimation2 = [];            
            forcing.DATA.melt = [];
            forcing.DATA.surfT = [];
            forcing.DATA.timeForcing = [];
            
            for i = data_full.timeForcing(1,1):post_proc.PARA.averaging_period:data_full.timeForcing(end,1)-post_proc.PARA.averaging_period
                range = find(data_full.timeForcing>=i & data_full.timeForcing < min(data_full.timeForcing(end,1), i + post_proc.PARA.averaging_period));
                forcing.DATA.timeForcing = [forcing.DATA.timeForcing; mean(data_full.timeForcing(range,1))];
                forcing.DATA.surfT = [forcing.DATA.surfT; mean(data_full.Tair(range,1)) + post_proc.PARA.absolute_change_Tair];
                sf = 0;
                rf = 0;
                for j=1:size(range,1)
                    precip = data_full.snowfall(range(j),1) + data_full.rainfall(range(j),1);
                    factor = max(0, min(1, (data_full.Tair(range(j),1) + post_proc.PARA.absolute_change_Tair - post_proc.PARA.all_snow_T) ./ max(1e-12, (post_proc.PARA.all_rain_T - post_proc.PARA.all_snow_T))));
                    sf = sf + precip.*(1 - factor);
                    rf = rf + precip.*factor;
                end
                forcing.DATA.snowfall = [forcing.DATA.snowfall; sf./size(range,1) .* post_proc.PARA.snow_fraction];
                forcing.DATA.rainfall = [forcing.DATA.rainfall; rf./size(range,1) .* post_proc.PARA.rain_fraction];
                
                %sublimation
                rho_a = rho_air(post_proc, data_full.p(range,1), min(0, data_full.Tair(range,1) + post_proc.PARA.absolute_change_Tair) + 273.15);
                rho_w = 1000;
                kappa = 0.4;
                q_surf = 0.622.*satPresIce(post_proc, min(0, data_full.Tair(range,1)+post_proc.PARA.absolute_change_Tair) + 273.15) ./ data_full.p(range,1);
                sublimation = rho_a ./ rho_w.*kappa.^2.*post_proc.PARA.wind_speed_class ./(log(2./1e-3)).^2 .* max(0, (-data_full.q(range,1) + q_surf));
                forcing.DATA.sublimation = [forcing.DATA.sublimation; mean(sublimation) .*1000 .* post_proc.CONST.day_sec];
                
                melt_depth = 0;
              %  sublimation2 = 0;
                for j = 0:post_proc.PARA.averaging_period-1 %loop over individual days
                    range = find(data_full.timeForcing>=i+j & data_full.timeForcing < min(data_full.timeForcing(end,1), i+j+1));
                    
                    % Ablation term
                    Lin = 0;
                    Sin = 0;
                    sf = 0;
                    for k=1:size(range,1)
                        sky_emissivity = data_full.Lin(range(k),1) ./ (data_full.Tair(range(k),1)+273.15).^4 ./ post_proc.CONST.sigma;
                        Lin = Lin + sky_emissivity .* post_proc.CONST.sigma .* (data_full.Tair(range(k),1) + 273.15 + post_proc.PARA.absolute_change_Tair).^4;
                        Sin = Sin + data_full.Sin(range(k),1) .*  post_proc.PARA.relative_change_Sin;
                        precip = data_full.snowfall(range(k),1) + data_full.rainfall(range(k),1);
                        factor = max(0, min(1, (data_full.Tair(range(k),1) - post_proc.PARA.all_snow_T) ./ max(1e-12, (post_proc.PARA.all_rain_T - post_proc.PARA.all_snow_T))));
                        sf = sf + precip.*(1 - factor);
                        
%                         %sublimation
%                         surf_T = get_surf_T(post_proc, forcing, Sin, Lin, data_full.Tair(range(k),1) + post_proc.PARA.absolute_change_Tair, post_proc.STATVAR.albedo);
%                         rho_a = rho_air(post_proc, data_full.p(range(k),1), min(0, data_full.Tair(range(k),1) + post_proc.PARA.absolute_change_Tair) + 273.15);
%                         q_surf = 0.622.*satPresIce(post_proc, min(0, surf_T) + 273.15) ./ data_full.p(range(k),1);
%                         sublimation2 = sublimation2 +  max(0, rho_a ./ rho_w.*kappa.^2.*post_proc.PARA.wind_speed_class ./(log(2./1e-3)).^2 .* max(0, (-data_full.q(range(k),1) + q_surf)));
                    end
                    
                    LW_net = post_proc.PARA.emissivity_snow .* Lin ./ size(range,1) - post_proc.STATVAR.Lupwelling; % Net  longwave
                    SW_net = (1-post_proc.STATVAR.albedo) .* Sin ./ size(range,1); % Net shortwave
                    SH_net = post_proc.PARA.relative_change_degree_day_factor .* post_proc.PARA.degree_day_factor .* mean(data_full.Tair(range,1)); % Warming through turbulent heat fluxes, parametrized using a restricted degree day approach.
                    
                    daily_melt_depth = (LW_net + SW_net + SH_net) .* post_proc.CONST.day_sec ./ post_proc.CONST.L_f .* 1000;
                    melt_depth = melt_depth + daily_melt_depth; % Melt depth over the time step.

                    
                    % Update snow albedo for next step.
                    % Latest ECMWF "continuous reset" snow albedo scheme (Dutra et al. 2010)
                    new_snow = sf./size(range,1) .* post_proc.PARA.snow_fraction; % mean(data_full.snowfall(range,1)); %in mm/day

                    net_acc = new_snow - max(0,daily_melt_depth); % Net accumulation for one day time-step.
                    constr = net_acc>0;
                    post_proc.STATVAR.albedo(1, constr) = post_proc.STATVAR.albedo(1, constr) + min(1,net_acc(1, constr)./(post_proc.PARA.taus .* 1000)) .* (post_proc.PARA.albsmax - post_proc.STATVAR.albedo(1, constr));
                    constr = net_acc==0; %"Steady" case (linear decay)
                    post_proc.STATVAR.albedo(1, constr) = post_proc.STATVAR.albedo(1, constr) - post_proc.PARA.taua;
                    constr = net_acc<0;
                    post_proc.STATVAR.albedo(1, constr) = (post_proc.STATVAR.albedo(1, constr) - post_proc.PARA.albsmin) .* exp(-post_proc.PARA.tauf) + post_proc.PARA.albsmin;
                    post_proc.STATVAR.albedo(post_proc.STATVAR.albedo < post_proc.PARA.albsmin) = post_proc.PARA.albsmin;

                end
               % melt_depth(melt_depth <0) = 0;
                forcing.DATA.melt = [forcing.DATA.melt; melt_depth ./ post_proc.PARA.averaging_period];  %in mm/day
           %     forcing.DATA.sublimation2 = [forcing.DATA.sublimation2; sublimation2 ./size(range,1) ./ post_proc.PARA.averaging_period .* 1000 .* post_proc.CONST.day_sec];
            end
            
            if size(forcing.PARA.heatFlux_lb,2) == 1
                tile.PARA.geothermal = repmat(forcing.PARA.heatFlux_lb, 1, post_proc.PARA.ensemble_size);
            end
            
            %overwrite target variables in TEMP in FORCING
            forcing.TEMP = [];
            forcing.TEMP.snowfall=0;
            forcing.TEMP.rainfall=0;
            forcing.TEMP.melt = 0;
            forcing.TEMP.sublimation = [];
            forcing.TEMP.surfT = 0;
        end
        
        function p = satPresIce(post_proc, T) %saturation pressure ice, Magnus formula
            p = 6.112.* 100.* exp(22.46.*(T-273.15)./(272.61-273.15+T));
        end
        
        function rho = rho_air(post_proc, p, T) %air density [kg m^(-3)]
            rho = p./(287.058.*T); 
        end
        
%         function surf_T = get_surf_T(post_proc, forcing, Sin, Lin, Tair,
%         albedo) %calculates analytic surface T taking only radiation and
%         sensible heat flux into account -> sublimation calulated hereof
%         looks too big, so it seems like the normal way of calculating
%         sublimation is more successful%             
%             
%             a = -post_proc.PARA.emissivity_snow .* forcing.CONST.sigma;
%             d = - post_proc.PARA.relative_change_degree_day_factor .* post_proc.PARA.degree_day_factor;
%             e = post_proc.PARA.emissivity_snow .* Lin + (1-albedo) .* Sin + post_proc.PARA.relative_change_degree_day_factor .* post_proc.PARA.degree_day_factor .* (Tair + 273.15);
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