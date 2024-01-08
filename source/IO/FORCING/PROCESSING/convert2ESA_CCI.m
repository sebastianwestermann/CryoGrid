%========================================================================
% CryoGrid FORCING processing class 
%
% compute ESA CCI forcing variables from "normal" forcing data, no ensemble
% generation
%
% Authors:
% S. Westermann, December 2022
%
%========================================================================

classdef convert2ESA_CCI < matlab.mixin.Copyable 
    
    properties
        
    end
    
    methods
        function proc = provide_PARA(proc)
            
            proc.PARA.averaging_period = [];  %in days

            
            proc.PARA.emissivity_snow = 0.99; % Snow emissivity (assumed known).
            
            proc.PARA.taus=0.0025; % Threshold snowfall for resetting to maximum [m w.e.].
            proc.PARA.taua=0.008; % Time constant for snow albedo change in non-melting conditions [/day].
            proc.PARA.tauf=0.24; % Time constant for snow albedo change in melting conditions [/day].
            
            proc.PARA.albsmax=0.85; % Maximum snow albedo.
            proc.PARA.albsmin=0.5; % Minimum snow albedo.
            
            proc.PARA.ar=0.2/1e2.* 3.34e8 ./ (24.*3600); % Restricted degree day factor (m*degC/day , value from Burbaker et al. 1996)
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
            proc.STATVAR.Lupwelling = proc.PARA.emissivity_snow.*proc.CONST.sigma.*proc.CONST.Tmfw.^4; % upwelling longwave radiation for melting snow, T=273.15K
            proc.STATVAR.albedo = proc.PARA.albsmax; %initialize 1st albedo values 
        end
        
        
        function forcing = process(proc, forcing, tile)
            
            data_full = forcing.DATA;
            forcing.DATA = [];
            forcing.DATA.snowfall = [];
            forcing.DATA.melt = [];
            forcing.DATA.surfT = [];
            forcing.DATA.timeForcing = [];
            
            for i = data_full.timeForcing(1,1):proc.PARA.averaging_period:data_full.timeForcing(end,1)-proc.PARA.averaging_period
                range = data_full.timeForcing>=i & data_full.timeForcing < min(data_full.timeForcing(end,1), i+ proc.PARA.averaging_period);
                forcing.DATA.timeForcing = [forcing.DATA.timeForcing; mean(data_full.timeForcing(range,1))];
                forcing.DATA.surfT = [forcing.DATA.surfT; mean(data_full.Tair(range,1))];
                forcing.DATA.snowfall = [forcing.DATA.snowfall; mean(data_full.snowfall(range,1))];
                
                melt_depth = 0;
                for j = 0:proc.PARA.averaging_period-1 %loop over individual days
                    range = data_full.timeForcing>=i+j & data_full.timeForcing < min(data_full.timeForcing(end,1), i+j+1);
                    
                    % Ablation term
                    LW_net = mean(proc.PARA.emissivity_snow.*data_full.Lin(range,1) - proc.STATVAR.Lupwelling); % Net  longwave
                    SW_net = mean((1-proc.STATVAR.albedo).*data_full.Sin(range,1)); % Net shortwave
                    SH_net = mean(proc.PARA.ar .* data_full.Tair(range,1)); % Warming through turbulent heat fluxes, parametrized using a restricted degree day approach.
                    
                    daily_melt_depth = (LW_net + SW_net + SH_net)  .* proc.CONST.day_sec ./ proc.CONST.L_f .* 1000;
                    melt_depth = melt_depth + daily_melt_depth; % Melt depth over the time step.

                    
                    % Update snow albedo for next step.
                    % Latest ECMWF "continuous reset" snow albedo scheme (Dutra et al. 2010)
                    new_snow = mean(data_full.snowfall(range,1)); %in mm/day

                    net_acc = new_snow - max(0,daily_melt_depth); % Net accumulation for one day time-step.
                    if net_acc>0
                        proc.STATVAR.albedo = proc.STATVAR.albedo + min(1,net_acc./(proc.PARA.taus .* 1000)) .* (proc.PARA.albsmax - proc.STATVAR.albedo);
                    elseif net_acc==0 %"Steady" case (linear decay)
                        proc.STATVAR.albedo = proc.STATVAR.albedo - proc.PARA.taua;
                    else
                        proc.STATVAR.albedo = (proc.STATVAR.albedo - proc.PARA.albsmin) .* exp(-proc.PARA.tauf) + proc.PARA.albsmin;
                        proc.STATVAR.albedo(proc.STATVAR.albedo < proc.PARA.albsmin) = proc.PARA.albsmin;
                    end

                end
                melt_depth(melt_depth <0) = 0;
                forcing.DATA.melt = [forcing.DATA.melt; melt_depth ./ proc.PARA.averaging_period];  %in mm/day
                
            end
            
            tile.PARA.geothermal = forcing.PARA.heatFlux_lb;
            
            %overwrite target variables in TEMP in FORCING
            forcing.TEMP = [];
            forcing.TEMP.snowfall=0;
            forcing.TEMP.melt = 0;
            forcing.TEMP.surfT = 0;
        end
        
        
%                 %-------------param file generation-----
%         function proc = param_file_info(post_proc)
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