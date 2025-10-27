%========================================================================
% CryoGrid FORCING processing class interpolate_sl_NORA3
% uses interpolation of Sin relying on S_TOA
%
% Authors:
% S. Westermann, December 2022
%
%========================================================================

classdef interpolate_sl_NORA3 < process_BASE
    
    
    methods
        function proc = provide_PARA(proc)
            
        end
        
        
        function proc = provide_CONST(proc)
            proc.CONST.Tmfw = [];
            proc.CONST.sigma = [];
        end
        
        
        function proc = provide_STATVAR(proc)
            
        end
        
        
        function proc = finalize_init(proc, tile)
            
        end
        
        
        function forcing = process(proc, forcing, tile)
            
            disp('interpolating surface level data')
            nora = forcing.TEMP.nora;
            
            [dist,~] = distance(forcing.SPATIAL.STATVAR.latitude,forcing.SPATIAL.STATVAR.longitude,nora.latitude,nora.longitude);
            
            [dist, ind] = sort(dist(:));
            
            dist=dist(1:4);
            ind = ind(1:4);
            [row, col] = ind2sub(size(nora.latitude), ind);
            weights = 1./dist ./ sum(1./dist);
            
            time_resolution = nora.t(2) - nora.t(1);


            wind = sqrt((double(nora.u10) .* nora.wind_sf).^2 + (double(nora.v10) .* nora.wind_sf).^2);
            wind = squeeze(wind(row(1), col(1), :) .* weights(1) + wind(row(2), col(2),:) .* weights(2) + wind(row(3), col(3),:) .* weights(3) + wind(row(4), col(4), :) .* weights(4));
            

            Sin = double(nora.SW) .* nora.rad_sf;
            
            S0=1370; % Solar constant (total TOA solar irradiance) [Wm^-2] used in ECMWF's IFS
            S_TOA  = []; 
            Sin2 = [];
            %compute actual S_TOA for NORA ccordinates
            for i=1:4
                [~, solar_zen] = solargeom(proc, [nora.t(1):1/240:nora.t(end) + time_resolution]', nora.latitude(row(i), col(i)), nora.longitude(row(i), col(i)));
                solar_zen = solar_zen(1:size(nora.t,1).*24.*time_resolution.*10,1);
                mu0=max(cos(solar_zen),0); % Trunacte negative values.
                mu0 = mean(reshape(mu0, 10*24.*time_resolution, size(solar_zen,1)/(10*24.*time_resolution)))';
                S_TOA = [S_TOA (S0.*mu0)];
                Sin2 = [Sin2 squeeze(Sin(row(i), col(i),:))];
            end
            S_TOA = S_TOA(1:size(Sin2,1),:);
            transmissivity = Sin2./S_TOA;
            transmissivity(transmissivity>=1) = NaN;           
            transmissivity = [transmissivity(1,:); transmissivity];    
            nan_transmissivity = double(~isnan(transmissivity));

            transmissivity(isnan(transmissivity)) = 0;
            transmissivity = (transmissivity(1:end-1,:) + transmissivity(2:end,:)) ./  (nan_transmissivity(1:end-1,:) + nan_transmissivity(2:end,:));  %average sky transmissivity at the timestamp


            transmissivity(isnan(transmissivity)) = nanmean((transmissivity(:)));
            weights_Sin = Sin2 ./ repmat(sum(Sin2,2), 1, 4);
            weights_Sin(isnan(weights_Sin)) = 0.25;
            weights_Sin = weights_Sin .* repmat(weights', size(weights_Sin,1),1);
            weights_Sin = weights_Sin ./ repmat(sum(weights_Sin, 2), 1, 4);

            transmissivity_point = sum(transmissivity .* weights_Sin,2);


            %compute actual S_TOA for target point
            [~, solar_zen]=solargeom(proc, nora.t, forcing.SPATIAL.STATVAR.latitude, forcing.SPATIAL.STATVAR.longitude);
            
            mu0=max(cos(solar_zen),0); % Trunacte negative values.
            mu0 = mu0(1:size(transmissivity_point,1),:);
            Sin_final = S0.*mu0 .* transmissivity_point;



            %calculate average TOA radiation in integration time window for the four cells
            %determine transmissivity through atmosphere, NaN if Sin < threshold
            %interpolate transmissivityin space, taking only the non-NaN-values
            %determine TOA at the target location  and multiply by
            %transmissivity

            
            Tair = double(nora.T2) .* nora.T_sf;

            Lin = double(nora.LW) .* nora.rad_sf;
            Tair2 = cat(3, Tair, Tair(:,:,end));
            Tair2 = (Tair2(:,:, 1:end-1) + Tair2(:,:, 2:end)) ./ 2; %average air T during the interval over which LW in is averaged

            sky_emissivity = Lin ./ proc.CONST.sigma ./ (Tair2 + proc.CONST.Tmfw).^4;
            sky_emissivity = cat(3, sky_emissivity(:,:,1), sky_emissivity);
            sky_emissivity = (sky_emissivity(:,:, 1:end-1) + sky_emissivity(:,:, 2:end)) ./ 2;  %average sky emissivity at the timestamp
            Lin = sky_emissivity .* proc.CONST.sigma .* (Tair + proc.CONST.Tmfw).^4;
            Lin = squeeze(Lin(row(1), col(1), :) .* weights(1) + Lin(row(2), col(2),:) .* weights(2) + Lin(row(3), col(3),:) .* weights(3) + Lin(row(4), col(4), :) .* weights(4));
            
            RH = double(nora.RH) .* nora.T_sf;
            RH = squeeze(RH(row(1), col(1), :) .* weights(1) + RH(row(2), col(2),:) .* weights(2) + RH(row(3), col(3),:) .* weights(3) + RH(row(4), col(4), :) .* weights(4));
            
            precip = double(nora.P) .* nora.P_sf;
            precip = squeeze(precip(row(1), col(1), :) .* weights(1) + precip(row(2), col(2),:) .* weights(2) + precip(row(3), col(3),:) .* weights(3) + precip(row(4), col(4), :) .* weights(4));
            precip = [precip(1); precip]; %precip accumulated in time interval after timestamp
            precip = (precip(1:end-1)+precip(2:end)) ./ 2;

            
            pressure = double(nora.ps) .* nora.ps_sf;
            pressure = squeeze(pressure(row(1), col(1), :) .* weights(1) + pressure(row(2), col(2),:) .* weights(2) + pressure(row(3), col(3),:) .* weights(3) + pressure(row(4), col(4), :) .* weights(4));

            
            Tair = squeeze(Tair(row(1), col(1), :) .* weights(1) + Tair(row(2), col(2),:) .* weights(2) + Tair(row(3), col(3),:) .* weights(3) + Tair(row(4), col(4), :) .* weights(4));

            
            forcing.DATA.Tair = Tair;
            forcing.DATA.RH = RH;
            forcing.DATA.wind = wind;
            forcing.DATA.Sin =  Sin_final;
            forcing.DATA.Lin = Lin;
            forcing.DATA.p = pressure;
            forcing.DATA.precip = precip;  %mm/hour to mm/day
            forcing.DATA.timeForcing = nora.t;
            
            forcing.TEMP.nora = [];
            
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

