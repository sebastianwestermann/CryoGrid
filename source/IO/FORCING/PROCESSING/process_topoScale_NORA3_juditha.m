%========================================================================
% CryoGrid FORCING processing class
%%uses "normal" interpolation of Sin, not the superior algroithm relying on
%S_TOA used in process_topoScale2
%
% Authors:
% S. Westermann, December 2022
%
%========================================================================

classdef process_topoScale_NORA3_juditha < process_BASE
    

    
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
            
                        
            disp('applying downscaling with TopoScale')
            nora = forcing.TEMP.nora;
            
            [dist,~] = distance(forcing.SPATIAL.STATVAR.latitude,forcing.SPATIAL.STATVAR.longitude,nora.latitude,nora.longitude);
            [dist, ind] = sort(dist(:));
            dist=dist(1:2);
            ind = ind(1:2);
            [row, col] = ind2sub(size(nora.latitude), ind);
            weights = 1./dist ./ sum(1./dist);
            
            row=[row; 1; 1];
            col = [col; 1;1];
            weights = [weights; 0;0];


            time_resolution = nora.t(2) - nora.t(1);

            % shortwave_in surface
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
            sunset=mu0<cosd(89); % Sunset switch.
            S_TOA_point = S0.*mu0;
            Sin =  transmissivity_point .* S_TOA_point;


            %longwave_in surface
            Tair = double(nora.T2) .* nora.T_sf;

            Lin = double(nora.LW) .* nora.rad_sf;
            Tair2 = cat(3, Tair, Tair(:,:,end));
            Tair2 = (Tair2(:,:, 1:end-1) + Tair2(:,:, 2:end)) ./ 2; %average air T during the interval over which LW in is averaged

            sky_emissivity = Lin ./ proc.CONST.sigma ./ (Tair2 + proc.CONST.Tmfw).^4;
            sky_emissivity = cat(3, sky_emissivity(:,:,1), sky_emissivity);
            sky_emissivity = (sky_emissivity(:,:, 1:end-1) + sky_emissivity(:,:, 2:end)) ./ 2;  %average sky emissivity at the timestamp
            Lin = sky_emissivity .* proc.CONST.sigma .* (Tair + proc.CONST.Tmfw).^4;
            Lin = squeeze(Lin(row(1), col(1), :) .* weights(1) + Lin(row(2), col(2),:) .* weights(2) + Lin(row(3), col(3),:) .* weights(3) + Lin(row(4), col(4), :) .* weights(4));
            
            

            
            Tair_sl = squeeze(Tair(row(1), col(1), :) .* weights(1) + Tair(row(2), col(2),:) .* weights(2) + Tair(row(3), col(3),:) .* weights(3) + Tair(row(4), col(4), :) .* weights(4));
            Tair_pl = double(nora.T) .* nora.T_sf;
            Tair_pl = squeeze(Tair_pl(row(1), col(1), :,:) .* weights(1) + Tair_pl(row(2), col(2),:,:) .* weights(2) + Tair_pl(row(3), col(3), :, :) .* weights(3) + Tair_pl(row(4), col(4), :, :) .* weights(4));
            Tair_pl = Tair_pl';

            wind_sl = sqrt((double(nora.u10) .* nora.wind_sf).^2 + (double(nora.v10) .* nora.wind_sf).^2);
            wind_sl = squeeze(wind_sl(row(1), col(1), :) .* weights(1) + wind_sl(row(2), col(2),:) .* weights(2) + wind_sl(row(3), col(3),:) .* weights(3) + wind_sl(row(4), col(4), :) .* weights(4));
            wind_pl = sqrt((double(nora.u) .* nora.wind_sf).^2 +(double(nora.v) .* nora.wind_sf).^2);
            wind_pl = squeeze(wind_pl(row(1), col(1), :,:) .* weights(1) + wind_pl(row(2), col(2),:,:) .* weights(2) + wind_pl(row(3), col(3), :, :) .* weights(3) + wind_pl(row(4), col(4), :, :) .* weights(4));
            wind_pl = wind_pl';

            RH_sl = double(nora.RH) .* nora.T_sf;
            RH_sl = squeeze(RH_sl(row(1), col(1), :) .* weights(1) + RH_sl(row(2), col(2),:) .* weights(2) + RH_sl(row(3), col(3),:) .* weights(3) + RH_sl(row(4), col(4), :) .* weights(4));
            RH_pl = double(nora.RH_pl) .* nora.T_sf;
            RH_pl = squeeze(RH_pl(row(1), col(1), :,:) .* weights(1) + RH_pl(row(2), col(2),:,:) .* weights(2) + RH_pl(row(3), col(3), :, :) .* weights(3) + RH_pl(row(4), col(4), :, :) .* weights(4));
            RH_pl = RH_pl';

            altitude_sl = double(nora.Zs(row(1), col(1)) .* weights(1) + nora.Zs(row(2), col(2)) .* weights(2) + nora.Zs(row(3), col(3)) .* weights(3) + nora.Zs(row(4), col(4)) .* weights(4));            
            altitude_pl = double(nora.Z);
            altitude_pl = squeeze(altitude_pl(row(1), col(1), :,:) .* weights(1) + altitude_pl(row(2), col(2),:,:) .* weights(2) + altitude_pl(row(3), col(3), :, :) .* weights(3) + altitude_pl(row(4), col(4), :, :) .* weights(4));            
            altitude_pl = altitude_pl';

            precip = double(nora.P) .* nora.P_sf;
            precip = squeeze(precip(row(1), col(1), :) .* weights(1) + precip(row(2), col(2),:) .* weights(2) + precip(row(3), col(3),:) .* weights(3) + precip(row(4), col(4), :) .* weights(4));
            precip = [precip(1); precip]; %precip accumulated in time interval after timestamp
            precip = (precip(1:end-1)+precip(2:end)) ./ 2;

            pressure_sl = double(nora.ps) .* nora.ps_sf;
            pressure_sl = squeeze(pressure_sl(row(1), col(1), :) .* weights(1) + pressure_sl(row(2), col(2),:) .* weights(2) + pressure_sl(row(3), col(3),:) .* weights(3) + pressure_sl(row(4), col(4), :) .* weights(4));
            pressure_pl = double(nora.p) .* 100;
            



            layer_below = altitude_pl.*0;
            layer_above = altitude_pl.*0;
                        
            for i=2:size(altitude_pl,2)
                layer_below(:,i) = double(altitude_pl(:,i) < forcing.SPATIAL.STATVAR.altitude & altitude_pl(:,i-1) >= forcing.SPATIAL.STATVAR.altitude & altitude_pl(:,i) > altitude_sl);
                layer_above(:,i-1) = double(altitude_pl(:,i) < forcing.SPATIAL.STATVAR.altitude & altitude_pl(:,i-1) >= forcing.SPATIAL.STATVAR.altitude & altitude_pl(:,i-1) > altitude_sl);
            end
            
            distance_Z_above = abs(sum(altitude_pl .* layer_above ,2) - forcing.SPATIAL.STATVAR.altitude) .* double(sum(layer_above,2) > 0);
            distance_Z_below = abs(sum(altitude_pl .* layer_below ,2) - forcing.SPATIAL.STATVAR.altitude) .* double(sum(layer_below,2) > 0);
            
            weights_Z_above = 1-distance_Z_above ./ max(1e-10, distance_Z_above + distance_Z_below);
            weights_Z_below = 1-distance_Z_below ./ max(1e-10, distance_Z_above + distance_Z_below);
            weights_Z_above = repmat(weights_Z_above, 1, size(layer_above,2),1) .* double(layer_above);
            weights_Z_below = repmat(weights_Z_below, 1, size(layer_below,2),1) .* double(layer_below);

            Tair = sum(Tair_pl .*  (weights_Z_above + weights_Z_below), 2);
            wind = sum(wind_pl .*  (weights_Z_above + weights_Z_below), 2);
            RH = sum(RH_pl .*   (weights_Z_above + weights_Z_below), 2);
            
            use_sl = sum(weights_Z_above + weights_Z_below,2) <1-1e-9  | forcing.SPATIAL.STATVAR.altitude < altitude_sl;
            merge_w_sl = sum(weights_Z_below,2) == 0;
            factor = min(1, max(0, distance_Z_above./100));
            
            Tair = double(~merge_w_sl) .* Tair + double(merge_w_sl) .* factor .* Tair + double(merge_w_sl) .* (1-factor) .* Tair_sl;
            Tair = double(~use_sl) .* Tair + double(use_sl) .* Tair_sl;
            wind = double(~merge_w_sl) .* wind + double(merge_w_sl) .* factor .* wind + double(merge_w_sl) .* (1-factor) .* wind_sl;
            wind = double(~use_sl) .* wind + double(use_sl) .* wind_sl;
            RH = double(~merge_w_sl) .* RH + double(merge_w_sl) .* factor .* RH + double(merge_w_sl) .* (1-factor) .* RH_sl;
            RH = double(~use_sl) .* RH + double(use_sl) .* RH_sl;

            g = 9.81;
            M = 0.0289644;
            R_gc = 8.3144598;
            pressure = sum( repmat(pressure_pl,size(nora.t,2),1) .*  (weights_Z_above + weights_Z_below), 2);
            pressure(pressure==0) = pressure_sl(pressure==0) .* exp(-g.*M.*(forcing.SPATIAL.STATVAR.altitude - altitude_sl)./R_gc./(Tair_sl(pressure==0) + 273.15));


            %scale Sin
            kd=0.952-1.041.*exp(-1.*exp(2.3-4.702.*transmissivity_point)); % Diffuse fraction.
            kd=max(kd,0);
            Sin_diff = kd.*Sin.*double(~sunset); % Diffuse shortwave radiation.
            Sin_dir = Sin - Sin_diff; % Direct component
            Sin_dir = S_TOA_point.*(Sin_dir./(max(1e-10, S_TOA_point))).^(pressure./pressure_sl);


            %scale Lin
            % Use the vapor pressure and temperature to calculate clear sky
            % emssivity at grid and subgrid. [also function]
            vp_sl = RH_sl .*satPresWater(proc, Tair_sl+273.15);
            vp_target = RH .*satPresWater(proc, Tair+273.15);
            x1=0.43; x2=5.7;
            cef=real(0.23+x1.*(vp_target./(Tair+273.15)).^(1/x2)); 
            cec=real(0.23+x1.*(vp_sl./(Tair_sl+273.15)).^(1/x2));
            
            % Diagnose the all sky emissivity at grid.
            aec=Lin./(proc.CONST.sigma.*(Tair_sl+273.15).^4);
            % Calculate the "cloud" emissivity at grid, assume this is the same at  subgrid.
            deltae=aec-cec;
            % Use the former cloud emissivity to compute the all sky emissivity at subgrid.
            aef=cef+deltae;
            Lin = aef.*proc.CONST.sigma.*(Tair + 273.15).^4;

            %Precipitation
            % threshold_precip = 0.1;
            %  Apply Liston & Elder (MicroMet) elevation-based precip adjustment
            dZ=forcing.SPATIAL.STATVAR.altitude - altitude_sl; % m
            dZ=dZ./1e3; % km
            dZ=min(dZ,2); % No larger that 2 km=3.3 adjustment
            dZ=max(dZ,-2);% For symmetry
            adjf=0.27; % Mean adjustment factor (following Fiddes and Gruber, 2014)
            adj=(1+adjf.*dZ)./(1-adjf.*dZ);
            precip = precip .*adj; 



% 
% 
% 
% 
% 
%             nora = forcing.TEMP.nora;
%             dist_lat = abs(forcing.SPATIAL.STATVAR.latitude - nora.lat);
%             dist_lon=abs(forcing.SPATIAL.STATVAR.longitude-nora.lon);
%             [dist_lat, ind_lat] = sort(dist_lat);
%             [dist_lon, ind_lon] = sort(dist_lon);
% 
%             dist_lat=dist_lat(1:2);
%             dist_lon=dist_lon(1:2);
% 
%             ind_lat = ind_lat(1:2);
%             ind_lon = ind_lon(1:2);
%             weights_lat = 1 - dist_lat./sum(dist_lat);
%             weights_lon = 1 - dist_lon./sum(dist_lon);
% 
%             nora_alt=double(nora.Z(ind_lon, ind_lat, :, :));
%             nora_alt=reshape(nora_alt, 4, size(nora.Z,3), size(nora.Z,4));
% 
%             nora_alt_sl  = double(nora.Zs(ind_lon, ind_lat));
%             nora_alt_sl = reshape(nora_alt_sl, 4,1);
%             nora_alt_sl = repmat(nora_alt_sl,1,1,size(nora_alt,3));
% 
%             nora_T=double(nora.T(ind_lon, ind_lat, :, :));
%             nora_T=reshape(nora_T, 4, size(nora.T,3), size(nora.T,4)).*nora.T_sf;
%             nora_u=double(nora.u(ind_lon, ind_lat, :, :)).* nora.wind_sf;
%             nora_u=reshape(nora_u, 4, size(nora.u,3), size(nora.u,4));
%             nora_v=double(nora.v(ind_lon, ind_lat, :, :)).* nora.wind_sf;
%             nora_v=reshape(nora_v, 4, size(nora.v,3), size(nora.v,4));
%             nora_RH=nora.RH_pl(ind_lon, ind_lat, :, :);
%             nora_RH=double(reshape(nora_RH, 4, size(nora.RH_pl,3), size(nora.RH_pl,4))).*nora.T_sf;
% %             nora_q=nora.q(ind_lon, ind_lat, :, :);
% %             nora_q=double(reshape(nora_q, 4, size(nora.q,3), size(nora.q,4))).*nora.q_sf;
% 
% 
%             nora_T_sl  = nora.T2(ind_lon, ind_lat, :);
%             nora_T_sl = reshape(double(nora_T_sl, 4,1, size(nora_T_sl,3))).*nora.T_sf;
%             nora_RH_sl  = nora.RH(ind_lon, ind_lat, :);
%             nora_RH_sl = reshape(double(nora_RH_sl, 4,1, size(nora_RH_sl,3))).*nora.T_sf;
%             nora_wind_sl = sqrt(double(nora.u10(ind_lon, ind_lat, :)).^2 + doubke(nora.v10(ind_lon, ind_lat, :)) .^2).* nora.wind_sf;
%             nora_wind_sl = reshape(nora_wind_sl, 4,1, size(nora_wind_sl,3));
%             nora_Lin_sl = double(nora.LW(ind_lon, ind_lat, :)).*nora.rad_sf;
%             nora_Lin_sl = reshape(nora_Lin_sl, 4,1, size(nora_Lin_sl,3));
%             nora_Sin_sl = double(nora.SW(ind_lon, ind_lat, :)).*nora.rad_sf;
%             nora_Sin_sl = reshape(nora_Sin_sl, 4,1, size(nora_Sin_sl,3));
% %             nora_S_TOA_sl = double(nora.S_TOA(ind_lon, ind_lat, :)).*nora.rad_sf;
% %             nora_S_TOA_sl = reshape(nora_S_TOA_sl, 4,1, size(nora_S_TOA_sl,3));
%             nora_precip_sl = double(nora.P(ind_lon, ind_lat, :)).*nora.P_sf;
%             nora_precip_sl = reshape(nora_precip_sl, 4,1, size(nora_precip_sl,3));
% 
%             %constants
%             R=287.05;  % Gas constant for dry air [JK^-1kg^-1]
%             g=9.81; % Accelnoration of gravity [ms^-1]
%             eps0=0.622; % Ratio of molecular weight of water and dry air [-]
%             S0=1370; % Solar constat (total TOA solar irradiance) [Wm^-2] used in ECMWF's IFS
% 
%             q2w=@(q) 0.5.*(1-sqrt(1-4.*q)); % Mixing ratio from specific humidity based on 2nd order Taylor series expansion.
%             wp2e=@(w,p) 0.5.*p.*(-1+sqrt(1+4.*w./eps0)); % Vapor pressure from mixing ratio based on 2nd order Taylor series expansion.
%             % AERK from Alduchov and Eskridge (1996).
%             A1=17.62; B1=243.12; C1=611.2; %A1=17.625; B1=243.04; C1=610.94;
%    %         Magnus=@(tc) C1.*exp(A1.*tc./(B1+tc)); % A version of the Magnus formula with the AERK parameters.
%             Magnus=@(tc) double(tc>=0) .* C1.*exp(A1.*tc./(B1+tc)) + double(tc<0) .* C1.*exp(22.46.*tc./(272.61+tc)); % A version of the Magnus formula with the AERK parameters.
% 
%             % Note, e=Magnus(tdc) and es=Magnus(tc)
% 
%             nora_p_sl = double(nora.ps(ind_lon, ind_lat, :)) .* nora.ps_sf;
%             nora_p_sl = reshape(nora_p_sl, 4,1, size(nora_p_sl,3));
%             nora_Td_sl = nora.Td2(ind_lon, ind_lat, :);
%             nora_Td_sl = double(reshape(nora_Td_sl, 4,1, size(nora_Td_sl,3))) .* nora.T_sf;
%             %vpsl=Magnus(K2C(nora_Td_sl));
%             vpsl=Magnus(nora_Td_sl);
% %             wsl=eps0.*vpsl./(nora_p_sl-vpsl);
% %             nora_q_sl=wsl./(1+wsl);
%             nora_q_sl = eps0.*vpsl./nora_p_sl;
% 
%             nora_RH_sl = Magnus(nora_Td_sl) ./ Magnus(double(nora_T_sl).*nora.T_sf);
% 
%             weights_lat = repmat(weights_lat', 2, 1, 1,size(nora_alt,3));
%             weights_lat=reshape(weights_lat, 4, 1 , size(nora_alt,3));
%             weights_lon = repmat(weights_lon, 1, 2, 1, size(nora_alt,3));
%             weights_lon=reshape(weights_lon, 4, 1 , size(nora_alt,3));
% 
%             layer_below = int16(nora_alt.*0);
%             layer_above = int16(nora_alt.*0);
% 
%             %do another one to get the lowermost pl above the orography
% 
%             for i=2:size(nora.Z,3)
%                 layer_below(:,i,:) = double(nora_alt(:,i,:) < forcing.SPATIAL.STATVAR.altitude & nora_alt(:,i-1,:) >= forcing.SPATIAL.STATVAR.altitude & nora_alt(:,i,:) > nora_alt_sl);
%                 layer_above(:,i-1,:) = double(nora_alt(:,i,:) < forcing.SPATIAL.STATVAR.altitude & nora_alt(:,i-1,:) >= forcing.SPATIAL.STATVAR.altitude & nora_alt(:,i-1,:) > nora_alt_sl);
%             end
% 
%             distance_Z_above = abs(sum(nora_alt .* layer_above ,2) - forcing.SPATIAL.STATVAR.altitude) .* double(sum(layer_above,2) > 0);
%             distance_Z_below = abs(sum(nora_alt .* layer_below ,2) - forcing.SPATIAL.STATVAR.altitude) .* double(sum(layer_below,2) > 0);
% 
%             weights_Z_above = 1-distance_Z_above ./ max(1e-10, distance_Z_above + distance_Z_below);
%             weights_Z_below = 1-distance_Z_below ./ max(1e-10, distance_Z_above + distance_Z_below);
%             weights_Z_above = repmat(weights_Z_above, 1, size(layer_above,2),1) .* double(layer_above);
%             weights_Z_below = repmat(weights_Z_below, 1, size(layer_below,2),1) .* double(layer_below);
% 
%             T_topoScale = sum(double(nora_T) .*  (weights_Z_above + weights_Z_below), 2);
%             wind_topoScale = sqrt(sum(double(nora_u) .*  (weights_Z_above + weights_Z_below), 2).^2 + sum(double(nora_v) .*  (weights_Z_above + weights_Z_below), 2).^2) ;
%             q_topoScale = sum(double(nora_q) .*  double(weights_Z_above + weights_Z_below), 2);
% 
%             use_sl = sum(weights_Z_above + weights_Z_below,2) <1-1e-9  | forcing.SPATIAL.STATVAR.altitude < nora_alt_sl;
%             merge_w_sl = sum(weights_Z_below,2) == 0;
%             factor = min(1, max(0, distance_Z_above./100));
% 
%             T_topoScale = double(~merge_w_sl) .* T_topoScale + double(merge_w_sl) .* factor .* T_topoScale + double(merge_w_sl) .* (1-factor) .* double(nora_T_sl);
%             T_topoScale = double(~use_sl) .* T_topoScale + double(use_sl) .* double(nora_T_sl);
%             wind_topoScale = double(~merge_w_sl) .* wind_topoScale + double(merge_w_sl) .* factor .* wind_topoScale + double(merge_w_sl) .* (1-factor) .* double(nora_wind_sl);
%             wind_topoScale = double(~use_sl) .* wind_topoScale + double(use_sl) .* double(nora_wind_sl);
%             q_topoScale = double(~merge_w_sl) .* q_topoScale + double(merge_w_sl) .* factor .* q_topoScale + double(merge_w_sl) .* (1-factor) .* double(nora_q_sl);
%             q_topoScale = double(~use_sl) .* q_topoScale + double(use_sl) .* double(nora_q_sl);
% 
% 
% 
%             T_topoScale = T_topoScale .* weights_lat;
%             T_topoScale = (T_topoScale(1:2,:,:) +T_topoScale(3:4,:,:));
%             wind_topoScale = wind_topoScale .* weights_lat;
%             wind_topoScale = (wind_topoScale(1:2,:,:) +wind_topoScale(3:4,:,:));
%             q_topoScale = q_topoScale .* double(weights_lat);
%             q_topoScale = (q_topoScale(1:2,:,:) + q_topoScale(3:4,:,:));
% 
%             RH_topoScale_sl = nora_RH_sl .* double(weights_lat);
%             RH_topoScale_sl = (RH_topoScale_sl(1:2,:,:) + RH_topoScale_sl(3:4,:,:));
% 
%             nora_Lin_sl = nora_Lin_sl .* weights_lat;
%             nora_Lin_sl = (nora_Lin_sl(1:2,:,:) + nora_Lin_sl(3:4,:,:));
%             nora_T_sl2 = double(nora_T_sl).*nora.T_sf .* weights_lat;
%             nora_T_sl2 = (nora_T_sl2(1:2,:,:) + nora_T_sl2(3:4,:,:));
%             nora_alt_sl2 = double(nora_alt_sl) .* weights_lat;
%             nora_alt_sl2 = (nora_alt_sl2(1:2,:,:) + nora_alt_sl2(3:4,:,:));
%             vp_sl2 = double(vpsl) .* weights_lat;
%             vp_sl2 = (vp_sl2(1:2,:,:) + vp_sl2(3:4,:,:));
%             nora_Sin_sl = nora_Sin_sl .* weights_lat;
%             nora_Sin_sl = (nora_Sin_sl(1:2,:,:) + nora_Sin_sl(3:4,:,:));
%             nora_S_TOA_sl = nora_S_TOA_sl .* weights_lat;
%             nora_S_TOA_sl = (nora_S_TOA_sl(1:2,:,:) + nora_S_TOA_sl(3:4,:,:));
%             nora_precip_sl = nora_precip_sl .* weights_lat;
%             nora_precip_sl = (nora_precip_sl(1:2,:,:) + nora_precip_sl(3:4,:,:));
% 
%             weights_lon = (weights_lon(1:2,:,:) + weights_lon(3:4,:,:))./2;
%             T_topoScale = squeeze(sum(T_topoScale .* weights_lon,1));
%             wind_topoScale = squeeze(sum(wind_topoScale .* weights_lon,1));
%             q_topoScale = squeeze(sum(q_topoScale .* double(weights_lon),1));
% 
%             RH_topoScale_sl = squeeze(sum(RH_topoScale_sl .* double(weights_lon),1));
% 
%             nora_Lin_sl = squeeze(sum(nora_Lin_sl .* weights_lon,1));
%             nora_alt_sl2 = squeeze(sum(nora_alt_sl2 .* weights_lon,1));
%             nora_T_sl2 = squeeze(sum(nora_T_sl2 .* weights_lon,1));
%             vp_sl2 = squeeze(sum(vp_sl2 .* weights_lon,1));
%             nora_Sin_sl = squeeze(sum(nora_Sin_sl .* weights_lon,1));
%             nora_S_TOA_sl = squeeze(sum(nora_S_TOA_sl .* weights_lon,1));            
%             nora_precip_sl = squeeze(sum(nora_precip_sl .* weights_lon,1));
% 
% 
% 
%             %pressure
%             for i=2:size(nora.Z,3)
%                 layer_below(:,i,:) = double(nora_alt(:,i,:) < forcing.SPATIAL.STATVAR.altitude & nora_alt(:,i-1,:) >= forcing.SPATIAL.STATVAR.altitude);
%                 layer_above(:,i-1,:) = double(nora_alt(:,i,:) < forcing.SPATIAL.STATVAR.altitude & nora_alt(:,i-1,:) >= forcing.SPATIAL.STATVAR.altitude);
%             end
% 
%             distance_Z_above = abs(sum(nora_alt .* layer_above ,2) - forcing.SPATIAL.STATVAR.altitude) .* double(sum(layer_above,2) > 0);
%             distance_Z_below = abs(sum(nora_alt .* layer_below ,2) - forcing.SPATIAL.STATVAR.altitude) .* double(sum(layer_below,2) > 0);
% 
%             weights_Z_above = 1-distance_Z_above ./ max(1e-10, distance_Z_above + distance_Z_below);
%             weights_Z_below = 1-distance_Z_below ./ max(1e-10, distance_Z_above + distance_Z_below);
%             weights_Z_above = repmat(weights_Z_above, 1, size(layer_above,2),1) .* double(layer_above);
%             weights_Z_below = repmat(weights_Z_below, 1, size(layer_below,2),1) .* double(layer_below);
% 
%             M = 0.0289644;
%             R_gc = 8.3144598;
%             p_topoScale = sum( repmat(nora.p,4,1,size(nora.t,2)) .*  (weights_Z_above + weights_Z_below), 2);
%             p_topoScale(p_topoScale==0) = nora_p_sl(p_topoScale==0) .* exp(-g.*M.*(forcing.SPATIAL.STATVAR.altitude - nora_alt_sl(p_topoScale==0))./R_gc./(double(nora_T_sl(p_topoScale==0)).*nora.T_sf + 273.15));
%             p_topoScale = p_topoScale .* weights_lat;
%             p_topoScale = (p_topoScale(1:2,:,:) + p_topoScale(3:4,:,:));
%             p_topoScale = squeeze(sum(p_topoScale .* weights_lon,1));
% 
%             %needed for Sin
%             nora_p_sl = nora_p_sl .* weights_lat;
%             nora_p_sl = (nora_p_sl(1:2,:,:) + nora_p_sl(3:4,:,:));
%             nora_p_sl = squeeze(sum(nora_p_sl .* weights_lon,1));
% 
% 
% 
%             %Long-wave
%             sbc=5.67e-8;
%             nora_Lin_sl(nora_Lin_sl<=0) = sbc .*(nora_T_sl2(nora_Lin_sl<=0) + 273.15).^4;
%             wf=q2w(q_topoScale); % Convert to mixing ratio at fine grid.
%             vpf=wp2e(wf,p_topoScale);
%             %disp(qout)
%             % Getting negative qout!
% 
%             % Use the vapor pressure and tempnorature to calculate clear sky
%             % emssivity at grid and subgrid. [also function]
%             x1=0.43; x2=5.7;
%             cef=real(0.23+x1.*(vpf./(T_topoScale+273.15)).^(1/x2)); % Obs! got negative vapor pressure-> imaginary number in LW calc
%             cec=real(0.23+x1.*(vp_sl2./(nora_T_sl2+273.15)).^(1/x2));
% 
%             % Diagnose the all sky emissivity at grid.
% 
%             aec=nora_Lin_sl./(sbc.*(nora_T_sl2+273.15).^4);
% 
%             % Calculate the "cloud" emissivity at grid, assume this is the same at
%             % subgrid.
%             deltae=aec-cec;
% 
%             % Use the former cloud emissivity to compute the all sky emissivity at subgrid.
%             aef=cef+deltae;
%             Lin_topoScale = aef.*sbc.*(T_topoScale+273.15).^4;
% 
% 
% 
%             %Sin
%             [solar_azm, solar_zen]=solargeom(proc, nora.t, forcing.SPATIAL.STATVAR.latitude, forcing.SPATIAL.STATVAR.longitude);
% 
% %             solar_azm=[0.5*(solar_azm(1:end-1)+solar_azm(2:end)) solar_azm(end)]';
% %             solar_zen=[0.5*(solar_zen(1:end-1)+solar_zen(2:end)) solar_zen(end)]';
% 
%             solar_azm = solar_azm';
%             solar_zen = solar_zen';
% 
%             % Compute downwelling TOA SW irradiance (i.e. the incoming shortwave
%             % incident on a horizontal plane at the TOA), by accounting for the
%             % solar zentih angle.
%             mu0=max(cos(solar_zen),0); % Trunacte negative values.
%             % May also want to consider treating values mu0<0 for prominent topography
%             % when the horizon  angles are less than 0.
%             sunset=mu0<cosd(89);%(mu0==0); % Sunset switch.
%             % Note, it would be better to use the true avnorage ((1/tau) integral_t^(t+tau) mu0 dt)
%             % but this approximation should be ok.
%             nora_S_TOA_point = S0.*mu0;
% 
%             kt=nora_Sin_sl./max(1e-10, nora_S_TOA_sl); % Clearness index.
%             kd=0.952-1.041.*exp(-1.*exp(2.3-4.702.*kt)); % Diffuse fraction.
%             kd=max(kd,0);
% 
%             % Diffuse component.
%             Sin_diff_topoScale = kd.*nora_Sin_sl.*double(~sunset); % Diffuse shortwave radiation.
%             % Direct component
%             Sin_dir_topoScale = nora_Sin_sl - Sin_diff_topoScale;
% 
%             Sin_dir_topoScale = nora_S_TOA_point.*(Sin_dir_topoScale./(max(1e-10, nora_S_TOA_point))).^(p_topoScale./nora_p_sl);
% 
%             % % Scale direct shortwave using Beer's law (see Aalstad 2019, Appendix A)
%             % ka=(g.*mu0./(psl)).*log(SWtoa./SWcdir); % Note, we don't get log(0) due to "if sunset" condition.
%             % SWfdir=SWtoa.*exp(-ka.*pout./(g*mu0));
% 
% 
%             %Precipitation
%             threshold_precip = 0.1;
%             %  Apply Liston & Elder (MicroMet) elevation-based precip adjustment
%             dZ=forcing.SPATIAL.STATVAR.altitude - nora_alt_sl2; % m
%             dZ=dZ./1e3; % km
%             dZ=min(dZ,2); % No larger that 2 km=3.3 adjustment
%             dZ=max(dZ,-2);% For symmetry
%             adjf=0.27; % Mean adjustment factor (following Fiddes and Gruber, 2014)
%             adj=(1+adjf.*dZ)./(1-adjf.*dZ);
% 
%             precip_topoScale = nora_precip_sl .*adj.*24; %in mm/day, check if timestep must be taken into account when not using 1h input data
            
            forcing.DATA.Tair = Tair;
            forcing.DATA.RH = RH;
            forcing.DATA.wind = wind;
            forcing.DATA.Sin_dir = Sin_dir;
            forcing.DATA.Sin_dif = Sin_diff;
            forcing.DATA.Sin =  forcing.DATA.Sin_dir + forcing.DATA.Sin_dif;
            forcing.DATA.S_TOA = S_TOA_point;
            forcing.DATA.Lin = Lin;
            forcing.DATA.p = pressure;
            forcing.DATA.precip = precip;
            forcing.DATA.timeForcing = nora.t;

            proc2 = calculate_q_from_RH();
            forcing = process(proc2, forcing, tile);
            
            forcing.TEMP.nora = [];
            
            %calculate q from RH

        end
        
        
        function p = satPresIce(proc, T)
            Tmfw = proc.CONST.Tmfw;
            p = 6.112.* 100.* exp(22.46.*(T-Tmfw)./(272.61+T-Tmfw));
        end
        
        function p = satPresWater(proc, T)
            Tmfw = proc.CONST.Tmfw;
            p = 6.112 .* 100 .* exp(17.62.*(T-Tmfw)./(243.12+T-Tmfw));
        end
        
%         function q = convertRelative2absoluteHumidity(forcing)
%             p=1005*100;
%             
%             q=(double((forcing.DATA.Tair)<=0).*(forcing.DATA.RH.*satPresIce(forcing, forcing.DATA.Tair+273.15)) + double(forcing.DATA.Tair>0).*(forcing.DATA.RH.*satPresWater(forcing, forcing.DATA.Tair+273.15)))./p;
%         end
        
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

