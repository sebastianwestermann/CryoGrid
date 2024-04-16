%========================================================================
% CryoGrid TIER1 INTERACTION (IA) class for functions related to surface
% energy balance extending through the top class (e.g. canopy/vegetation)
% R. B. Zweigel, August 2021
%========================================================================

classdef IA_SEB < IA_BASE

    methods

        % ---------- upper boundary conditions for surface ----------
        function ia_seb_water = get_boundary_condition_m_GROUND(ia_seb_water, tile)
            ia_seb_water = Q_h_CLM5_below_canopy(ia_seb_water, tile);
            ia_seb_water = Q_e_CLM5_below_canopy(ia_seb_water, tile);
            ia_seb_water = get_boundary_condition_m_RichardsEq_below_canopy(ia_seb_water, tile);

            % add fluxes to uppermost cell (ratiative fluxes are added by penetration)
            ia_seb_water.NEXT.TEMP.F_ub = (-ia_seb_water.NEXT.STATVAR.Qh - ia_seb_water.NEXT.STATVAR.Qe) .* ia_seb_water.NEXT.STATVAR.area(1);
            ia_seb_water.NEXT.TEMP.d_energy(1) = ia_seb_water.NEXT.TEMP.d_energy(1) + ia_seb_water.NEXT.TEMP.F_ub; % cant skip F_ub, it is required for the ground-snow IA_class
        end

        function ia_seb_water = get_boundary_condition_m_GROUND_Xice(ia_seb_water, tile)
            ia_seb_water = Q_h_CLM5_below_canopy (ia_seb_water, tile);
            ia_seb_water = Q_e_CLM5_below_canopy_Xice(ia_seb_water, tile);
            ia_seb_water = get_boundary_condition_m_RichardsEq_Xice_below_canopy(ia_seb_water, tile);

            % add fluxes to uppermost cell (ratiative fluxes are added by penetration)
            ia_seb_water.NEXT.TEMP.F_ub = (-ia_seb_water.NEXT.STATVAR.Qh - ia_seb_water.NEXT.STATVAR.Qe) .* ia_seb_water.NEXT.STATVAR.area(1);
            ia_seb_water.NEXT.TEMP.d_energy(1) = ia_seb_water.NEXT.TEMP.d_energy(1) + ia_seb_water.NEXT.TEMP.F_ub; % cant skip F_ub, it is required for the ground-snow IA_class
        end

        function snow = get_boundary_condition_m_create_CHILD(ia_seb_water, tile)
            % below-canopy equivalent of get_boundary_condition_u_create_CHILD in SNOW_crocus_bucketW_seb
            snow = ia_seb_water.NEXT.CHILD;
            forcing = tile.FORCING;
            ia_seb_water = get_boundary_condition_allSNOW_throughfall_m(ia_seb_water, tile); %add full snow, but rain only for snow-covered part

            snow_property_function = str2func(snow.PARA.snow_property_function);
            snow = snow_property_function(snow,forcing);

            snow.TEMP.F_ub = 0;
            snow.TEMP.F_lb = 0;
            snow.TEMP.F_ub_water = 0;
            snow.TEMP.F_lb_water = 0;
            snow.TEMP.F_ub_water_energy = 0;
            snow.TEMP.F_lb_water_energy = 0;
            snow.STATVAR.sublimation = 0;
            snow.TEMP.sublimation_energy = 0;
            snow.STATVAR.evaporation = 0;
            snow.TEMP.evaporation_energy = 0;
            snow.TEMP.rain_energy = 0;
            snow.TEMP.rainfall = 0;

            snow.TEMP.d_energy = 0;
            snow.TEMP.d_water = 0;
            snow.TEMP.d_water_energy = 0;

            snow.STATVAR.d = 0;
            snow.STATVAR.s = 0;
            snow.STATVAR.gs = 0;
            snow.STATVAR.time_snowfall = 0;
            %new
            snow.STATVAR.top_snow_date = forcing.TEMP.t;
            snow.STATVAR.bottom_snow_date = forcing.TEMP.t-0.1./24;

            snow.TEMP.metam_d_d = 0;
            snow.TEMP.wind_d_d = 0;
            snow.TEMP.metam_d_s = 0;
            snow.TEMP.wind_d_s = 0;
            snow.TEMP.metam_d_gs = 0;
            snow.TEMP.compact_d_D = 0;
            snow.TEMP.wind_d_D = 0;
            snow.TEMP.wind = forcing.TEMP.wind;
            snow.TEMP.wind_surface = ia_seb_water.PREVIOUS.STATVAR.u_star;

            %start with some non-zero values for area and layerThick
            snow.STATVAR.area = 1;
            snow.STATVAR.layerThick = 0.5 .* snow.PARA.swe_per_cell ./ (snow.TEMP.newSnow.STATVAR.density ./1000); %[m] layerThick adjusted so that always 0.5 .* snow.PARA.swe_per_cell is contained
            snow.STATVAR.energy = 0;
            snow.STATVAR.waterIce = 0;
            snow.STATVAR.ice = 0;
            snow.STATVAR.excessWater = 0;
            snow.STATVAR.upperPos = snow.PARENT.STATVAR.upperPos;
        end

        function ia_seb_water = get_boundary_condition_m_CHILD(ia_seb_water, tile)
            % below-canopy equivalent of get_boundary_condition_u_CHILD in SNOW_crocus_bucketW_seb
            forcing = tile.FORCING;
            snow = ia_seb_water.NEXT.CHILD;
            ia_seb_water = get_boundary_condition_allSNOW_throughfall_m(ia_seb_water, tile); %add full snow, but rain only for snow-covered part
            ia_seb_water = get_boundary_condition_m_SNOW_CHILD_below_canopy(ia_seb_water, tile);

            snow_property_function = str2func(snow.PARA.snow_property_function);
            snow = snow_property_function(snow,forcing);

            ia_seb_water = Q_h_CLM5_below_canopy_CHILD(ia_seb_water, tile);
            ia_seb_water = Q_e_CLM5_below_canopy_CHILD(ia_seb_water, tile);

            snow.TEMP.d_energy(1) = snow.TEMP.d_energy(1) + (-snow.STATVAR.Qh - snow.STATVAR.Qe).*snow.STATVAR.area(1);

            snow.TEMP.wind = forcing.TEMP.wind;
            snow.TEMP.wind_surface = ia_seb_water.PREVIOUS.STATVAR.u_star;
        end

        % ---------- Sensible heat flux for surface ---------
        function ia_seb_water = Q_h_CLM5_below_canopy(ia_seb_water, tile)
            % ground sensible heat flux as described in CLM5
            stratigraphy1 = ia_seb_water.PREVIOUS; % vegetation
            stratigraphy2 = ia_seb_water.NEXT; % ground/snow
            Tg = stratigraphy2.STATVAR.T(1)+tile.FORCING.CONST.Tmfw;
            Cp = stratigraphy1.CONST.cp;
            Ts = stratigraphy1.STATVAR.Ts+tile.FORCING.CONST.Tmfw; % canopy air temperature (K)
            rho_atm = stratigraphy1.TEMP.rho_atm; % moist air density
            r_ah_prime = stratigraphy1.TEMP.r_a_prime; % aerodynamic resistance to heat transfer soil - canopy air

            stratigraphy2.STATVAR.Qh  = -rho_atm.*Cp.*(Ts-Tg)./r_ah_prime; % Eq. 5.90
        end

        function ia_seb_water = Q_h_CLM5_below_canopy_CHILD(ia_seb_water, tile)
            % ground sensible heat flux as described in CLM5
            stratigraphy1 = ia_seb_water.PREVIOUS; % vegetation
            stratigraphy2 = ia_seb_water.NEXT.CHILD; % snow CHILD
            Tg = stratigraphy2.STATVAR.T(1)+tile.FORCING.CONST.Tmfw;
            Cp = stratigraphy1.CONST.cp;
            Ts = stratigraphy1.STATVAR.Ts+tile.FORCING.CONST.Tmfw; % canopy air temperature (K)
            rho_atm = stratigraphy1.TEMP.rho_atm; % moist air density
            r_ah_prime = stratigraphy1.TEMP.r_a_prime; % aerodynamic resistance to heat transfer soil - canopy air

            stratigraphy2.STATVAR.Qh  = -rho_atm.*Cp.*(Ts-Tg)./r_ah_prime; % Eq. 5.90
        end

        % ---------- Latent heat flux for surface ---------
        function ia_seb_water = Q_e_CLM5_below_canopy(ia_seb_water, tile)
            % ground sensible heat flux as described in CLM5
            stratigraphy1 = ia_seb_water.PREVIOUS; % vegetation
            stratigraphy2 = ia_seb_water.NEXT; % ground
            q_g = stratigraphy2.STATVAR.q_g; % ground specific humidity
            q_s = stratigraphy1.STATVAR.q_s; % canopy air specific humidity
            Tg = stratigraphy2.STATVAR.T(1);
            Tmfw = tile.FORCING.CONST.Tmfw;
            r_aw_prime = stratigraphy1.TEMP.r_a_prime;
            r_soil = stratigraphy1.TEMP.r_soil;
            rho_atm = stratigraphy1.TEMP.rho_atm; % moist air density

            water_fraction = stratigraphy2.STATVAR.water(1) ./ stratigraphy2.STATVAR.waterIce(1);
            ice_fraction = stratigraphy2.STATVAR.ice(1) ./ stratigraphy2.STATVAR.waterIce(1);
            latent_heat = water_fraction.*latent_heat_evaporation(stratigraphy2, Tg+Tmfw) + ice_fraction.*latent_heat_sublimation(stratigraphy2, Tg+Tmfw);

            stratigraphy2.STATVAR.Qe = -rho_atm.*latent_heat.*(q_s-q_g)./(r_aw_prime+r_soil);

            stratigraphy2.STATVAR.evaporation = -water_fraction .* stratigraphy2.STATVAR.Qe ./ (latent_heat .* stratigraphy2.CONST.rho_w);
            stratigraphy2.STATVAR.sublimation = -ice_fraction .* stratigraphy2.STATVAR.Qe ./ (latent_heat .* stratigraphy2.CONST.rho_w);
            stratigraphy2.TEMP.evaporation_energy =  stratigraphy2.STATVAR.evaporation.*  (double(Tg>=0) .* stratigraphy2.CONST.c_w .* Tg + ...
                double(Tg<0) .* stratigraphy2.CONST.c_i .* Tg);
            stratigraphy2.TEMP.sublimation_energy =  stratigraphy2.STATVAR.sublimation .* (stratigraphy2.CONST.c_i .* stratigraphy2.STATVAR.T(1,1) - stratigraphy2.CONST.L_f);

            stratigraphy2.TEMP.d_water_ET(1) = stratigraphy2.TEMP.d_water_ET(1) + stratigraphy2.STATVAR.evaporation.*stratigraphy2.STATVAR.area(1);
            stratigraphy2.TEMP.d_water_ET_energy(1) = stratigraphy2.TEMP.d_water_ET_energy(1) + stratigraphy2.TEMP.evaporation_energy.*stratigraphy2.STATVAR.area(1);

        end

        function ia_seb_water = Q_e_CLM5_below_canopy_Xice(ia_seb_water, tile)
            % ground sensible heat flux as described in CLM5
            stratigraphy1 = ia_seb_water.PREVIOUS; % vegetation
            stratigraphy2 = ia_seb_water.NEXT; % ground
            q_g = stratigraphy2.STATVAR.q_g; % ground specific humidity
            q_s = stratigraphy1.STATVAR.q_s; % canopy air specific humidity
            Tg = stratigraphy2.STATVAR.T(1);
            Tmfw = tile.FORCING.CONST.Tmfw;
            r_aw_prime = stratigraphy1.TEMP.r_a_prime;
            r_soil = stratigraphy1.TEMP.r_soil;
            rho_atm = stratigraphy1.TEMP.rho_atm; % moist air density

            water_fraction = (stratigraphy2.STATVAR.water(1,1)+stratigraphy2.STATVAR.Xwater(1,1) ) ./ (stratigraphy2.STATVAR.waterIce(1,1) + stratigraphy2.STATVAR.XwaterIce(1,1));
            ice_fraction = (stratigraphy2.STATVAR.ice(1,1) + stratigraphy2.STATVAR.Xice(1,1)) ./ (stratigraphy2.STATVAR.waterIce(1,1) + stratigraphy2.STATVAR.XwaterIce(1,1));
            water_fraction(stratigraphy2.STATVAR.waterIce == 0) = double(Tg>=0);
            ice_fraction(stratigraphy2.STATVAR.waterIce == 0) = double(Tg<0);
            latent_heat = water_fraction.*latent_heat_evaporation(stratigraphy2, Tg+Tmfw) + ice_fraction.*latent_heat_sublimation(stratigraphy2, Tg+Tmfw);

            Qe_evap = -rho_atm.*latent_heat_evaporation(stratigraphy2, Tg+Tmfw).*(q_s-q_g)./(r_aw_prime+r_soil) .* water_fraction;
            Qe_sublim = -rho_atm.*latent_heat_sublimation(stratigraphy2, Tg+Tmfw).*(q_s-q_g)./(r_aw_prime+r_soil) .* ice_fraction;

            stratigraphy2.STATVAR.Qe = Qe_evap + Qe_sublim;

            stratigraphy2.STATVAR.evaporation = -Qe_evap ./ (latent_heat .* stratigraphy2.CONST.rho_w);
            stratigraphy2.STATVAR.sublimation = -Qe_sublim ./ (latent_heat .* stratigraphy2.CONST.rho_w);
            stratigraphy2.TEMP.evaporation_energy =  stratigraphy2.STATVAR.evaporation.*  (double(Tg>=0) .* stratigraphy2.CONST.c_w .* Tg + ...
                double(Tg<0) .* stratigraphy2.CONST.c_i .* Tg);
            stratigraphy2.TEMP.sublimation_energy =  stratigraphy2.STATVAR.sublimation .* (stratigraphy2.CONST.c_i .* Tg - stratigraphy2.CONST.L_f);

            stratigraphy2.TEMP.d_water_ET(1) = stratigraphy2.TEMP.d_water_ET(1) + stratigraphy2.STATVAR.evaporation.*stratigraphy2.STATVAR.area(1);
            stratigraphy2.TEMP.d_water_ET_energy(1) = stratigraphy2.TEMP.d_water_ET_energy(1) + stratigraphy2.TEMP.evaporation_energy.*stratigraphy2.STATVAR.area(1);
        end

        function ia_seb_water = Q_e_CLM5_SNOW_below_canopy(ia_seb_water, tile)
            % below-canopy equivalent to Q_eq_potET_snow in SNOW
            stratigraphy1 = ia_seb_water.PREVIOUS; % vegetation
            stratigraphy2 = ia_seb_water.NEXT; % snow
            q_g = stratigraphy2.STATVAR.q_g; % snow specific humidity
            q_s = stratigraphy1.STATVAR.q_s; % canopy air specific humidity
            Tsurf = stratigraphy2.STATVAR.T(1);
            Tmfw = tile.FORCING.CONST.Tmfw;
            r_aw_prime = stratigraphy1.TEMP.r_a_prime;
            rho_atm = stratigraphy1.TEMP.rho_atm; % moist air density

            porespace = stratigraphy2.STATVAR.layerThick(1).*stratigraphy2.STATVAR.area(1) - stratigraphy2.STATVAR.ice(1);
            saturation = stratigraphy2.STATVAR.water(1)./porespace;
            saturation(porespace == 0) = 0;
            L_w = latent_heat_evaporation(stratigraphy2, Tsurf+Tmfw);
            L_i = latent_heat_sublimation(stratigraphy2, Tsurf+Tmfw);

            if Tsurf < 0 || saturation <= 0 % All sublimation
                stratigraphy2.STATVAR.Qe = -rho_atm.*L_i.*(q_s-q_g)./r_aw_prime;
                stratigraphy2.STATVAR.sublimation = - stratigraphy2.STATVAR.Qe ./ (L_i .* stratigraphy2.CONST.rho_w) .*stratigraphy2.STATVAR.area(1);
                stratigraphy2.TEMP.sublimation_energy =  stratigraphy2.STATVAR.sublimation .* (stratigraphy2.CONST.c_i .* stratigraphy2.STATVAR.T(1,1) - stratigraphy2.CONST.L_f);
                stratigraphy2.STATVAR.evaporation = 0;
                stratigraphy2.TEMP.evaporation_energy = 0;
            else
                evap_fraction = max(0, min(saturation./(2.*stratigraphy2.PARA.field_capacity),1)); % all evap when saturation is more than twice of field capacity (in fraction of porespace)
                sublim_fraction = 1 - evap_fraction;
                Qe_sublim = -sublim_fraction.*rho_atm.*L_i.*(q_s-q_g)./r_aw_prime;
                Qe_evap = -evap_fraction.*rho_atm.*L_w.*(q_s-q_g)./r_aw_prime;
                stratigraphy2.STATVAR.sublimation = - Qe_sublim ./(stratigraphy2.CONST.rho_w .* L_i) .* stratigraphy2.STATVAR.area(1);
                stratigraphy2.TEMP.sublimation_energy = stratigraphy2.STATVAR.sublimation .* (Tsurf .* stratigraphy2.CONST.c_i - stratigraphy2.CONST.L_f);
                stratigraphy2.STATVAR.evaporation = - Qe_evap ./(stratigraphy2.CONST.rho_w .* L_w) .* stratigraphy2.STATVAR.area(1);
                stratigraphy2.TEMP.evaporation_energy = stratigraphy2.STATVAR.evaporation .* stratigraphy2.CONST.c_w;

                stratigraphy2.STATVAR.Qe = Qe_evap + Qe_sublim;
            end
            % Note: no d_water_ET for SNOW, evap/sublimation is added to water fluxes in advance_prognostics(..)
        end

        function ia_seb_water = Q_e_CLM5_below_canopy_CHILD(ia_seb_water, tile)
            % below-canopy equivalent to Q_eq_potET_snow in SNOW for SNOW in CHILD phase
            stratigraphy1 = ia_seb_water.PREVIOUS; % vegetation
            stratigraphy2 = ia_seb_water.NEXT.CHILD; % snow CHILD
            q_g = stratigraphy2.STATVAR.q_g; % ground specific humidity
            q_s = stratigraphy1.STATVAR.q_s; % canopy air specific humidity
            Tsurf = stratigraphy2.STATVAR.T(1);
            Tmfw = tile.FORCING.CONST.Tmfw;
            r_aw_prime = stratigraphy1.TEMP.r_a_prime;
            rho_atm = stratigraphy1.TEMP.rho_atm; % moist air density

            porespace = stratigraphy2.STATVAR.layerThick(1).*stratigraphy2.STATVAR.area(1) - stratigraphy2.STATVAR.ice(1);
            saturation = stratigraphy2.STATVAR.water(1)./porespace;
            saturation(porespace == 0) = 0;
            L_w = latent_heat_evaporation(stratigraphy2, Tsurf+Tmfw);
            L_i = latent_heat_sublimation(stratigraphy2, Tsurf+Tmfw);

            if Tsurf < 0 || saturation <= 0 % All sublimation
                stratigraphy2.STATVAR.Qe = -rho_atm.*L_i.*(q_s-q_g)./r_aw_prime;
                stratigraphy2.STATVAR.sublimation = - stratigraphy2.STATVAR.Qe ./ (L_i .* stratigraphy2.CONST.rho_w) .*stratigraphy2.STATVAR.area(1);
                stratigraphy2.TEMP.sublimation_energy =  stratigraphy2.STATVAR.sublimation .* (stratigraphy2.CONST.c_i .* stratigraphy2.STATVAR.T(1,1) - stratigraphy2.CONST.L_f);
                stratigraphy2.STATVAR.evaporation = 0;
                stratigraphy2.TEMP.evaporation_energy = 0;
            else
                evap_fraction = max(0, min(saturation./(2.*stratigraphy2.PARA.field_capacity),1)); % all evap when saturation is more than twice of field capacity (in fraction of porespace)
                sublim_fraction = 1 - evap_fraction;
                Qe_sublim = -sublim_fraction.*rho_atm.*L_i.*(q_s-q_g)./r_aw_prime;
                Qe_evap = -evap_fraction.*rho_atm.*L_w.*(q_s-q_g)./r_aw_prime;
                stratigraphy2.STATVAR.sublimation = - Qe_sublim ./(stratigraphy2.CONST.rho_w .* L_i) .* stratigraphy2.STATVAR.area(1);
                stratigraphy2.TEMP.sublimation_energy = stratigraphy2.STATVAR.sublimation .* (Tsurf .* stratigraphy2.CONST.c_i - stratigraphy2.CONST.L_f);
                stratigraphy2.STATVAR.evaporation = - Qe_evap ./(stratigraphy2.CONST.rho_w .* L_w) .* stratigraphy2.STATVAR.area(1);
                stratigraphy2.TEMP.evaporation_energy = stratigraphy2.STATVAR.evaporation .* stratigraphy2.CONST.c_w;

                stratigraphy2.STATVAR.Qe = Qe_evap + Qe_sublim; % Note this is not like in Q_eq_potET_snow! RBZ
            end
        end

    end
end
