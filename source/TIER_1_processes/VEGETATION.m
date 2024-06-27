% %========================================================================
% CryoGrid TIER1 library class, functions related to vegetation
% R. B. Zweigel, June 2021
%========================================================================

classdef VEGETATION < BASE

    methods

        function canopy = PFT_PARA(canopy)
            % Assign PFT specific parameters from: https://www2.cesm.ucar.edu/models/cesm2/land/CLM50_Tech_Note.pdf
            if strcmp(canopy.PARA.PFT,'NDT boreal') || strcmp(canopy.PARA.PFT,'NET boreal')
                canopy.PARA.Khi_L           = 0.01;
                canopy.PARA.alpha_leaf_vis  = 0.07;
                canopy.PARA.alpha_leaf_nir  = 0.35;
                canopy.PARA.alpha_stem_vis  = 0.16;
                canopy.PARA.alpha_stem_nir  = 0.39;
                canopy.PARA.tau_leaf_vis    = 0.05;
                canopy.PARA.tau_leaf_nir    = 0.10;
                canopy.PARA.tau_stem_vis    = 0.001;
                canopy.PARA.tau_stem_nir    = 0.001;
                canopy.PARA.R_z0            = 0.055;
                canopy.PARA.R_d             = 0.67;
                canopy.PARA.d_leaf          = 0.04;
            elseif strcmp(canopy.PARA.PFT,'BDT boreal')
                canopy.PARA.Khi_L           = 0.25;
                canopy.PARA.alpha_leaf_vis  = 0.10;
                canopy.PARA.alpha_leaf_nir  = 0.45;
                canopy.PARA.alpha_stem_vis  = 0.16;
                canopy.PARA.alpha_stem_nir  = 0.39;
                canopy.PARA.tau_leaf_vis    = 0.05;
                canopy.PARA.tau_leaf_nir    = 0.25;
                canopy.PARA.tau_stem_vis    = 0.001;
                canopy.PARA.tau_stem_nir    = 0.001;
                canopy.PARA.R_z0            = 0.055;
                canopy.PARA.R_d             = 0.67;
                canopy.PARA.d_leaf          = 0.04;
            elseif strcmp(canopy.PARA.PFT,'BDS boreal')
                canopy.PARA.Khi_L           = 0.25;
                canopy.PARA.alpha_leaf_vis  = 0.10;
                canopy.PARA.alpha_leaf_nir  = 0.45;
                canopy.PARA.alpha_stem_vis  = 0.16;
                canopy.PARA.alpha_stem_nir  = 0.39;
                canopy.PARA.tau_leaf_vis    = 0.05;
                canopy.PARA.tau_leaf_nir    = 0.25;
                canopy.PARA.tau_stem_vis    = 0.001;
                canopy.PARA.tau_stem_nir    = 0.001;
                canopy.PARA.R_z0            = 0.120;
                canopy.PARA.R_d             = 0.68;
                canopy.PARA.d_leaf          = 0.04;
            elseif strcmp(canopy.PARA.PFT,'C3 grass')
                canopy.PARA.Khi_L           = -0.30;
                canopy.PARA.alpha_leaf_vis  = 0.11;
                canopy.PARA.alpha_leaf_nir  = 0.35;
                canopy.PARA.alpha_stem_vis  = 0.31;
                canopy.PARA.alpha_stem_nir  = 0.53;
                canopy.PARA.tau_leaf_vis    = 0.05;
                canopy.PARA.tau_leaf_nir    = 0.34;
                canopy.PARA.tau_stem_vis    = 0.12;
                canopy.PARA.tau_stem_nir    = 0.25;
                canopy.PARA.R_z0            = 0.12;
                canopy.PARA.R_d             = 0.68;
                canopy.PARA.d_leaf          = 0.04;
            else
                error('PFT-specific parameters must be assigned')
            end

        end
        
% ------- diagnose canopy properties ------
        function canopy = get_heat_capacity_leaves(canopy)
            % Heat capacities from Bonan 2018, assuming stems have same
            % heat capacity as leaves.
            L = canopy.STATVAR.LAI;
            S = canopy.STATVAR.SAI;
            SLA = canopy.PARA.SLA;
            f_c = canopy.PARA.f_carbon;
            f_w = canopy.PARA.f_water;
            c_w = canopy.CONST.c_w./canopy.CONST.rho_w; % [J/(kg*K)]
            c_dry = c_w ./ 3; % specific heat capacity of dry biomass
            
            Ma = 1/SLA; % leaf dry mass per unit area
            c_areal = c_dry.*Ma./f_c + c_w.*(f_w./(1-f_w)).*Ma./f_c;

            canopy.STATVAR.c_canopy = (L+S).*c_areal.*canopy.STATVAR.area(1); % [J/K]
        end

        function canopy = get_heat_capacity_canopy(canopy)
            % As get_get_heat_capacity_canopy_leaves(..), but with trunk heat
            % capacity from Swenson et al. (2018)
            kv = canopy.PARA.kv;
            D_bh = canopy.PARA.D_bh;
            h_tree = canopy.STATVAR.layerThick;
            N_tree = canopy.PARA.N_tree;
            rho_wood = canopy.PARA.rho_wood;
            L = canopy.STATVAR.LAI;
            SLA = canopy.PARA.SLA;
            f_c = canopy.PARA.f_carbon;
            f_w = canopy.PARA.f_water;
            c_w = canopy.CONST.c_w./canopy.CONST.rho_w; % [J/(kg*K)]
            c_dry = c_w ./ 3; % specific heat capacity of dry biomass
            
            Ma = 1/SLA; % leaf dry mass per unit area
            c_leaf_areal = c_dry.*Ma./f_c + c_w.*(f_w./(1-f_w)).*Ma./f_c;
            
            V_tree = kv.*pi.*(D_bh/2).^2.*h_tree; % trunk volume
            M_tree =  N_tree.*rho_wood.*V_tree; % tree dry mass per area 
            c_stem = ( c_dry + c_w.*(f_w./(1-f_w)) ).*M_tree;
            
            canopy.STATVAR.c_canopy = (L.*c_leaf_areal + c_stem).*canopy.STATVAR.area(1); % [J/K]
        end
        
        function canopy = get_z0_d_vegetation(canopy) % roughness length of vegetated surface (CLM5)
            L = canopy.STATVAR.LAI; % Leaf area index
            S = canopy.PARA.SAI; % Stem area index
            z0g = get_z0_surface(canopy.NEXT); % Roughness lenght of ground/snow surface
            R_z0 = canopy.PARA.R_z0; % Ratio of momentum roughness length to canopy height
            R_d = canopy.PARA.R_d; % Ratio of displacement height to canopy height
            z_top = sum(canopy.STATVAR.layerThick); % canopy height
            
            V = ( 1-exp(-1.*min(L+S,2))) ./ (1-exp(-2)); % Eq. 5.127
            z0 = exp( V.*log(z_top.*R_z0) + (1-V).*log(z0g) ); % Eq. 5.125
            d = z_top.*R_d.*V; % Eq. 5.126
            canopy.STATVAR.z0 = z0;
            canopy.STATVAR.d = d;
        end
        
        function canopy = get_E_water_vegetation(canopy)
            canopy.STATVAR.water = double(canopy.STATVAR.T>=0).*canopy.STATVAR.waterIce; % [m3]
            canopy.STATVAR.ice = double(canopy.STATVAR.T<0).*canopy.STATVAR.waterIce; % [m3]
            c_canopy = canopy.STATVAR.c_canopy; % [J/K]
            cp_waterIce = (canopy.STATVAR.water.*canopy.CONST.c_w + canopy.STATVAR.ice.*canopy.CONST.c_i); % [J/K]
            
            canopy.STATVAR.energy = canopy.STATVAR.T.*(c_canopy + cp_waterIce) - double(canopy.STATVAR.T<0).*canopy.STATVAR.waterIce.*canopy.CONST.L_f; % [J]
        end
        
        function canopy = get_T_water_vegetation(canopy)
            Lf = canopy.CONST.L_f; % [J/m3]
            c_i = canopy.CONST.c_i; % [J/m3/K]
            c_w = canopy.CONST.c_w; % [J/m3/K]
            c_canopy = canopy.STATVAR.c_canopy; % [J/K]
            L = canopy.STATVAR.LAI;
            S = canopy.STATVAR.SAI;
            W_max = canopy.PARA.Wmax; 
            
            e_frozen = -Lf.*canopy.STATVAR.waterIce;
            
            canopy.STATVAR.T = double(canopy.STATVAR.energy < e_frozen).*(canopy.STATVAR.energy - e_frozen)./(c_i.*canopy.STATVAR.waterIce + c_canopy) ...
                + double(canopy.STATVAR.energy > 0).* canopy.STATVAR.energy./(c_w.*canopy.STATVAR.waterIce + c_canopy);
            canopy.STATVAR.ice = double(canopy.STATVAR.energy <= e_frozen).*canopy.STATVAR.waterIce + double(canopy.STATVAR.energy > e_frozen & canopy.STATVAR.energy < 0) .* canopy.STATVAR.energy./(-Lf);
            canopy.STATVAR.water = double(canopy.STATVAR.energy >= 0).*canopy.STATVAR.waterIce + double(canopy.STATVAR.energy > e_frozen & canopy.STATVAR.energy < 0).*(canopy.STATVAR.energy - e_frozen)./Lf;
            
            canopy.STATVAR.f_wet = ( canopy.STATVAR.waterIce ./ (W_max.*canopy.STATVAR.area(1).*(L+S)) ).^(2/3);
            canopy.STATVAR.f_wet = min(1,canopy.STATVAR.f_wet);
            canopy.STATVAR.f_dry = (1-canopy.STATVAR.f_wet).*L./(L+S);
        end
        
        % NOT in use now
        % function canopy = get_T_simpleVegetatation(canopy)
        %     % Disregards water in canopy
        %     canopy.STATVAR.T = canopy.STATVAR.energy ./ canopy.STATVAR.c_canopy;
        % end
        % 
        % function canopy = get_E_simpleVegetation(canopy)
        %     canopy.STATVAR.energy = canopy.STATVAR.T .* canopy.STATVAR.c_canopy; % [J]
        % end

% -------- GET TIMESTEP FUNCTIONS ---------
        function timestep = get_timestep_canopy_T(canopy)
            d_energy = canopy.TEMP.d_energy;
            c_canopy =  canopy.STATVAR.c_canopy; % [J/m3]
            cp_waterIce = (canopy.STATVAR.water.*canopy.CONST.c_w + canopy.STATVAR.ice.*canopy.CONST.c_i); % [J/K]
            
            timestep = canopy.PARA.dT_max ./ ( abs(d_energy)./(c_canopy + cp_waterIce) );
        end

        function timestep = get_timestep_snow_vegetation(canopy)
            timestep = canopy.STATVAR.ice ./ -canopy.TEMP.d_snow;
            timestep(canopy.TEMP.d_snow == 0) = canopy.PARA.dt_max;
        end
        
        function timestep = get_timestep_water_vegetation(canopy)
            timestep(canopy.TEMP.d_water ~= 0) = double(canopy.TEMP.d_water < 0).* canopy.STATVAR.waterIce ./ -canopy.TEMP.d_water + ...
                double(canopy.TEMP.d_water > 0).* 0.1 .* canopy.PARA.Wmax.*canopy.STATVAR.area ./ canopy.TEMP.d_water;
            timestep(canopy.TEMP.d_water == 0) = canopy.PARA.dt_max;
            timestep(timestep<=0) = canopy.PARA.dt_max;
        end
        
% ------- modify canopy structure ------ % 
        function canopy = build_canopy(canopy) 
            canopy.STATVAR.emissivity = 1 - exp(-canopy.STATVAR.LAI-canopy.STATVAR.SAI); % my_bar = 1 for longwave!
            if isfield(canopy.PARA, 'heat_capacity_function')
                if ~isempty(canopy.PARA.heat_capacity_function)
                    heat_capacity_function = str2func(canopy.PARA.heat_capacity_function);
                    canopy = heat_capacity_function(canopy);
                else
                    canopy = get_heat_capacity_canopy(canopy);
                end
            else
                canopy = get_heat_capacity_canopy(canopy);
            end
            canopy = get_E_water_vegetation(canopy); % derive energy from temperature
            canopy = get_z0_d_vegetation(canopy);
        end
        
        function canopy = add_canopy(canopy)
            canopy.STATVAR.LAI = canopy.PARA.LAI;
            canopy.STATVAR.emissivity = 1 - exp(-canopy.STATVAR.LAI-canopy.STATVAR.SAI); % my_bar = 1 for longwave!
            if isfield(canopy.PARA, 'heat_capacity_function')
                if ~isempty(canopy.PARA.heat_capacity_function)
                    heat_capacity_function = str2func(canopy.PARA.heat_capacity_function);
                    canopy = heat_capacity_function(canopy);
                else
                    canopy = get_heat_capacity_canopy(canopy);
                end
            else
                canopy = get_heat_capacity_canopy(canopy);
            end
            canopy = get_E_water_vegetation(canopy); % derive energy from temperature
            canopy = get_z0_d_vegetation(canopy);
        end
        
        function canopy = add_canopy_linear(canopy, doy)
            doy = min(doy, canopy.PARA.t_leafsprout + canopy.PARA.leafsprout_period); % Avoid overshoot
            canopy.STATVAR.LAI = canopy.PARA.LAI.* max(canopy.PARA.leafsprout_period, doy - canopy.PARA.t_leafsprout + canopy.PARA.leafsprout_period) ./ canopy.PARA.leafsprout_period;
            canopy.STATVAR.emissivity = 1 - exp(-canopy.STATVAR.LAI-canopy.STATVAR.SAI);  % my_bar = 1 for longwave!
            if isfield(canopy.PARA, 'heat_capacity_function')
                if ~isempty(canopy.PARA.heat_capacity_function)
                    heat_capacity_function = str2func(canopy.PARA.heat_capacity_function);
                    canopy = heat_capacity_function(canopy);
                else
                    canopy = get_heat_capacity_canopy(canopy);
                end
            else
                canopy = get_heat_capacity_canopy(canopy);
            end
            canopy = get_E_water_vegetation(canopy); % derive energy from temperature
            canopy = get_z0_d_vegetation(canopy);
        end
        
        function canopy = remove_canopy(canopy)
            canopy.STATVAR.LAI = 0;
            canopy.STATVAR.emissivity = 1 - exp(-canopy.STATVAR.LAI-canopy.STATVAR.SAI); % my_bar = 1 for longwave!
            if isfield(canopy.PARA, 'heat_capacity_function')
                if ~isempty(canopy.PARA.heat_capacity_function)
                    heat_capacity_function = str2func(canopy.PARA.heat_capacity_function);
                    canopy = heat_capacity_function(canopy);
                else
                    canopy = get_heat_capacity_canopy(canopy);
                end
            else
                canopy = get_heat_capacity_canopy(canopy);
            end
            canopy = get_E_water_vegetation(canopy); % derive energy from temperature            
            canopy = get_z0_d_vegetation(canopy);
        end
        
        function canopy = remove_canopy_linear(canopy, doy)
            doy = min(doy,canopy.PARA.t_leaffall + canopy.PARA.leaffall_period);
            canopy.STATVAR.LAI = canopy.PARA.LAI.* max(0,canopy.PARA.t_leaffall - doy) ./ canopy.PARA.leaffall_period;
            canopy.STATVAR.emissivity = 1 - exp(-canopy.STATVAR.LAI-canopy.STATVAR.SAI); % my_bar = 1 for longwave!
            if isfield(canopy.PARA, 'heat_capacity_function')
                if ~isempty(canopy.PARA.heat_capacity_function)
                    heat_capacity_function = str2func(canopy.PARA.heat_capacity_function);
                    canopy = heat_capacity_function(canopy);
                else
                    canopy = get_heat_capacity_canopy(canopy);
                end
            else
                canopy = get_heat_capacity_canopy(canopy);
            end
            canopy = get_E_water_vegetation(canopy); % derive energy from temperature
        end

    end
end

