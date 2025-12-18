%========================================================================
% CryoGrid TIER1 library class for functions related to air convection
% CAUTION: this is highly experimental and should not be used!!
% The theory will likely not withstand a peer-review process!
% S. Westermann, October 2020
%========================================================================


classdef AIR_CONVECTION < BASE
    
    methods
        
        function ground = get_boundary_condition_u_convection(ground, forcing)
            
            rho_air = ground.PARA.pressure ./ (ground.CONST.R_spec.*([forcing.TEMP.Tair; ground.STATVAR.T(1,1)] + ground.CONST.Tmfw));
            
            delta_p = double(rho_air(1,1) > rho_air(2,1)) .* (rho_air(1,1) - rho_air(2,1)) ./ 4 .* ground.CONST.g ;

            %rough pipe turbulent regime, constant Darcy friction factor
            velocity_turbulent = sqrt(2 .* ground.STATVAR.diamater_pipe(1,1) ./ mean(rho_air,1) ./ground.CONST.Darcy_friction_factor .* delta_p ./ (ground.STATVAR.layerThick(1,1)./2)); %Darcy Weisbach   

            %laminar flow
            velocity_laminar = ground.STATVAR.diamater_pipe(1,1).^2 ./ (32.*ground.CONST.viscosity_air) .* delta_p ./ (ground.STATVAR.layerThick(1,1)./2);
           
            air_flux = pi() ./ 4 .* ground.STATVAR.diamater_pipe(1,1).^2 .* min(velocity_turbulent, velocity_laminar)  .* ground.STATVAR.number_of_pipes(1,1);        
     
            energy_convection_down =  air_flux.* ground.CONST.cp .*(rho_air(1,1) + rho_air(2,1))./2 .* forcing.TEMP.Tair .* ground.STATVAR.area(1,1);
            energy_convection_up = air_flux.* ground.CONST.cp .* (rho_air(1,1) + rho_air(2,1))./2 .* ground.STATVAR.T(1,1) .* ground.STATVAR.area(1,1);
            
            ground.TEMP.d_energy(1,1) = ground.TEMP.d_energy(1,1) - energy_convection_up + energy_convection_down;
        end
        
        %-----derivatives----------
        
        function ground = get_derivative_air_convection_Darcy(ground) %dont use, this assumes laminar flow and overestimates flow at high block sizes

            permeability_air = ground.STATVAR.permeability_air(1:end-1,1) .* ground.STATVAR.permeability_air(2:end,1) ./ ...
                (ground.STATVAR.permeability_air(1:end-1,1) .* ground.STATVAR.layerThick(2:end,1) ./2 + ground.STATVAR.permeability_air(2:end,1) .* ground.STATVAR.layerThick(1:end-1,1) ./ 2);

            permeability_air(isnan(permeability_air)) = 0;
            
            rho_air = ground.PARA.pressure ./ (ground.CONST.R_spec.*(ground.STATVAR.T + ground.CONST.Tmfw));
            
            delta_p = double(rho_air(1:end-1,1) > rho_air(2:end,1)) .* (rho_air(1:end-1,1) - (rho_air(1:end-1,1)+rho_air(2:end,1))./2) .* ground.STATVAR.layerThick(2:end) ./ 2 .* ground.CONST.g ;

            air_flux =  permeability_air ./ ground.CONST.viscosity_air .* delta_p; %m3/sec m2
                        
            energy_convection_down =  air_flux.* ground.CONST.cp .*(rho_air(1:end-1,1)+rho_air(2:end,1))./2 .* ground.STATVAR.T(1:end-1,1) .* ground.STATVAR.area(1:end-1,1);
            energy_convection_up = air_flux.* ground.CONST.cp .* (rho_air(1:end-1,1)+rho_air(2:end,1))./2 .* ground.STATVAR.T(2:end,1) .* ground.STATVAR.area(2:end,1);
            
            ground.TEMP.d_energy(1:end-1,1) = ground.TEMP.d_energy(1:end-1,1) + energy_convection_up - energy_convection_down;
            ground.TEMP.d_energy(2:end,1) = ground.TEMP.d_energy(2:end,1) - energy_convection_up + energy_convection_down;
        end
        
                
        function ground = get_derivative_air_convection_Darcy_Weisbach(ground)
            
            rho_air = ground.PARA.pressure ./ (ground.CONST.R_spec.*(ground.STATVAR.T + ground.CONST.Tmfw));
            
            delta_p = double(rho_air(1:end-1,1) > rho_air(2:end,1)) .* (rho_air(1:end-1,1) - rho_air(2:end,1)) .* ground.STATVAR.layerThick(1:end-1,1).* ground.STATVAR.layerThick(2:end,1) ./ (ground.STATVAR.layerThick(1:end-1,1) + ground.STATVAR.layerThick(2:end,1)).^2 .* ground.CONST.g; 
%this line makes no sense, leads to /4 when both cells are equal - why?

            %rough pipe turbulent regime, constant Darcy friction factor
            delta_p_upperCell = ground.STATVAR.diamater_pipe(2:end,1).^5 .* ground.STATVAR.number_of_pipes(2:end,1).^2 .* delta_p ./ ground.STATVAR.layerThick(2:end,1) .* ...
                (ground.STATVAR.diamater_pipe(2:end,1).^5 .* ground.STATVAR.number_of_pipes(2:end,1).^2 ./ ground.STATVAR.layerThick(2:end,1) +  ground.STATVAR.diamater_pipe(1:end-1,1).^5 .* ground.STATVAR.number_of_pipes(1:end-1,1).^2 ./ ground.STATVAR.layerThick(1:end-1,1)).^-1;
            delta_p_upperCell(isnan(delta_p_upperCell)) = 0;
%somehow does weighing analogous to series of resistor, but  bit unsure why it is done like that -> attempts to construct a permeability and then combine for different cells            
            
            velocity_upper_cell_turbulent = sqrt(2 .* ground.STATVAR.diamater_pipe(1:end-1,1) ./ mean(rho_air,1) ./ground.CONST.Darcy_friction_factor .* delta_p_upperCell ./ (ground.STATVAR.layerThick(1:end-1,1)./2));%Darcy Weisbach   

            %laminar flow
            delta_p_upperCell = ground.STATVAR.diamater_pipe(2:end,1).^4 .* ground.STATVAR.number_of_pipes(2:end,1) .* delta_p ./ ground.STATVAR.layerThick(2:end,1) .* ...
                (ground.STATVAR.diamater_pipe(2:end,1).^4 .* ground.STATVAR.number_of_pipes(2:end,1) ./ ground.STATVAR.layerThick(2:end,1) +  ground.STATVAR.diamater_pipe(1:end-1,1).^4 .* ground.STATVAR.number_of_pipes(1:end-1,1) ./ ground.STATVAR.layerThick(1:end-1,1)).^-1;
            delta_p_upperCell(isnan(delta_p_upperCell)) = 0;
            
            velocity_upper_cell_laminar = ground.STATVAR.diamater_pipe(1:end-1,1).^2 ./ (32.*ground.CONST.viscosity_air) .* delta_p_upperCell ./ (ground.STATVAR.layerThick(1:end-1,1)./2);
            
            
            ground.TEMP.air_flux = pi() ./ 4 .* ground.STATVAR.diamater_pipe(1:end-1).^2 .* min(velocity_upper_cell_turbulent, velocity_upper_cell_laminar)  .* ground.STATVAR.number_of_pipes(1:end-1);        
            
            energy_convection_down =  ground.TEMP.air_flux.* ground.CONST.cp .*(rho_air(1:end-1,1)+rho_air(2:end,1))./2 .* ground.STATVAR.T(1:end-1,1) .* ground.STATVAR.area(1:end-1,1);
            energy_convection_up = ground.TEMP.air_flux.* ground.CONST.cp .* (rho_air(1:end-1,1)+rho_air(2:end,1))./2 .* ground.STATVAR.T(2:end,1) .* ground.STATVAR.area(2:end,1);
            
            ground.TEMP.d_energy(1:end-1,1) = ground.TEMP.d_energy(1:end-1,1) + energy_convection_up - energy_convection_down;
            ground.TEMP.d_energy(2:end,1) = ground.TEMP.d_energy(2:end,1) - energy_convection_up + energy_convection_down;
        end
        
        
        function timestep = get_timestep_air_convection(ground) %leads to super-small timesteps
            air_flux_nonzero = ground.TEMP.air_flux > 0;
            timestep_up = min((ground.STATVAR.layerThick(air_flux_nonzero) ./ ground.STATVAR.area(air_flux_nonzero) - ground.STATVAR.waterIce(air_flux_nonzero) - ground.STATVAR.mineral(air_flux_nonzero) - ground.STATVAR.organic(air_flux_nonzero)) ./ 2 ./ ground.TEMP.air_flux(air_flux_nonzero)); %no more than half the void space up

            timestep_down  = (ground.STATVAR.layerThick(1:end-1,1) ./ ground.STATVAR.area(1:end-1,1) - ground.STATVAR.waterIce(1:end-1,1) - ground.STATVAR.mineral(1:end-1,1) - ground.STATVAR.organic(1:end-1,1)) ./ 2 ./ ground.TEMP.air_flux(2:end,1);
            
            timestep_down = nanmin(timestep_down);
            timestep = min(timestep_up, timestep_down);
            timestep = max(1, timestep);
        end
        
        
        %---permeability air--------------
        function ground = pipes_Darcy_Weisbach(ground)
            ground.STATVAR.porosity = 1 - (ground.STATVAR.waterIce + ground.STATVAR.mineral + ground.STATVAR.organic) ./ ground.STATVAR.layerThick ./ ground.STATVAR.area; 
            
            ground.STATVAR.diamater_pipe = 2./3 .*ground.STATVAR.porosity ./(1-ground.STATVAR.porosity) .* ground.STATVAR.grain_size;
            ground.STATVAR.number_of_pipes = ground.STATVAR.porosity ./ (pi()./4 .* ground.STATVAR.diamater_pipe.^2 .* ground.CONST.tortuosity_air); 
            ground.STATVAR.number_of_pipes(isnan(ground.STATVAR.number_of_pipes)) = 0; %if diameter is zero
        end
        
        function ground = permeability_air_Carman_Kozeny(ground)
            porosity = 1 - (ground.STATVAR.waterIce + ground.STATVAR.mineral + ground.STATVAR.organic) ./ ground.STATVAR.layerThick ./ ground.STATVAR.area;
            ground.STATVAR.permeability_air =  ground.STATVAR.grain_size.^2 ./ 180 .* porosity .^3 ./ (1 - porosity).^2;
            ground.STATVAR.porosity = porosity;
        end
        
        function ground = permeability_air_Rumpf_Gupte(ground)
            porosity = 1 - (ground.STATVAR.waterIce + ground.STATVAR.mineral + ground.STATVAR.organic) ./ ground.STATVAR.layerThick ./ ground.STATVAR.area;
            ground.STATVAR.permeability_air = grain_size.^2 ./ 5.6 .* porosity .^5.5 ./ 1.05; %Rumpf-Gupte model [m2]
        end



        %new functions based on Darcy-Forchheimer equation, S. Westermann Dec 2025
        function ground = permeability_air_Carman_Kozeny2(ground)

            porosity = 1 - (ground.STATVAR.waterIce + ground.STATVAR.mineral + ground.STATVAR.organic) ./ ground.STATVAR.layerThick ./ ground.STATVAR.area;
            ground.STATVAR.porosity = porosity;

            %permeability for each grid cell
            ground.STATVAR.permeability_air =  ground.STATVAR.effective_grain_size.^2 ./ ground.PARA.CozenyKarman_constant .* porosity .^3 ./ (1 - porosity).^2;
            
            %combine permeabilities  
            ground.TEMP.permeability_air = ground.STATVAR.permeability_air(1:end-1,1) .* ground.STATVAR.permeability_air(2:end,1).* (ground.STATVAR.layerThick(1:end-1,1) + ground.STATVAR.layerThick(2:end,1)) ./ ...
                (ground.STATVAR.permeability_air(1:end-1,1) .* ground.STATVAR.layerThick(2:end,1) + ground.STATVAR.permeability_air(2:end,1).* ground.STATVAR.layerThick(1:end-1,1));
        end

        function ground = Forchheimer_facor(ground)
            %calculate average porosity
            ground.TEMP.porosity =  (ground.STATVAR.porosity(1:end-1,1) + ground.STATVAR.porosity(2:end,1))./2;
            ground.TEMP.Forchheimer_factor = ground.PARA.Forchheimer_constant./sqrt(ground.PARA.CozenyKarman_constant .* ground.TEMP.porosity.^3);
        end

        function ground = calculate_air_density_viscosity_convection(ground, tile)
            %calculate average T
            T = (ground.STATVAR.T(1:end-1,1) + ground.STATVAR.T(2:end,1))./2;
            T= T+273.15;
            %dry air density
            ground.TEMP.density_air = tile.FORCING.TEMP.p./(287.058.*T); 
            %Sutherlands formula
            ground.TEMP.viscosity_air = 1.72e-5 .* (T./273.15).^1.5 .* (273.15+110.4)./(T+110.4);
        end

        function ground = calculate_velocity_convection(ground)

            %using Darcy-Forchheimer equation as in
            %https://www.mdpi.com/2227-7390/12/23/3653, omitting Brinkman term
            delta_T = (ground.STATVAR.T(1:end-1,1) - ground.STATVAR.T(2:end,1)); 
            delta_T = double(delta_T<0) .* delta_T;

            %factors a,b,c in quadratic equation
            a=-ground.TEMP.Forchheimer_factor./sqrt(ground.TEMP.permeability_air);
            b = -ground.TEMP.viscosity_air ./ ground.TEMP.permeability_air ./ ground.TEMP.density_air;
            thermal_expansion_coefficient_air = (273.15+(ground.STATVAR.T(1:end-1,1) + ground.STATVAR.T(2:end,1))./2).^-1;
            c = -thermal_expansion_coefficient_air .* delta_T/2 .* ground.CONST.g; %the 2 is because it is T-T_mean

            p=b./a .* double(delta_T<0); %last condition set to avoid numerical problems
            q=c./a;

            ground.STATVAR.velocity = - p./2 + real(sqrt((p./2).^2 - q)); %velocity [m/ec] is positive
            ground.STATVAR.velocity(isnan(ground.STATVAR.velocity)) = 0;

        end

        function ground = get_derivative_air_convection_DarcyForchheimer(ground, tile) 
            % 
            % if tile.t > datenum(2013,4,20,8,0,0)
            % 
            %     a=1;
            % end
            
            air_flux =  0.5 .* ground.STATVAR.velocity .*  ground.TEMP.porosity .* (ground.STATVAR.area(1:end-1,1)+ground.STATVAR.area(2:end,1))./2; %covering half the area for up and down, respectively

            energy_convection_down =  air_flux.* ground.CONST.cp .*ground.TEMP.density_air .* ground.STATVAR.T(1:end-1,1);
            energy_convection_up = air_flux.* ground.CONST.cp .*ground.TEMP.density_air .* ground.STATVAR.T(2:end,1);
            
            ground.TEMP.d_energy(1:end-1,1) = ground.TEMP.d_energy(1:end-1,1) + energy_convection_up - energy_convection_down;
            ground.TEMP.d_energy(2:end,1) = ground.TEMP.d_energy(2:end,1) - energy_convection_up + energy_convection_down;

            %water
            density_air = tile.FORCING.TEMP.p./(287.058.*(ground.STATVAR.T+273.15)); 
            q_vol = 0.622 .* 6.112.* 100.* exp(22.46.*ground.STATVAR.T./(272.61+ground.STATVAR.T)) ./ tile.FORCING.TEMP.p .* density_air; %[kg water/ m3 air]
            ice_equivalent = q_vol ./ ground.CONST.rho_i; %m3 ice equivalent / m3 air

            energy_convection_latent_down =  air_flux.* q_vol(1:end-1,1) .* ground.CONST.L_s; %m3/sec .* kg/m3 .* J/kg = J/sec
            energy_convection_latent_up = air_flux.* q_vol(2:end,1) .* ground.CONST.L_s; %only latent part
            
            ground.TEMP.d_energy(1:end-1,1) = ground.TEMP.d_energy(1:end-1,1) + energy_convection_latent_up - energy_convection_latent_down;
            ground.TEMP.d_energy(2:end,1) = ground.TEMP.d_energy(2:end,1) - energy_convection_latent_up + energy_convection_latent_down;
            
            water_convection_down =  air_flux.* ice_equivalent(1:end-1,1);  %m3/sec
            water_convection_up = air_flux.* ice_equivalent(2:end,1);

            ground.TEMP.d_water_vapour(1:end-1,1) = ground.TEMP.d_water_vapour(1:end-1,1) + water_convection_up - water_convection_down;
            ground.TEMP.d_water_vapour(2:end,1) = ground.TEMP.d_water_vapour(2:end,1) - water_convection_up + water_convection_down;

            %vapour diffusion
            diffusivity_vapour_in_air = 1.7e-5 .* (tile.FORCING.TEMP.p ./ 1e5) .* ((ground.STATVAR.T+273.15)./293).^2.072 .* ground.STATVAR.porosity; %based on https://www.sciencedirect.com/science/article/pii/S1352231097003919
            
            diffusivity_vapour_in_air_eff = diffusivity_vapour_in_air(1:end-1,1) .* diffusivity_vapour_in_air(2:end,1) ./ ...
                (diffusivity_vapour_in_air(1:end-1,1) .* ground.STATVAR.layerThick(2:end,1) ./ 2 + diffusivity_vapour_in_air(2:end,1) .* ground.STATVAR.layerThick(1:end-1,1) ./ 2);
            diffusivity_vapour_in_air_eff(isnan(diffusivity_vapour_in_air_eff)) = 0;

            diffusion_flux = -diffusivity_vapour_in_air_eff .* (q_vol(1:end-1,1) - q_vol(2:end,1)) .* (ground.STATVAR.area(1:end-1,1)+ground.STATVAR.area(2:end,1))./2; 

            ground.TEMP.d_water_vapour(1:end-1,1) = ground.TEMP.d_water_vapour(1:end-1,1) + diffusion_flux./ ground.CONST.rho_i;
            ground.TEMP.d_water_vapour(2:end,1) = ground.TEMP.d_water_vapour(2:end,1) - diffusion_flux./ ground.CONST.rho_i;

            ground.TEMP.d_energy(1:end-1,1) = ground.TEMP.d_energy(1:end-1,1) + diffusion_flux.* ground.CONST.L_s;
            ground.TEMP.d_energy(2:end,1) = ground.TEMP.d_energy(2:end,1) - diffusion_flux.* ground.CONST.L_s;
        end
        
    end
end

