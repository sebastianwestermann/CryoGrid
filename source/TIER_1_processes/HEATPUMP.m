%========================================================================
% CryoGrid TIER1 library class for functions related to heat pumps
% NOTE: this class also contains code related to the free water freeze curve, 
% as well as functions computing thermal conductivity  
% S. Westermann, August 2025
%========================================================================

classdef HEATPUMP < BASE
    
    methods
        
        %-----derivatives----------
        %conductive heat flux between grid cells
        function ground = get_derivative_energy_heatPump(ground)
            %vertical fluxes
            fluxes = (ground.STATVAR.T(1:end-1,:) - ground.STATVAR.T(2:end,:)) .* ground.STATVAR.thermCond(1:end-1,:) .* ground.STATVAR.thermCond(2:end,:) ./...
                (ground.STATVAR.thermCond(1:end-1,:).* ground.STATVAR.layerThick(2:end,:)./2 +  ground.STATVAR.thermCond(2:end,:).* ground.STATVAR.layerThick(1:end-1,:)./2 );
            fluxes = fluxes .*(ground.STATVAR.area(1:end-1,:) + ground.STATVAR.area(2:end,:))./2; %[J/sec]
            
            d_energy=ground.STATVAR.energy.*0;
            d_energy(1,:) =  - fluxes(1,:);
            d_energy(2:end-1,:) = fluxes(1:end-1,:) - fluxes(2:end,:);
            d_energy(end,:) =  fluxes(end,:);
            
            ground.TEMP.d_energy = ground.TEMP.d_energy + d_energy;

            %lateral fluxes
            fluxes = (ground.STATVAR.T(:,1:end-1) - ground.STATVAR.T(:,2:end)) .* ground.STATVAR.thermCond(:,1:end-1) .* ground.STATVAR.thermCond(:,2:end) ./...
                (ground.STATVAR.thermCond(:,1:end-1).* ground.STATVAR.lateral_thickness(:,2:end)./2 +  ground.STATVAR.thermCond(:,2:end).* ground.STATVAR.lateral_thickness(:,1:end-1)./2 );
            fluxes = fluxes .* ground.STATVAR.layerThick(:,1:end-1) .* ground.STATVAR.contact_length; %[J/sec]
            
            d_energy=ground.STATVAR.energy.*0;
            d_energy(:,1) =  - fluxes(:,1);
            d_energy(:,2:end-1) = fluxes(:,1:end-1) - fluxes(:,2:end);
            d_energy(:,end) =  fluxes(:,end);
            
            ground.TEMP.d_energy = ground.TEMP.d_energy + d_energy;
        end
        
        
        %----diagnostic functions---------
        %free water freeze curve
        function ground = get_T_water_freeW_heatPump(ground)
            
            Lf = ground.CONST.L_f;
            c_w = ground.CONST.c_w;
            c_i = ground.CONST.c_i;
            c_o = ground.CONST.c_o;
            c_m = ground.CONST.c_m;
            c_cf = ground.PARA.heatCapacity_fluid;
            
            E_frozen = -Lf.*ground.STATVAR.waterIce;
            
            ground.STATVAR.T = double(ground.STATVAR.energy < E_frozen) .* (ground.STATVAR.energy - E_frozen) ./ (c_cf.*ground.STATVAR.cooling_fluid + c_i.*ground.STATVAR.waterIce + c_m.*ground.STATVAR.mineral + c_o.*ground.STATVAR.organic) + ...
                double(ground.STATVAR.energy >0) .* ground.STATVAR.energy ./ (c_cf.*ground.STATVAR.cooling_fluid + c_w.*ground.STATVAR.waterIce + c_m.*ground.STATVAR.mineral + c_o.*ground.STATVAR.organic);% T = 0 if energy between E_frozen and 0
            
            ground.STATVAR.ice = double(ground.STATVAR.energy <= E_frozen) .*ground.STATVAR.waterIce + double(ground.STATVAR.energy > E_frozen & ground.STATVAR.energy < 0) .* ground.STATVAR.energy ./ (-Lf);
            ground.STATVAR.water = double(ground.STATVAR.energy >= 0) .*ground.STATVAR.waterIce + double(ground.STATVAR.energy > - Lf.*ground.STATVAR.waterIce & ground.STATVAR.energy < 0) .* (ground.STATVAR.energy + Lf.*ground.STATVAR.waterIce) ./ Lf;
        end
        
        %calculate energy from temperature and water contents, free water
        %freeze curve
        function ground = get_E_freeW_heatPump(ground) %required for initialization

            T = ground.STATVAR.T;
            mineral= ground.STATVAR.mineral ./ (ground.STATVAR.layerThick .* ground.STATVAR.area);
            organic = ground.STATVAR.organic ./ (ground.STATVAR.layerThick .* ground.STATVAR.area);
            waterIce = ground.STATVAR.waterIce ./ (ground.STATVAR.layerThick .* ground.STATVAR.area);
            cooling_fluid = ground.STATVAR.cooling_fluid ./ (ground.STATVAR.layerThick .* ground.STATVAR.area);
            layerThick = ground.STATVAR.layerThick;
            area = ground.STATVAR.area;

            energy = T.*(mineral .* ground.CONST.c_m + organic .* ground.CONST.c_o + + cooling_fluid .* ground.PARA.heatCapacity_fluid + double(T>=0).*(waterIce .* ground.CONST.c_w) + ...
                double(T<0).*(waterIce .* ground.CONST.c_i )) - double(T<0) .* (waterIce) .* ground.CONST.L_f;
            
            ground.STATVAR.waterIce = waterIce .* layerThick .* area; % [m3]
            ground.STATVAR.mineral = mineral .* layerThick .* area; % [m3]
            ground.STATVAR.organic = organic .* layerThick .* area; % [m3]
            ground.STATVAR.cooling_fluid = cooling_fluid .* layerThick .* area;
            ground.STATVAR.energy = energy .* layerThick .* area;  % [J]
            
            ground.STATVAR.water = double(T>=0) .* waterIce .* layerThick .* area;  % [m3]
            ground.STATVAR.ice = double(T<0) .* waterIce .* layerThick .* area; %[m3]
            ground = conductivity(ground);
            
        end
        
        %---thermal conductivities--------------
        function ground = conductivity_mixing_squares_heatPump(ground)
            
            water = ground.STATVAR.water./ground.STATVAR.layerThick ./ ground.STATVAR.area;
            ice = ground.STATVAR.ice ./ ground.STATVAR.layerThick./ ground.STATVAR.area;
            mineral = ground.STATVAR.mineral./ground.STATVAR.layerThick./ ground.STATVAR.area;
            organic = ground.STATVAR.organic./ground.STATVAR.layerThick./ ground.STATVAR.area;
            cooling_fluid = ground.STATVAR.cooling_fluid ./ (ground.STATVAR.layerThick .* ground.STATVAR.area);
            air = 1 - mineral - organic - water - ice;
            
            ground.STATVAR.thermCond = ((cooling_fluid+water).* ground.CONST.k_w.^0.5 + ice.* ground.CONST.k_i.^0.5 ...
                + mineral.* ground.CONST.k_m.^0.5 + organic.* ground.CONST.k_o.^0.5 + air.* ground.CONST.k_a.^0.5).^2;
        end

    end
end

