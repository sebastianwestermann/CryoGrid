%========================================================================
% CryoGrid INTERACTION (IA) class handling interactions between
% VEGETATION_CLM5_seb and
% GROUND (RichardsEqW, Xice)
% with CHILD SNOW_crocus2_bucketW_seb
% R. B. Zwegiel, February 2022
%========================================================================

classdef IA_SEB_vegetation_CLM5_GROUND_Xice_snow < IA_SEB_vegetation_CLM5
   
    methods
        
        function ia_seb_water = get_boundary_condition_m(ia_seb_water, tile)
            % below-canopy equivalent to get_boundary_condition_u in GROUND_freezeC_RichardsEqW_Xice_seb_snow
            ground = ia_seb_water.NEXT; % for easy access
            
            if ground.CHILD == 0 % No CHILD exists
                ia_seb_water = get_boundary_condition_m_GROUND_Xice(ia_seb_water, tile);
                
                if tile.FORCING.TEMP.snowfall > 0 % Create CHILD
                    ground.CHILD = copy(tile.STORE.SNOW);
                    ground.CHILD.PARENT = ground;
                    ground.CHILD.NEXT = ground; 
                    ground.IA_CHILD = get_IA_class(class(ground.CHILD), class(ground));
                    ground.IA_CHILD.PREVIOUS = ground.CHILD;
                    ground.IA_CHILD.NEXT = ground; %SNOW CHILD created
                    
                    ground.CHILD = get_boundary_condition_m_create_CHILD(ia_seb_water, tile);  %initialize with fresh snowfall
                end
            else % CHILD exists
                total.area = ground.STATVAR.area;
                total.waterIce  = ground.STATVAR.waterIce;
                total.water = ground.STATVAR.water;
                total.ice = ground.STATVAR.ice;
                total.mineral = ground.STATVAR.mineral;
                total.organic = ground.STATVAR.organic;
                total_XwaterIce = ground.STATVAR.XwaterIce;
                total_Xice = ground.STATVAR.Xice;
                total_Xwater = ground.STATVAR.Xwater;
                
                ground.STATVAR.area = ground.STATVAR.area - ground.CHILD.STATVAR.area(1); % use snow free area from now on
                reduction = ground.STATVAR.area(1) ./ total.area(1);
                ground.STATVAR.waterIce = ground.STATVAR.waterIce .* reduction;
                ground.STATVAR.water = ground.STATVAR.water .* reduction;
                ground.STATVAR.ice = ground.STATVAR.ice .* reduction;
                ground.STATVAR.mineral = ground.STATVAR.mineral .* reduction;
                ground.STATVAR.organic = ground.STATVAR.organic .* reduction;
                ground.STATVAR.XwaterIce = ground.STATVAR.XwaterIce .* reduction;
                ground.STATVAR.Xice = ground.STATVAR.Xice .* reduction;
                ground.STATVAR.Xwater = ground.STATVAR.Xwater .* reduction;
                
                %-------------
                
                ground.CHILD.STATVAR.Lstar = ground.STATVAR.Lstar;
                
                ia_seb_water = get_boundary_condition_m_CHILD(ia_seb_water, tile); % snowfall, throughfall, Qh, Qe & d_energy
                
                ia_seb_water = get_boundary_condition_m_GROUND_Xice(ia_seb_water, tile); % throughfall, Qh, Qe & d_energy
                
                get_IA_CHILD_boundary_condition_u(ground.IA_CHILD, tile); %call designated mandatory function for CHILD-PARENT interactions in the IA class governing IA between SNOW and GROUND
                
                % Mix the SEB fluxes from snow and ground
                ground.STATVAR.Qh = (ground.STATVAR.area(1,1) .* ground.STATVAR.Qh + ground.CHILD.STATVAR.area .* ground.CHILD.STATVAR.Qh) ./ total.area(1,1);
                ground.STATVAR.Qe = (ground.STATVAR.area(1,1) .* ground.STATVAR.Qe + ground.CHILD.STATVAR.area .* ground.CHILD.STATVAR.Qe) ./ total.area(1,1);
                
                %reassign the true totals of ground
                ground.STATVAR.area = total.area;
                ground.STATVAR.waterIce = total.waterIce;
                ground.STATVAR.water = total.water;
                ground.STATVAR.ice = total.ice;
                ground.STATVAR.mineral = total.mineral;
                ground.STATVAR.organic = total.organic;
                ground.STATVAR.XwaterIce = total_XwaterIce;
                ground.STATVAR.Xice = total_Xice;
                ground.STATVAR.Xwater = total_Xwater;
            end
        
        end   
        
        function ia_seb_water = canopy_drip(ia_seb_water, tile)
            if ia_seb_water.NEXT.CHILD ~= 0
                ia_seb_water = canopy_drip_Xice_SNOW_CHILD(ia_seb_water, tile);
            else
                ia_seb_water = canopy_drip_Xice(ia_seb_water, tile);
            end
        end

        function r_soil = ground_resistance_evap(ia_seb_water, tile)
            ground = ia_seb_water.NEXT;
            vol_water_first_cell = (ground.STATVAR.waterIce(1) + ground.STATVAR.XwaterIce(1)) ./ (ground.STATVAR.layerThick(1,1) .* ground.STATVAR.area(1,1));
            reduce_yes_no = vol_water_first_cell < ground.STATVAR.field_capacity(1,1);
            betaCLM4_5 = 1 +  double(reduce_yes_no) .* (-1 +  0.25 .* (1-(cos(pi() .* vol_water_first_cell ./ ground.STATVAR.field_capacity(1,1)))).^2);

            r_soil = min(1e10, 250.*((1./betaCLM4_5).^0.5 -1));
            if ground.CHILD ~= 0
                if ground.CHILD.STATVAR.area/ground.STATVAR.area(1) > 0.5
                    % Simple swich for now - use r_soil for ground until child covers 50%
                    r_soil = 0;
                end
            end
        end
    
    end 
    
end