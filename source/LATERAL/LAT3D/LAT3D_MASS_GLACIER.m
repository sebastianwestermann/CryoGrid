%========================================================================
% CryoGrid LATERAL_IA class LAT3D_MASS_GLACIER 
% simulates lateral mass flow between pairs of CryoGrid stratigraphies
% for the topmost unconfined aquifer
% S. Westermann, Apr 2024
%========================================================================


classdef LAT3D_MASS_GLACIER < BASE_LATERAL

    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        
        function lateral = provide_CONST(lateral)
            lateral.CONST.day_sec = []; %24 .* 3600;
            lateral.CONST.rho_i = [];
            lateral.CONST.g = [];
        end
        
        function lateral = provide_PARA(lateral)
            lateral.PARA.ia_time_increment = []; %0.25; %must be a multiple of the time increment of the main lateral class
            %lateral.PARA.ia_time_next = [];
        end
        
        function lateral = provide_STATVAR(lateral)
%             lateral.STATVAR.subsurface_run_off = [];
%             lateral.STATVAR.surface_runoff = [];
        end

        function lateral = finalize_init(lateral, tile)
%             lateral.STATVAR.subsurface_run_off = 0;
%             lateral.STATVAR.surface_runoff = 0;
        end
        
        %--- time integration---------------
        
        function lateral = pull(lateral, tile)
            
            %to pull: all the variables that are also used in the merging
%                         extensive_variables = {'layerThick'; 'layerThick_glacier'; 'layerThick_sediment'; 'waterIce'; 'XwaterIce'; 'Xwater'; 'mineral'; 'organic'; 'energy'};
%             intensive_variables = {'area'; 'ice_fraction'; 'field_capacity_sediment'; 'satHydraulicConductivity_sediment'}; 
            %top and bottom altitude for glacier, bottom defined as
            %ice-rich layer, 'layerThick_glacier' > 0 
            
            %big difference in surface elevation -> large flow speed
            %ice thickness -> thicker ice moves faster
            %average ice T - Glens flow law
            
            %make a unified layerThick data set that is then used to get
            %cell overlap: mean glacier thickness = (g1+g2)/2. s1 =
            %(g1+g2)/(2 g1); sf2 = (g1+g2)/(2 g2)

                   
            CURRENT = lateral.PARENT.TOP.NEXT;
            
            while ~(strcmp(class(CURRENT), 'Bottom'))
                CURRENT = lateral3D_pull_mass_GLACIER(CURRENT, lateral);
                CURRENT = CURRENT.NEXT;
            end
        end
        
        function lateral = get_derivatives(lateral, tile) %no need to loop through stratigraphy, all the information is in lateral.PARENT

            lateral.PARENT = get_overlap_cells_mass(lateral.PARENT);
            
            lateral.PARENT.STATVAR.d_glacier = lateral.PARENT.STATVAR.glacier_fraction .* 0;
            lateral.PARENT.STATVAR.d_sediment = lateral.PARENT.STATVAR.glacier_fraction .* 0;
            lateral.PARENT.STATVAR.d_Xwater = lateral.PARENT.STATVAR.glacier_fraction .* 0;
            lateral.PARENT.STATVAR.d_waterIce = lateral.PARENT.STATVAR.glacier_fraction .* 0;
            lateral.PARENT.STATVAR.d_ice = lateral.PARENT.STATVAR.glacier_fraction .* 0;
            lateral.PARENT.STATVAR.d_mineral = lateral.PARENT.STATVAR.glacier_fraction .* 0;
            lateral.PARENT.STATVAR.d_organic = lateral.PARENT.STATVAR.glacier_fraction .* 0;
            lateral.PARENT.STATVAR.d_energy = lateral.PARENT.STATVAR.glacier_fraction .* 0;
            lateral.PARENT.STATVAR.d_glacier_ice_fraction = lateral.PARENT.STATVAR.glacier_fraction .* 0;
            lateral.PARENT.STATVAR.d_field_capacity_sediment = lateral.PARENT.STATVAR.glacier_fraction .* 0;
            lateral.PARENT.STATVAR.d_satHydraulicConductivity_sediment = lateral.PARENT.STATVAR.glacier_fraction .* 0;
            
            
            for j=1:size(lateral.PARENT.ENSEMBLE,1)
                if lateral.PARENT.PARA.connected(lateral.PARENT.STATVAR.index, lateral.PARENT.ENSEMBLE{j,1}.index)  %add if water is available at all
                    %flow between saturated cells
                    contact_length = lateral.PARENT.PARA.contact_length(lateral.PARENT.STATVAR.index, lateral.PARENT.ENSEMBLE{j,1}.index);
                    distance = lateral.PARENT.PARA.distance(lateral.PARENT.STATVAR.index, lateral.PARENT.ENSEMBLE{j,1}.index);
                    
                    ice_thickness = (sum(lateral.PARENT.STATVAR.layerThick) + sum(lateral.PARENT.ENSEMBLE{j,1}.layerThick))./2;
                    average_T = (sum(lateral.PARENT.STATVAR.layerThick .* lateral.PARENT.STATVAR.T) + sum(lateral.PARENT.ENSEMBLE{j,1}.layerThick + lateral.PARENT.ENSEMBLE{j,1}.T)) ./ (sum(lateral.PARENT.STATVAR.layerThick) + sum(lateral.PARENT.ENSEMBLE{j,1}.layerThick));

                    %Glens flow law for shallow ice approximation
                    glensFlowParameter = 4.1e-24 .* exp(0.23 .* min(0, average_T)); %[Pa-3 sec-1]
                    glacier_flow_velocity = - 0.4 .* glensFlowParameter .* (lateral.CONST.rho_i .* lateral.CONST.g).^3 .*  ice_thickness.^4 .* ((lateral.PARENT.STATVAR.upperPos - lateral.PARENT.ENSEMBLE{j,1}.upperPos) ./ distance).^3;
                    lateral.PARENT.ENSEMBLE{j,1}.glacier_flow_velocity = glacier_flow_velocity;

                    %                   %in [m/sec], must be multipled with contact length and vertical overlap

                    mf = contact_length .* lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec;
                    if glacier_flow_velocity <= 0 %own tile looses
                        loosing_tile = lateral.PARENT.STATVAR;
                        ind = 1; 
                    else %own tile gains from neighboring tile
                        loosing_tile = lateral.PARENT.ENSEMBLE{j,1};
                        ind = 2;
                    end
                    
                    for i=1:size(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass,1) %does not do anything when overlap is empty!
                        
                        lateral.PARENT.STATVAR.d_glacier(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,1),1) = lateral.PARENT.STATVAR.d_glacier(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,1),1) ...
                            + mf .* lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,3) .* glacier_flow_velocity .* loosing_tile.glacier_fraction(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(ind,1),1);
                        lateral.PARENT.STATVAR.d_sediment(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,1),1) = lateral.PARENT.STATVAR.d_sediment(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,1),1) ...
                            + mf .* lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,3) .* glacier_flow_velocity .* loosing_tile.sediment_fraction(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(ind,1),1);
                        lateral.PARENT.STATVAR.d_Xwater(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,1),1) = lateral.PARENT.STATVAR.d_Xwater(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,1),1) ...
                            + mf .* lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,3) .* glacier_flow_velocity .* loosing_tile.Xwater_fraction(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(ind,1),1);
                        lateral.PARENT.STATVAR.d_waterIce(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,1),1) = lateral.PARENT.STATVAR.d_waterIce(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,1),1) ...
                            + mf .* lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,3) .* glacier_flow_velocity .* loosing_tile.waterIce_fraction(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(ind,1),1);
                        lateral.PARENT.STATVAR.d_ice(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,1),1) = lateral.PARENT.STATVAR.d_ice(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,1),1) ...
                            + mf .* lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,3) .* glacier_flow_velocity .* loosing_tile.ice_fraction(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(ind,1),1);
                        lateral.PARENT.STATVAR.d_mineral(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,1),1) = lateral.PARENT.STATVAR.d_mineral(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,1),1) ...
                            + mf .* lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,3) .* glacier_flow_velocity .* loosing_tile.mineral_fraction(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(ind,1),1);
                        lateral.PARENT.STATVAR.d_organic(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,1),1) = lateral.PARENT.STATVAR.d_organic(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,1),1) ...
                            + mf .* lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,3) .* glacier_flow_velocity .* loosing_tile.organic_fraction(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(ind,1),1);
                        lateral.PARENT.STATVAR.d_energy(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,1),1) = lateral.PARENT.STATVAR.d_energy(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,1),1) ...
                            + mf .* lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,3) .* glacier_flow_velocity .* loosing_tile.energy(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(ind,1),1);
                        
                        lateral.PARENT.STATVAR.d_glacier_ice_fraction(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,1),1) = lateral.PARENT.STATVAR.d_glacier_ice_fraction(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,1),1) ...
                            + mf .* lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,3) .* glacier_flow_velocity .* loosing_tile.glacier_ice_fraction(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(ind,1),1) .* loosing_tile.glacier_fraction(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(ind,1),1);
                       lateral.PARENT.STATVAR.d_field_capacity_sediment(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,1),1) = lateral.PARENT.STATVAR.d_field_capacity_sediment(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,1),1) ...
                            + mf .* lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,3) .* glacier_flow_velocity .* loosing_tile.field_capacity_sediment(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(ind,1),1) .* loosing_tile.sediment_fraction(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(ind,1),1);
                        lateral.PARENT.STATVAR.d_satHydraulicConductivity_sediment(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,1),1) = lateral.PARENT.STATVAR.d_satHydraulicConductivity_sediment(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,1),1) ...
                            + mf .* lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,3) .* glacier_flow_velocity .* loosing_tile.satHydraulicConductivity_sediment(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(ind,1),1) .* loosing_tile.sediment_fraction(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(ind,1),1);
%                       lateral.PARENT.STATVAR.d_glacier_ice_fraction_scaled(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,1),1) = lateral.PARENT.STATVAR.d_glacier_ice_fraction_scaled(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,1),1) ...
%                             + mf .* lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,3) .* glacier_flow_velocity .* loosing_tile.glacier_ice_fraction_scaled(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(ind,1),1);
%                        lateral.PARENT.STATVAR.d_field_capacity_sediment_scaled(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,1),1) = lateral.PARENT.STATVAR.d_field_capacity_sediment_scaled(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,1),1) ...
%                             + mf .* lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,3) .* glacier_flow_velocity .* loosing_tile.field_capacity_sediment_scaled(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(ind,1),1);
%                         lateral.PARENT.STATVAR.d_satHydraulicConductivity_sediment_scaled(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,1),1) = lateral.PARENT.STATVAR.d_satHydraulicConductivity_sediment_scaled(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,1),1) ...
%                             + mf .* lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(i,3) .* glacier_flow_velocity .* loosing_tile.satHydraulicConductivity_sediment_scaled(lateral.PARENT.ENSEMBLE{j,1}.overlap_mass(ind,1),1);

                    end

                end
            end
        end

        function lateral = push(lateral, tile)

                CURRENT = lateral.PARENT.TOP.NEXT; %find correct stratigraphy class
                while ~(strcmp(class(CURRENT), 'Bottom'))
                    CURRENT = lateral3D_push_mass_GLACIER(CURRENT, lateral);
                    CURRENT = CURRENT.NEXT;
                end
        end
       
        
        function lateral = set_ACTIVE(lateral, i, t)
            lateral.PARENT.ACTIVE(i,1) = 0;
            if t + lateral.PARENT.IA_TIME_INCREMENT >= lateral.PARA.ia_time_next - 1e-7
                lateral.PARENT.ACTIVE(i,1) = 1;
                %lateral.PARA.ia_time_next = t + lateral.PARENT.IA_TIME_INCREMENT + lateral.PARA.ia_time_increment;
                lateral.PARA.ia_time_next = lateral.PARA.ia_time_next + lateral.PARA.ia_time_increment;
                %disp(lateral.PARA.ia_time_next-floor(lateral.PARA.ia_time_next));
            end
        end
        
        function lateral = set_ia_time(lateral, t)
            lateral.PARA.ia_time_next = t;
        end
        
        
        %-------------param file generation-----
%         function ground = param_file_info(ground)
%             ground = param_file_info@BASE_LATERAL(ground);
%             
%             ground.PARA.class_category = 'LATERAL_IA';
%             
%             ground.PARA.options = [];
%             ground.PARA.STATVAR = [];
%             
%             ground.PARA.default_value.hardBottom_cutoff = {0.03};
%             ground.PARA.comment.hardBottom_cutoff = {'hard bottom  = no water flow if saturated and water content below [vol. water content, -]'};
%             
%             ground.PARA.default_value.GaMa_coefficient = {15};
%             ground.PARA.comment.GaMa_coefficient = {'Gauckler-Manning coefficient, https://en.wikipedia.org/wiki/Manning_formula'};
%             
%             ground.PARA.default_value.tortuosity = {1};
%             ground.PARA.comment.tortuosity = {'multiply direct distance with this factor to get flow path length'};
%             
%             ground.PARA.default_value.ia_time_increment = {0.25};
%             ground.PARA.comment.ia_time_increment ={'time step [days], must be multiple of LATERAL class timestep'};
%         end
    end
    
end


