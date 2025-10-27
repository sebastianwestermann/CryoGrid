%========================================================================
% CryoGrid LATERAL_IA class LAT3D_WATER_UNCONFINED_AQUIFER 
% simulates lateral water flow between pairs of CryoGrid stratigraphies
% for the topmost unconfined aquifer.
% NOTE: no flow if there is no unconfined aquifer for one of the two stratigraphies,
% e.g. if the first cell is saturated with ice. Use LAT_3D_WATER instead. 
% S. Westermann, Oct 2020
%========================================================================


% Required variables:
% Depths: absolute elevation of each grid cell
% Gravitational_potential: same as depth for unsaturated, hydrostatic potential for saturated, based on water_table_elevation, or the first saturated cell in case of Richards Eq class
% When getting to unsaturated cell, go back to depth and start new (assumes that water can always drain, also below a hard bottom)
% Water_potential: zero when saturated, and for classes w.o. water potential in hydrology
% Make two variables out of water-status
% Water_available: water available for drainage in a cell
% Water_deficit: max water that can be taken up in a cell 
% Hydraulic conductivity: always has a value, set to zero when frozen/hard_bottom
% Ground_surface_elevation: as before;
% 
% Potential = Gravitational_potential – depths + Water_potential.
% Flux = -hyd_cond * (potential1+ ground_surface_elevation1 - potential2-ground_surface_elevation2) / distance .* cross_section;
% Then reduce with the two water variables
% 
% To push: take out/add water flux



classdef LAT3D_WATER_UNCONFINED_AQUIFER_RICHARDS_EQ_OVERLAND_FLOW < BASE_LATERAL

    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        
        
        function lateral = provide_CONST(lateral)
            lateral.CONST.day_sec = []; %24 .* 3600;
        end
        
        function lateral = provide_PARA(lateral)
            lateral.PARA.hardBottom_cutoff = []; %0.03; %hard bottom if saturated and water content below
            lateral.PARA.GaMa_coefficient = [];
            lateral.PARA.ia_time_increment = []; %0.25; %must be a multiple of the time increment of the main lateral class
            lateral.PARA.tortuosity = 1;
        end
        
        function lateral = provide_STATVAR(lateral)
            lateral.STATVAR.subsurface_run_off = [];
            lateral.STATVAR.surface_runoff = [];
        end

        function lateral = finalize_init(lateral, tile)
            lateral.STATVAR.subsurface_run_off = 0;
            lateral.STATVAR.surface_runoff = 0;
        end
        
        %--- time integration---------------
        
        function lateral = pull(lateral, tile)
            %calculate as in ground = lateral_push_water_reservoir_RichardsEq_simple(ground, lateral)
            lateral.PARENT.STATVAR.depths = []; %not sure if needed?
            lateral.PARENT.STATVAR.water_status = [];
            lateral.PARENT.STATVAR.hydraulicConductivity = [];
            lateral.PARENT.STATVAR.matric_potential_head = [];
            lateral.PARENT.STATVAR.hydrostatic_head = [];
            lateral.PARENT.STATVAR.T_water = [];
            lateral.PARENT.STATVAR.water_depth = 0; %depth of the free water at the surface flowing according to Gauckler-Manning 
        %    lateral.PARENT.STATVAR.ground_surface_elevation = []; %NEW SEBASTIAN for hillslope flow between grid cells
            
            lateral.TEMP.open_system = 1; %start with open system
            
            CURRENT = lateral.PARENT.TOP.NEXT;
            CURRENT = lateral3D_pull_water_overland_flow(CURRENT, lateral);

            CURRENT = lateral.PARENT.TOP.NEXT;
            while ~(strcmp(class(CURRENT), 'Bottom')) && lateral.TEMP.open_system == 1
                CURRENT = lateral3D_pull_water_unconfined_aquifer_RichardsEq(CURRENT, lateral);
                CURRENT = CURRENT.NEXT;
            end
        end
        
        function lateral = get_derivatives(lateral, tile) %no need to loop through stratigraphy, all the information is in lateral.PARENT

            lateral.PARENT = get_overlap_cells2(lateral.PARENT, 'depths', 'overlap');
            
            
            %calculate fluxes
            flux = lateral.PARENT.STATVAR.hydraulicConductivity .* 0;
            flux_energy = flux;
            for j=1:size(lateral.PARENT.ENSEMBLE,1)
                if lateral.PARENT.PARA.connected(lateral.PARENT.STATVAR.index, lateral.PARENT.ENSEMBLE{j,1}.index)  %add if water is available at all
                    %flow between saturated cells
                    
                    %calculate offset from grid cell midpoint of overlap
                    %section midpoint
                    offset_from_midpoint = [];
                    overlap = lateral.PARENT.ENSEMBLE{j,1}.overlap;
                    depth1=lateral.PARENT.STATVAR.depths;
                    depth2 = lateral.PARENT.ENSEMBLE{j,1}.depths;
                    for k=1:size(overlap,1)
                        offset_from_midpoint=[offset_from_midpoint; ...
                            [-(depth1(overlap(k,1),1)+depth1(overlap(k,1)+1,1))/2+(min(depth1(overlap(k,1),1), depth2(overlap(k,2),1))-overlap(k,3)./2)...
                            -(depth2(overlap(k,2),1)+depth2(overlap(k,2)+1,1))/2+(min(depth1(overlap(k,1),1), depth2(overlap(k,2),1))-overlap(k,3)./2)]];
                    end
                    

                    contact_length = lateral.PARENT.PARA.contact_length(lateral.PARENT.STATVAR.index, lateral.PARENT.ENSEMBLE{j,1}.index);
                    distance = lateral.PARENT.PARA.distance(lateral.PARENT.STATVAR.index, lateral.PARENT.ENSEMBLE{j,1}.index);
                    for i=1:size(lateral.PARENT.ENSEMBLE{j,1}.overlap,1) %does not do anything when overlap is empty!
                        
                        hc1 = lateral.PARENT.STATVAR.hydraulicConductivity(lateral.PARENT.ENSEMBLE{j,1}.overlap(i,1),1);
                        hc2 = lateral.PARENT.ENSEMBLE{j,1}.hydraulicConductivity(lateral.PARENT.ENSEMBLE{j,1}.overlap(i,2),1);
                        hydraulicConductivity = hc1 .* hc2 ./ (hc1 .* distance./2 + hc2 .* distance./2); %change to different distances later
                        hydraulicConductivity(isnan(hydraulicConductivity)) = 0;

                        flux_i = contact_length .* lateral.PARENT.ENSEMBLE{j,1}.overlap(i,3) .* hydraulicConductivity .* ...
                            -(lateral.PARENT.STATVAR.matric_potential_head(lateral.PARENT.ENSEMBLE{j,1}.overlap(i,1),1) - lateral.PARENT.ENSEMBLE{j,1}.matric_potential_head(lateral.PARENT.ENSEMBLE{j,1}.overlap(i,2),1) + ...
                            double(lateral.PARENT.STATVAR.hydrostatic_head(lateral.PARENT.ENSEMBLE{j,1}.overlap(i,1),1)~=0).*(lateral.PARENT.STATVAR.hydrostatic_head(lateral.PARENT.ENSEMBLE{j,1}.overlap(i,1),1) - offset_from_midpoint(i,1)) - ...
                            double(lateral.PARENT.ENSEMBLE{j,1}.hydrostatic_head(lateral.PARENT.ENSEMBLE{j,1}.overlap(i,2),1)~=0).*(lateral.PARENT.ENSEMBLE{j,1}.hydrostatic_head(lateral.PARENT.ENSEMBLE{j,1}.overlap(i,2),1)  - offset_from_midpoint(i,2)) + ...
                            double(lateral.PARENT.PARA.hill_slope) .* (lateral.PARENT.STATVAR.ground_surface_elevation - lateral.PARENT.ENSEMBLE{j,1}.ground_surface_elevation)); 
%           flux_i = contact_length .* lateral.PARENT.ENSEMBLE{j,1}.overlap(i,3) .* hydraulicConductivity .* ...
%                             -(lateral.PARENT.STATVAR.matric_potential_head(lateral.PARENT.ENSEMBLE{j,1}.overlap(i,1),1) - lateral.PARENT.ENSEMBLE{j,1}.matric_potential_head(lateral.PARENT.ENSEMBLE{j,1}.overlap(i,2),1) + ...
%                             (lateral.PARENT.STATVAR.hydrostatic_head(lateral.PARENT.ENSEMBLE{j,1}.overlap(i,1),1) - lateral.PARENT.ENSEMBLE{j,1}.hydrostatic_head(lateral.PARENT.ENSEMBLE{j,1}.overlap(i,2),1)) + ...
%                             double(lateral.PARENT.PARA.hill_slope) .* (lateral.PARENT.STATVAR.ground_surface_elevation - lateral.PARENT.ENSEMBLE{j,1}.ground_surface_elevation)); 
                            
                        flux(lateral.PARENT.ENSEMBLE{j,1}.overlap(i,1),1) = flux(lateral.PARENT.ENSEMBLE{j,1}.overlap(i,1),1) + flux_i;
                        flux_energy(lateral.PARENT.ENSEMBLE{j,1}.overlap(i,1),1) = flux_energy(lateral.PARENT.ENSEMBLE{j,1}.overlap(i,1),1) + flux_i .* lateral.PARENT.CONST.c_w .* ...
                            (double(flux_i<0).*lateral.PARENT.STATVAR.T_water(lateral.PARENT.ENSEMBLE{j,1}.overlap(i,1),1) + ...
                            double(flux_i>=0) .* lateral.PARENT.ENSEMBLE{j,1}.T_water(lateral.PARENT.ENSEMBLE{j,1}.overlap(i,2),1));
                        %flux <0 outflow of own worker
                        %flux >0 inflow into own worker
                    end
                    
                    
                    %overland flow according to Gauckler-Manning
                    if (lateral.PARENT.STATVAR.water_depth > 1e-4 && lateral.PARENT.ENSEMBLE{j,1}.water_depth > 1e-4  && ~isempty(lateral.PARENT.STATVAR.T_water) && ~isempty(lateral.PARENT.ENSEMBLE{j,1}.T_water)) ...
                            || (lateral.PARENT.STATVAR.water_depth <= 1e-4 && lateral.PARENT.ENSEMBLE{j,1}.water_depth > 1e-4  && ~isempty(lateral.PARENT.ENSEMBLE{j,1}.T_water)  && ~isempty(lateral.PARENT.STATVAR.T_water) && lateral.PARENT.STATVAR.depths(1,1) < lateral.PARENT.ENSEMBLE{j,1}.depths(1,1)) ...
                            || (lateral.PARENT.STATVAR.water_depth > 1e-4 && lateral.PARENT.ENSEMBLE{j,1}.water_depth <= 1e-4  && ~isempty(lateral.PARENT.ENSEMBLE{j,1}.T_water)  && ~isempty(lateral.PARENT.STATVAR.T_water) && lateral.PARENT.STATVAR.depths(1,1) > lateral.PARENT.ENSEMBLE{j,1}.depths(1,1))
                        if lateral.PARENT.STATVAR.depths(1,1) < lateral.PARENT.ENSEMBLE{j,1}.depths(1,1) && lateral.PARENT.STATVAR.water_depth < 1e-4
                            lateral.PARENT.STATVAR.water_depth = 1e-4;
                            lateral.PARENT.STATVAR.max_flow = 1e-2 .* lateral.PARENT.STATVAR.area_flow;
                        end
                        if lateral.PARENT.STATVAR.depths(1,1) > lateral.PARENT.ENSEMBLE{j,1}.depths(1,1) && lateral.PARENT.ENSEMBLE{j,1}.water_depth < 1e-4
                            lateral.PARENT.ENSEMBLE{j,1}.water_depth = 1e-4;
                            lateral.PARENT.ENSEMBLE{j,1}.max_flow = 1e-2 .* lateral.PARENT.ENSEMBLE{j,1}.area_flow;
                        end
                        
                        
                        %could be necessary to  limit flow for numerical
                        %stability?
                        gradient = -(lateral.PARENT.STATVAR.depths(1,1) - lateral.PARENT.ENSEMBLE{j,1}.depths(1,1)) ./ (distance.*lateral.PARA.tortuosity);
                        if gradient <= 0  %own realization higher
                            T_water =  lateral.PARENT.STATVAR.T_water(1,1);
                            depth1 = min(lateral.PARENT.STATVAR.depths(1,1) - (lateral.PARENT.ENSEMBLE{j,1}.depths(1,1) - lateral.PARENT.ENSEMBLE{j,1}.water_depth), lateral.PARENT.STATVAR.water_depth);
                            depth2 = lateral.PARENT.ENSEMBLE{j,1}.water_depth;
                        elseif gradient > 0  %ensemble realization higher
                            T_water = lateral.PARENT.ENSEMBLE{j,1}.T_water(1,1);
                            depth1 = min(lateral.PARENT.ENSEMBLE{j,1}.depths(1,1) - (lateral.PARENT.STATVAR.depths(1,1) - lateral.PARENT.STATVAR.water_depth), lateral.PARENT.ENSEMBLE{j,1}.water_depth);
                            depth2 = lateral.PARENT.STATVAR.water_depth;
                        end
                        
                        
                        %velocity = lateral.PARA.GaMa_coefficient .* real(min(lateral.PARENT.STATVAR.water_depth.^(2/3), lateral.PARENT.ENSEMBLE{j,1}.water_depth.^(2/3)) .* abs(gradient).^0.5);
                        velocity = lateral.PARA.GaMa_coefficient .* real(((depth1+depth2)./2).^(2/3) .* abs(gradient).^0.5);
                        %flow =  velocity .* min(lateral.PARENT.STATVAR.water_depth, lateral.PARENT.ENSEMBLE{j,1}.water_depth) .* contact_length;
                        flow =  velocity .* (depth1+depth2)./2 .* contact_length;
                        %                         if flow <0
                        %                             flow = max(flow, -lateral.PARENT.STATVAR.max_flow ./ (lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec));
                        %                         end
                        %flow that would lead to the same water level in
                        %both tiles
                        max_flow2same_level = lateral.PARENT.STATVAR.area_flow .* lateral.PARENT.ENSEMBLE{j,1}.area_flow .* (lateral.PARENT.STATVAR.depths(1,1) - lateral.PARENT.ENSEMBLE{j,1}.depths(1,1)) ./ (lateral.PARENT.STATVAR.area_flow + lateral.PARENT.ENSEMBLE{j,1}.area_flow) ./ (lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec);
                        
                        %TEST Sebastian
                        %flow = sign(gradient) .* min(abs(max_flow2same_level./8), min(min(flow, lateral.PARENT.STATVAR.max_flow ./ (lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec)), lateral.PARENT.ENSEMBLE{j,1}.max_flow ./ (lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec)));
                        
                        
                        if gradient < 0 %own realization higher, looses water
                            flow = sign(gradient) .* min(abs(max_flow2same_level./8), min(flow, lateral.PARENT.STATVAR.max_flow ./ (lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec)));
                        else
                            flow = sign(gradient) .* min(abs(max_flow2same_level./8), min(flow, lateral.PARENT.ENSEMBLE{j,1}.max_flow ./ (lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec)));
                        end
                        %end TEST Sebastian
                        
                        
                        flow_energy = flow .* lateral.PARENT.CONST.c_w .* T_water;

                        if ~isempty(flux)
                        flux(1,1) = flux(1,1) + flow;
                        flux_energy(1,1) = flux_energy(1,1) + flow_energy;
                        else
                        flux = flow;
                        flux_energy = flow_energy;
                        end
                        
                        lateral.STATVAR.flow = flow;
                        lateral.STATVAR.max_flow = -lateral.PARENT.STATVAR.max_flow ./ (lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec);
                        lateral.STATVAR.water_depth = lateral.PARENT.STATVAR.water_depth;
                        
                    end
                    
                end
            end
            
            
            
            %a->inflow of seepage face water
            %1. find the uppermost cell of the own worker that is saturated
            %water (water_status>0)
            %2. go down the cells of the other worker with water_status >0,
            %calculate seepage head and the flux
            %until either the bottom is reached or the cell's lower level
            %is below the water level of the own worker
            %4. add the fluxes up and add the flux to cell found in 1.
            
            %b -> outflow of seepgae face water
            %steps 2 and 3, but subrtact the water from each cell

            lateral.PARENT.STATVAR.water_flux = flux .* lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec;
            lateral.PARENT.STATVAR.water_flux_energy = flux_energy .* lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec;
            lateral.STATVAR.subsurface_run_off = sum(lateral.PARENT.STATVAR.water_flux);
            
        end

        function lateral = push(lateral, tile)
            if ~isempty(lateral.PARENT.STATVAR.water_flux)
                
                lateral.PARENT.STATVAR.water_up = 0;
                lateral.PARENT.STATVAR.water_up_energy = 0;
                
                CURRENT = lateral.PARENT.TOP.NEXT; %find correct stratigraphy class
                size_to_here = size(CURRENT.STATVAR.energy,1);
                while ~(strcmp(class(CURRENT), 'Bottom')) && size_to_here < size(lateral.PARENT.STATVAR.water_flux,1)
                    CURRENT = CURRENT.NEXT;
                    size_to_here = size_to_here + size(CURRENT.STATVAR.energy,1);
                end

                while ~(strcmp(class(CURRENT), 'Top'))
                    CURRENT = lateral3D_push_water_unconfined_aquifer_RichardsEq(CURRENT, lateral);
                    CURRENT = CURRENT.PREVIOUS;
                end
            end
        end
        
        function lateral = set_ACTIVE(lateral, i, t)
            lateral.PARENT.ACTIVE(i,1) = 0;
            if t + lateral.PARENT.IA_TIME_INCREMENT >= lateral.PARA.ia_time_next - 1e-7
                lateral.PARENT.ACTIVE(i,1) = 1;
                lateral.PARA.ia_time_next = t + lateral.PARENT.IA_TIME_INCREMENT + lateral.PARA.ia_time_increment;
                %disp(lateral.PARA.ia_time_next-floor(lateral.PARA.ia_time_next));
            end
        end
        
        function lateral = set_ia_time(lateral, t)
            lateral.PARA.ia_time_next = t;
        end
        
        
        %-------------param file generation-----
        function ground = param_file_info(ground)
            ground = param_file_info@BASE_LATERAL(ground);
            
            ground.PARA.class_category = 'LATERAL_IA';
            
            ground.PARA.options = [];
            ground.PARA.STATVAR = [];
            
            ground.PARA.default_value.hardBottom_cutoff = {0.03};
            ground.PARA.comment.hardBottom_cutoff = {'hard bottom  = no water flow if saturated and water content below [vol. water content, -]'};
            
            ground.PARA.default_value.ia_time_increment = {0.25};
            ground.PARA.comment.ia_time_increment ={'time step [days], must be multiple of LATERAL class timestep'};
        end
        
    end
    
end


