%========================================================================
% CryoGrid LATERAL_IA class LAT3D_WATER_OVERLAND_FLOW2
% designed to be called as stand-alone
% instead of LAT3D_WATER_UNCONFINED_AQUIFER_OVERLAND_FLOW, enabling flow
% between tiles without concurrent subsurface flow; drainage "out of the
% system" can be realized by using LAT3D_WATER_OVERLAND_FLOW after this
% class
% setting GaMa_coefficient = 0 prevents flow between tiles, so that only 1D drainage out of the system is enabled  
% S. Westermann, Apr 2024
%========================================================================


classdef LAT3D_WATER_OVERLAND_FLOW2 < BASE_LATERAL

    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        
        
        function lateral = provide_CONST(lateral)
            lateral.CONST.day_sec = []; %24 .* 3600;
        end
        
        function lateral = provide_PARA(lateral)
            lateral.PARA.GaMa_coefficient = []; %Gauckler-Manning coefficient, https://en.wikipedia.org/wiki/Manning_formula
            lateral.PARA.tortuosity = [];
            lateral.PARA.ia_time_increment = []; %0.25; %must be a multiple of the time increment of the main lateral class
        end
        
        function lateral = provide_STATVAR(lateral)
            lateral.STATVAR.surface_run_off = [];
        end

        function lateral = finalize_init(lateral, tile)
            


        end
        
        %-----time integration-------
        
        function lateral = pull(lateral, tile)
            

            lateral.PARENT.STATVAR.water_depth = 0; %depth of the free water at the surface flowing according to Gauckler-Manning
            CURRENT = lateral.PARENT.TOP.NEXT;
            CURRENT = lateral3D_pull_water_overland_flow(CURRENT, lateral);
        end
        
        function lateral = get_derivatives(lateral, tile) %no need to loop through stratigraphy, all the information is in lateral.PARENT

            lateral.PARENT.STATVAR.water_flux(1,1) = 0;
            lateral.PARENT.STATVAR.water_flux_energy(1,1) = 0;
            
            for j=1:size(lateral.PARENT.ENSEMBLE,1)
                if lateral.PARENT.PARA.connected(lateral.PARENT.STATVAR.index, lateral.PARENT.ENSEMBLE{j,1}.index)  %add if water is available at all
                    contact_length = lateral.PARENT.PARA.contact_length(lateral.PARENT.STATVAR.index, lateral.PARENT.ENSEMBLE{j,1}.index);
                    distance = lateral.PARENT.PARA.distance(lateral.PARENT.STATVAR.index, lateral.PARENT.ENSEMBLE{j,1}.index);
                    
                    % overland flow according to Gauckler-Manning
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
                        
                        gradient = -(lateral.PARENT.STATVAR.depths(1,1) - lateral.PARENT.ENSEMBLE{j,1}.depths(1,1)) ./ (distance.*lateral.PARA.tortuosity);
                        if gradient < 0  %own realization higher
                            T_water =  lateral.PARENT.STATVAR.T_water(1,1);
                            depth1 = min(lateral.PARENT.STATVAR.depths(1,1) - (lateral.PARENT.ENSEMBLE{j,1}.depths(1,1) - lateral.PARENT.ENSEMBLE{j,1}.water_depth), lateral.PARENT.STATVAR.water_depth);
                            depth2 = lateral.PARENT.ENSEMBLE{j,1}.water_depth;
                        elseif gradient > 0  %ensemble realization higher
                            T_water = lateral.PARENT.ENSEMBLE{j,1}.T_water(1,1);
                            depth1 = min(lateral.PARENT.ENSEMBLE{j,1}.depths(1,1) - (lateral.PARENT.STATVAR.depths(1,1) - lateral.PARENT.STATVAR.water_depth), lateral.PARENT.ENSEMBLE{j,1}.water_depth);
                            depth2 = lateral.PARENT.STATVAR.water_depth;
                        end
                        
                        velocity = lateral.PARA.GaMa_coefficient .* real(((depth1+depth2)./2).^(2/3) .* abs(gradient).^0.5);
                        flow =  velocity .* (depth1+depth2)./2 .* contact_length; 
                        max_flow2same_level = lateral.PARENT.STATVAR.area_flow .* lateral.PARENT.ENSEMBLE{j,1}.area_flow .* (lateral.PARENT.STATVAR.depths(1,1) - lateral.PARENT.ENSEMBLE{j,1}.depths(1,1)) ./ (lateral.PARENT.STATVAR.area_flow + lateral.PARENT.ENSEMBLE{j,1}.area_flow) ./ (lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec);

                        if gradient < 0 %own realization higher, looses water
                            flow = sign(gradient) .* min(abs(max_flow2same_level./8), min(flow, lateral.PARENT.STATVAR.max_flow ./ (lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec)));
                        else
                            flow = sign(gradient) .* min(abs(max_flow2same_level./8), min(flow, lateral.PARENT.ENSEMBLE{j,1}.max_flow ./ (lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec)));
                        end
                    
                        flow_energy = flow .* lateral.PARENT.CONST.c_w .* T_water;
                        

                        lateral.STATVAR.flow = flow;
                        
                        lateral.PARENT.STATVAR.water_flux(1,1) = lateral.PARENT.STATVAR.water_flux(1,1) + flow .* lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec;
                        lateral.PARENT.STATVAR.water_flux_energy(1,1) = lateral.PARENT.STATVAR.water_flux_energy(1,1) + flow_energy .* lateral.PARA.ia_time_increment .* lateral.PARENT.CONST.day_sec;
                    end
                    
                end
            end
            
        end
        
        
        function lateral = push(lateral, tile)
            CURRENT = lateral.PARENT.TOP.NEXT; 
            CURRENT = lateral3D_push_water_overland_flow(CURRENT, lateral);
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
        function ground = param_file_info(ground)
            ground = param_file_info@BASE_LATERAL(ground);
            
            ground.PARA.class_category = 'LATERAL_IA';
            
            ground.PARA.options = [];
            ground.PARA.STATVAR = [];
            
            ground.PARA.default_value.gradient = {'0.01'};
            ground.PARA.comment.gradient = {'gradient of surface [vertial m/ horizontal m]'};
            
            ground.PARA.default_value.GaMa_coefficient = {15};
            ground.PARA.comment.GaMa_coefficient = {'Gauckler-Manning coefficient, https://en.wikipedia.org/wiki/Manning_formula'};
            
            ground.PARA.default_value.overland_flow_contact_length = {1};
            ground.PARA.comment.overland_flow_contact_length = {'lateral contact length for overland flow = width of channel [m]'};
            
            ground.PARA.default_value.overflow_threshold_elevation = {0};
            ground.PARA.comment.overflow_threshold_elevation = {'threshold elevation, no overland flow when water level is below [m a.s.l.]'};
            
            ground.PARA.default_value.ia_time_increment = {0.25};
            ground.PARA.comment.ia_time_increment ={'time step [days]'};
        end
    end
    
end


