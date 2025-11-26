
classdef SNOW_MULTITILE_BUCKET_GROUND < BASE
    

    methods
        
        %-----initialize-----------------
        
        function ground = provide_PARA(ground)
            ground.PARA.timestep = [];
            ground.PARA.heat_capacity = [];
            ground.PARA.depth_below_ground = [];
            ground.PARA.conductivity_ground = [];
            ground.PARA.snow_density = [];
            ground.PARA.conductivity_snow = [];
        end
        
        function ground = provide_CONST(ground)
            
            ground.CONST.day_sec = []; %24*3600;
            ground.CONST.L_f = [];
            
        end
        

        function ground = provide_STATVAR(ground)

            ground.STATVAR.SWE = []; % thickness of grid cells [m]
        end
        
        function ground = finalize_init(ground, tile)

           
            ground.STATVAR.SWE = 0;
            ground.STATVAR.SWE_water = 0;
            ground.STATVAR.T = 0;
            ground.STATVAR.energy = 0;
            ground.TEMP.SWE2SD = 920 ./ ground.PARA.snow_density;
       end
        
        
        %-----mandatory functions------------------------
        function ground = get_boundary_condition_u(ground, tile)
            
            melting = ground.STATVAR.SWE_water > 0;
            ground.STATVAR.SWE = ground.STATVAR.SWE + tile.FORCING.TEMP.snowfall ./1000 ./ ground.CONST.day_sec .* ground.PARA.timestep;
            ground.STATVAR.SWE_water = ground.STATVAR.SWE_water + tile.FORCING.TEMP.rainfall ./1000 ./ ground.CONST.day_sec .* ground.PARA.timestep;
            
            ground.STATVAR.SWE = ground.STATVAR.SWE -  tile.FORCING.TEMP.sublimation ./1000 ./ ground.CONST.day_sec .* ground.PARA.timestep;
            %ground heat fllux 
            freeze_melt = tile.FORCING.TEMP.melt ./ 1000 ./  ground.CONST.day_sec + ground.STATVAR.T .* ground.PARA.conductivity_ground ./ ground.PARA.depth_below_ground ./ ground.CONST.L_f;

            %melt
            melt = min(ground.STATVAR.SWE, double(freeze_melt>0) .* freeze_melt .* ground.PARA.timestep);  %in [m], constant timestep
            ground.STATVAR.SWE = ground.STATVAR.SWE - melt;
            ground.STATVAR.SWE_water = ground.STATVAR.SWE_water + melt;
            
            %refreezing
            refreeze =  min(ground.STATVAR.SWE_water, - double(freeze_melt<0) .* freeze_melt .* ground.PARA.timestep);  %in [m], constant timestep
            ground.STATVAR.SWE = ground.STATVAR.SWE + refreeze;
            ground.STATVAR.SWE_water = max(0, ground.STATVAR.SWE_water - refreeze);

            ground.STATVAR.SWE = max(0, ground.STATVAR.SWE);
            ground.STATVAR.SWE_water = min(ground.STATVAR.SWE ./ 0.5, ground.STATVAR.SWE_water);

            %update of subsurface ground temperature  
            ground.STATVAR.energy = ground.STATVAR.energy + ground.PARA.timestep .* (-double(~melting & ground.STATVAR.SWE > 0) .* (ground.STATVAR.T - tile.FORCING.TEMP.surfT) .* ground.PARA.conductivity_ground .* ground.PARA.conductivity_snow  ./ ...
                (ground.PARA.conductivity_ground .* ground.STATVAR.SWE .* ground.TEMP.SWE2SD + ground.PARA.conductivity_snow .* ground.PARA.depth_below_ground) - ...
                double(melting & ground.STATVAR.SWE > 0) .* ground.STATVAR.T .* ground.PARA.conductivity_ground ./ ground.PARA.depth_below_ground - ...
                double(ground.STATVAR.SWE <= 0) .* (ground.STATVAR.T - tile.FORCING.TEMP.surfT) .* ground.PARA.conductivity_ground ./ ground.PARA.depth_below_ground);
            ground.STATVAR.T = ground.STATVAR.energy ./ (2 .* ground.PARA.heat_capacity .* ground.PARA.depth_below_ground);

        end
                
        function ground = get_boundary_condition_l(ground,  tile)

        end
        
        %calculate spatial derivatives
        function ground = get_derivatives_prognostic(ground, tile)
            
            
        end
        
        %prognostic step - integrate prognostic variables in time
        function ground = advance_prognostic(ground, tile)
             
     
        end
        
        %diagnostic step - compute diagnostic variables
        function ground = compute_diagnostic(ground, tile)
            
        end
        
        %triggers
        function ground = check_trigger(ground, tile)
         
        end
     
        
    end
    
end

