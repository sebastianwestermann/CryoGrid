
classdef SNOW_MULTITILE_BUCKET_SIMPLE < BASE
    

    methods
        
        %-----initialize-----------------
        
        function ground = provide_PARA(ground)
            ground.PARA.timestep = [];
            %ground.PARA.wind_speed_class = []; %= multiplication factor for sulimation
            %ground.PARA.snowfall_factor = [];

        end
        
        function ground = provide_CONST(ground)
            
            ground.CONST.day_sec = []; %24*3600;
            
        end
        

        function ground = provide_STATVAR(ground)

            ground.STATVAR.SWE = []; % thickness of grid cells [m]
        end
        
        function ground = finalize_init(ground, tile)

           
            ground.STATVAR.SWE = 0;
            ground.STATVAR.SWE_water = 0;
       end
        
        
        %-----mandatory functions------------------------
        function ground = get_boundary_condition_u(ground, tile)
            
            ground.STATVAR.SWE = ground.STATVAR.SWE + tile.FORCING.TEMP.snowfall ./1000 ./ ground.CONST.day_sec .* ground.PARA.timestep;
            ground.STATVAR.SWE_water = ground.STATVAR.SWE_water + tile.FORCING.TEMP.rainfall ./1000 ./ ground.CONST.day_sec .* ground.PARA.timestep;
            
            ground.STATVAR.SWE = ground.STATVAR.SWE -  tile.FORCING.TEMP.sublimation ./1000 ./ ground.CONST.day_sec .* ground.PARA.timestep;
            %melt
            melt = min(ground.STATVAR.SWE, double(tile.FORCING.TEMP.melt>0) .* tile.FORCING.TEMP.melt ./ 1000 ./ground.CONST.day_sec .* ground.PARA.timestep);  %in [m], constant timestep
            ground.STATVAR.SWE = ground.STATVAR.SWE - melt;
            ground.STATVAR.SWE_water = ground.STATVAR.SWE_water + melt;
            
            %refreezing
            refreeze =  min(ground.STATVAR.SWE_water, - double(tile.FORCING.TEMP.melt<0) .* tile.FORCING.TEMP.melt ./ 1000 ./ground.CONST.day_sec .* ground.PARA.timestep);  %in [m], constant timestep
            ground.STATVAR.SWE = ground.STATVAR.SWE + refreeze;
            ground.STATVAR.SWE_water = ground.STATVAR.SWE_water - refreeze;

            ground.STATVAR.SWE = max(0, ground.STATVAR.SWE);
            ground.STATVAR.SWE_water = min(ground.STATVAR.SWE./2, ground.STATVAR.SWE_water);

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

