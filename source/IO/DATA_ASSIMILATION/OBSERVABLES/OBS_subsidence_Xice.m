classdef OBS_subsidence_Xice < matlab.mixin.Copyable
    
    properties
        PARA
        CONST
        STATVAR
        TEMP
    end
    
    methods

        function obs = provide_PARA(obs)
          end
        
        function obs = provide_CONST(obs)
            
        end
        
        function obs = provide_STATVAR(obs)
            
        end
        
        function obs = finalize_init(obs, tile)
            obs.TEMP.initialized = 0;
            obs.TEMP.initial_ground_surface_elevation = 0;
            %this could also load the obs-files, and then automatcally
            %determine the value for the first observation -> or simply set the
            %target to be the difference between consecutive ground surface
            %levations
        end
        

        function result = observable_operator(obs, tile) 

            CURRENT = tile.TOP.NEXT;
            while ~is_ground_surface(CURRENT) && ~(strcmp(class(CURRENT), 'Bottom'))
                CURRENT = CURRENT.NEXT;
            end

            if ~obs.TEMP.initialized

                obs.TEMP.initialized = 1;
                obs.TEMP.initial_ground_surface_elevation = sum(CURRENT.STATVAR.layerThick) - CURRENT.STATVAR.XwaterIce(1,1) ./ CURRENT.STATVAR.area(1,1);

                result = 0;
            else
                ground_surface_elevation = sum(CURRENT.STATVAR.layerThick) - CURRENT.STATVAR.XwaterIce(1,1) ./ CURRENT.STATVAR.area(1,1);
                
                result = ground_surface_elevation - obs.TEMP.initial_ground_surface_elevation;
            end
        end
        

        function obs = reset_new_stratigraphy(obs, tile)
            obs = finalize_init(obs, tile);
        end
        
        
    end
end

