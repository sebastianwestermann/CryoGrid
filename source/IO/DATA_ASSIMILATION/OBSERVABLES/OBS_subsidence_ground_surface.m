classdef OBS_subsidence_ground_surface < matlab.mixin.Copyable
    
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
            obs.TEMP.current_year = 1e20; %set some unrealistic year so that ice height is set with the first observation
            obs.TEMP.ice_height = 0;
        end
        
        function result = observable_operator(obs, tile) 

            if year(tile.t) == obs.TEMP.current_year

                CURRENT = tile.TOP.NEXT;
                while ~is_ground_surface(CURRENT) && ~(strcmp(class(CURRENT), 'Bottom'))
                    CURRENT = CURRENT.NEXT;
                end
                current_ice_height = 0;
                while ~(strcmp(class(CURRENT), 'Bottom'))
                    current_ice_height = current_ice_height + sum(CURRENT.STATVAR.ice ./ CURRENT.STATVAR.area);
                    CURRENT = CURRENT.NEXT;
                end
                
                result = (current_ice_height - obs.TEMP.ice_height) .* 0.08; 
            else
                
                obs.TEMP.current_year = year(tile.t);
                %total ice
                result = 0;
                
                CURRENT = tile.TOP.NEXT;
                while ~is_ground_surface(CURRENT) && ~(strcmp(class(CURRENT), 'Bottom'))
                    CURRENT = CURRENT.NEXT;
                end
                
                obs.TEMP.ice_height = 0;
                while ~(strcmp(class(CURRENT), 'Bottom'))
                    obs.TEMP.ice_height = obs.TEMP.ice_height + sum(CURRENT.STATVAR.ice ./ CURRENT.STATVAR.area);
                    CURRENT = CURRENT.NEXT;
                end
                
            end
        end
        
        function obs = reset_new_stratigraphy(obs, tile)
            obs = finalize_init(obs, tile);
        end
        
        
    end
end

