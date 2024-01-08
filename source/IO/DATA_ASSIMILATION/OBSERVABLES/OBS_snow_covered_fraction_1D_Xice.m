classdef OBS_snow_covered_fraction_1D_Xice < matlab.mixin.Copyable
    
    properties
        PARA
        CONST
        STATVAR
    end
    
    methods

        function obs = provide_PARA(obs)
            obs.PARA.reference_SWE = []; 
          end
        
        function obs = provide_CONST(obs)
            
        end
        
        function obs = provide_STATVAR(obs)
            
        end
        
        function obs = finalize_init(obs, tile)
        end
        
        function result = observable_operator(obs, tile) 
            result = 0;
            CURRENT = tile.TOP.NEXT;
            while ~(strcmp(class(CURRENT), 'Bottom'))
                class_name = class(CURRENT);
                if strcmp(class_name(1:4), 'SNOW')
                    result = 1;
                    break
                end
                if any(strcmp(properties(CURRENT), 'CHILD'))
                    if isempty(CURRENT.CHILD) || CURRENT.CHILD == 0
                        ice_height = CURRENT.STATVAR.Xice(1,1) ./ CURRENT.STATVAR.area(1,1);
                        result = min(1,max(0, ice_height./ obs.PARA.reference_SWE));  
                    else
                        ice_height = (CURRENT.CHILD.STATVAR.ice + CURRENT.STATVAR.Xice(1,1)) ./ CURRENT.STATVAR.area(1,1);
                        result = min(1,max(0, ice_height./ obs.PARA.reference_SWE));    
                    end
                    break
                end

                CURRENT = CURRENT.NEXT;
            end
            
        end
        
        function obs = reset_new_stratigraphy(obs, tile)
            
        end
        
    end
end

