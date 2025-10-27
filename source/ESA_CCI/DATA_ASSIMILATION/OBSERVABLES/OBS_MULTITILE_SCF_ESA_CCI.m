classdef OBS_MULTITILE_SCF_ESA_CCI < matlab.mixin.Copyable
    
    properties
        PARA
        CONST
        STATVAR
    end
    
    methods

        function obs = provide_PARA(obs)
          end
        
        function obs = provide_CONST(obs)
            
        end
        
        function obs = provide_STATVAR(obs)
            
        end
        
        function obs = finalize_init(obs, tile)
        end
        
        function result = observable_operator(obs, tile) 
            
            SWE = reshape(tile.SUBSURFACE_CLASS.STATVAR.SWE, tile.PARA.number_of_realizations, tile.ENSEMBLE.PARA.subgrid_ensemble_size, tile.ENSEMBLE.PARA.grid_ensemble_size);
            result=squeeze(mean(double(SWE>1e-5),2));
            result = result(:)';
 
        end
        
        function obs = reset_new_stratigraphy(obs, tile)
            
        end
        
    end
end

