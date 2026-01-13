classdef SPATIAL2ENSEMBLE < matlab.mixin.Copyable

    properties
        PARA
        CONST
        STATVAR
        TEMP
        PARENT
    end
    
    methods
        
        function ensemble = provide_PARA(ensemble)
            
        end
        
        function ensemble = provide_CONST(ensemble)

        end
        
        function ensemble = provide_STATVAR(ensemble)

        end 
                
        function ensemble = finalize_init(ensemble, tile)
            ensemble.PARA.id_variable = [];
        end

        function ensemble_class = generate_ensemble_from_existing(ensemble_class, proj)
            
            ensemble_class.STATVAR = ensemble_class.PARENT.RUN_INFO.SPATIAL.STATVAR;
            
        end

    end
end

