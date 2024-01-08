%========================================================================
% CryoGrid FORCING processing class
%
%
% Authors:
% S. Westermann, December 2022
%
%========================================================================

classdef calculate_precip < matlab.mixin.Copyable %makes the TRANSFORM object
    
    properties
        CONST
        PARA
        STATVAR
        TEMP
    end
    
    methods
        
        function proc = provide_PARA(proc)
            
        end
        
        
        function proc = provide_CONST(proc)

        end
        
        
        function proc = provide_STATVAR(proc)
            
        end
        
        
        function proc = finalize_init(proc, tile)

        end
        
        
        function forcing = process(proc, forcing, tile)
            
                forcing.DATA.precip = forcing.DATA.snowfall + forcing.DATA.rainfall;
        end
        
    end
    
end

