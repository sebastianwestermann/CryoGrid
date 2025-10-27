%TRANSFORM class q_constant_RH

% S. Westermann Dec 2023

classdef q_constant_RH < matlab.mixin.Copyable
    
    properties
        PARA
        STATVAR
        CONST
        TEMP
    end
    
    methods
        
        function transform = provide_PARA(transform)

        end
        
        function transform = provide_CONST(transform)
            
        end
        
        function transform = provide_STATVAR(transform)

        end
        
        function transform = finalize_init(transform, tile)

        end
        
        function transform = fit_transform(transform, forcing, tile)
            
        end
        
        function forcing_corrected = apply_transform(transform, forcing, tile)
            
            forcing_corrected = forcing.CARRIER.DATA.q;
            forcing_T_old =  forcing.CARRIER.DATA.Tair + 273.15;
            forcing_T_new = forcing.DATA.Tair + 273.15;
            
            range = find(forcing_T_old>=273.15);
            forcing_corrected(range) = forcing_corrected(range) .* exp(17.62.*(forcing_T_new(range)-273.15)./(243.12-273.15+forcing_T_new(range))) ./ exp(17.62.*(forcing_T_old(range)-273.15)./(243.12-273.15+forcing_T_old(range)));
            range = find(forcing_T_old<273.15);
            forcing_corrected(range) = forcing_corrected(range) .*  exp(22.46.*(forcing_T_new(range)-273.15)./(272.61-273.15+forcing_T_new(range))) ./ exp(22.46.*(forcing_T_old(range)-273.15)./(272.61-273.15+forcing_T_old(range)));
        end
        
        
    end
end