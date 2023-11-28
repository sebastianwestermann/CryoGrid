classdef q_constant_RH < matlab.mixin.Copyable
    
    properties
        PARA
        STATVAR
        CONST
    end
    
    methods
        
        function bias_correct = provide_PARA(bias_correct)

        end
        
        function bias_correct = provide_CONST(bias_correct)
            
        end
        
        function bias_correct = provide_STATVAR(bias_correct)

        end
        
        function bias_correct = finalize_init(bias_correct, tile)

        end
        
        function [forcing_corrected, overlap_pairs] = bias_correct(bias_correct, forcing, carrier_class, tile)
            forcing_corrected = carrier_class.DATA.q;
            forcing_T_old = carrier_class.DATA.Tair + 273.15;
            forcing_T_new = forcing.DATA.Tair + 273.15;
            
            range = find(forcing_T_old>=273.15);
            forcing_corrected(range) = forcing_corrected(range) .* exp(17.62.*(forcing_T_new(range)-273.15)./(243.12-273.15+forcing_T_new(range))) ./ exp(17.62.*(forcing_T_old(range)-273.15)./(243.12-273.15+forcing_T_old(range)));
            range = find(forcing_T_old<273.15);
            forcing_corrected(range) = forcing_corrected(range) .*  exp(22.46.*(forcing_T_new(range)-273.15)./(272.61-273.15+forcing_T_new(range))) ./ exp(22.46.*(forcing_T_old(range)-273.15)./(272.61-273.15+forcing_T_old(range)));
            overlap_pairs = NaN;
        end
        
        
    end
end