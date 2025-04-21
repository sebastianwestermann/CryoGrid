classdef ReLU

    methods
        
        function out = propagate(af, in)
            out = max(in,0);
        end

    end
end