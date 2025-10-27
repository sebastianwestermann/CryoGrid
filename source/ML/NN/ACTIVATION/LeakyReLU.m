classdef LeakyReLU

    methods

        function out = propagate(af, in)
            out = in;
            out(in<0) = 0.01 .*  out(in<0); 
        end

    end
end