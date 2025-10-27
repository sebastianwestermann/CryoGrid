classdef sigmoid

    methods

        function out = propagate(af, in)
            out = 1./(1+exp(-in));
        end

    end
end


        