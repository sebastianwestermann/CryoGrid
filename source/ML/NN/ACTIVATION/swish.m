classdef swish

    methods

        function out = propagate(af, in)
            out = in./(1+exp(-in));
        end

    end
end

        