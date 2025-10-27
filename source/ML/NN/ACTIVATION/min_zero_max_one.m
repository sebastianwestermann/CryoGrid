classdef min_zero_max_one

    methods

        function out = propagate(af, in)
            out = min(1, max(0, in));
        end

    end
end


        