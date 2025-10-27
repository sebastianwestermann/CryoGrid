classdef normalize_multicolumn_std_reshape
 < matlab.mixin.Copyable
    
    properties
        PARA
        STATVAR
        CONST
    end

    methods
        function in_out = provide_PARA(in_out)

        end

        function in_out = provide_CONST(in_out)

        end
        function in_out = provide_STATVAR(in_out)

        end

        function in_out = finalize_init(in_out, var, tile)
            in_out.STATVAR.var = var;
            in_out.STATVAR.var_mean = mean(var,1);
            in_out.STATVAR.var_std = std(var,[],1);
            in_out.STATVAR.var_size = size(var);          
        end


        function out = real_world2NN_normalize_reshape(in_out, in, tile)
            out = (in - in_out.STATVAR.var_mean) ./ in_out.STATVAR.var_std;
            out = out(:);
        end

        function out = NN2real_world_normalize(in_out, in, tile)
            in = reshape(in, in_out.STATVAR.var_size);
            out = in .* in_out.STATVAR.var_std + in_out.STATVAR.var_mean;
        end

        function out = real_world2NN_no_normalize(in_out, in, tile)
            out = in(:);
        end

        function out = NN2real_world_no_normalize(in_out, in, tile)
            in = reshape(in, in_out.STATVAR.var_size);
        end
    end
end

