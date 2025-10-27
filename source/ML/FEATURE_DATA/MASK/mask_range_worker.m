classdef mask_range_worker < matlab.mixin.Copyable
    
    properties
        PARA
        STATVAR
        CONST
        TEMP
    end

    methods
        function mask = provide_PARA(mask)
            mask.PARA.additive = 0; 
        end


        function mask = provide_CONST(mask)

        end

        function mask = provide_STATVAR(mask)

        end

        function mask = finalize_init(mask, tile)

        end

        function valid = apply_mask2(mask, tile, valid)
            valid_vector = zeros(size(tile.RUN_INFO.SPATIAL.STATVAR.latitude,1),1);
            valid_vector(tile.PARA.range,1) = 1;
            valid_vector = logical(valid_vector);
            if mask.PARA.additive
                valid = valid | valid_vector;
            else
                valid = valid & valid_vector;
            end
            
        end

    end
end

