classdef features_cosinus < matlab.mixin.Copyable
    
    properties
        PARA
        STATVAR
        CONST
        TEMP
    end

    methods
        function spatial2features = provide_PARA(spatial2features)
            spatial2features.PARA.variables = []; 
        end


        function spatial2features = provide_CONST(spatial2features)

        end

        function spatial2features = provide_STATVAR(spatial2features)

        end

        function spatial2features = finalize_init(spatial2features, tile)

        end

        function [out, data_groups] = features_from_spatial(spatial2features, tile, valid)
            out = [];
            data_groups = [];
            for j=1:size(spatial2features.PARA.variables, 1)
                out = [out cosd(tile.RUN_INFO.SPATIAL.STATVAR.(spatial2features.PARA.variables{j,1})(valid,:))];
                data_groups = [data_groups j];
            end
        end

    end
end

