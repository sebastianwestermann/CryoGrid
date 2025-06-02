classdef features_annual_time_series < matlab.mixin.Copyable
    
    properties
        PARA
        STATVAR
        CONST
        TEMP
    end

    methods
        function forcing2features = provide_PARA(forcing2features)
            forcing2features.PARA.variables = []; 
        end


        function forcing2features = provide_CONST(forcing2features)

        end

        function forcing2features = provide_STATVAR(forcing2features)

        end

        function forcing2features = finalize_init(forcing2features, tile)

        end

        function out = features_from_forcing(forcing2features, features, forcing)

            out=[];
            for v=1:size(forcing2features.PARA.variables,1)
                for j=1:size(features.DATA.timestamp,1)-1
                    start_date = datenum(year(features.DATA.timestamp(j,1)), month(features.DATA.timestamp(j,1)), day(features.DATA.timestamp(j,1)));
                    end_date = datenum(year(features.DATA.timestamp(j,1))+1, month(features.DATA.timestamp(j,1)), day(features.DATA.timestamp(j,1)));
                    range = find(forcing.DATA.timeForcing>=start_date & forcing.DATA.timeForcing<end_date);
                    out = [out mean(forcing.DATA.(forcing2features.PARA.variables{v,1})(range,1))];
                end
            end
        
        end

    end
end

