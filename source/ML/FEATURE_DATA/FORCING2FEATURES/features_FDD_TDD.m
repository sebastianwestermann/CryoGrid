classdef features_FDD_TDD < matlab.mixin.Copyable
    
    properties
        PARA
        STATVAR
        CONST
        TEMP
    end

    methods
        function forcing2features = provide_PARA(forcing2features)
            forcing2features.PARA.variables = []; 
            forcing2features.PARA.number_of_years = [];
            forcing2features.PARA.number_of_periods = [];
        end


        function forcing2features = provide_CONST(forcing2features)

        end

        function forcing2features = provide_STATVAR(forcing2features)

        end

        function forcing2features = finalize_init(forcing2features, tile)

        end

        function [out, data_groups] = features_from_forcing(forcing2features, features, forcing)
            out = [];
            data_groups = [];
            data_group_index = 1;
            for j=1:size(features.DATA.timestamp,1)
                out_i = [];
                for v=1:size(forcing2features.PARA.variables,1)
                    for i=forcing2features.PARA.number_of_periods:-1:1
                        start_date = datenum(year(features.DATA.timestamp(j,1))-i.*forcing2features.PARA.number_of_years, month(features.DATA.timestamp(j,1)), day(features.DATA.timestamp(j,1)));
                        end_date = datenum(year(features.DATA.timestamp(j,1))-(i-1).*forcing2features.PARA.number_of_years, month(features.DATA.timestamp(j,1)), day(features.DATA.timestamp(j,1)));
                        range = find(forcing.DATA.timeForcing>=start_date & forcing.DATA.timeForcing<end_date);
                        out_i = [out_i sum(forcing.DATA.(forcing2features.PARA.variables{v,1})(range,1) .* double(forcing.DATA.(forcing2features.PARA.variables{v,1})(range,1)>0) )];
                        out_i = [out_i sum(forcing.DATA.(forcing2features.PARA.variables{v,1})(range,1) .* double(forcing.DATA.(forcing2features.PARA.variables{v,1})(range,1)<0) )];
                        if j==1
                            data_groups = [data_groups data_group_index data_group_index+1];
                        end
                    end
                    data_group_index = data_group_index + 2;

                end
                out = [out; out_i];

            end
        end

    end
end

