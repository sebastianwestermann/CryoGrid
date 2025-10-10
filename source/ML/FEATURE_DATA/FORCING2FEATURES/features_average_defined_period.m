classdef features_average_defined_period < matlab.mixin.Copyable
    
    properties
        PARA
        STATVAR
        CONST
        TEMP
    end

    methods
        function forcing2features = provide_PARA(forcing2features)
            forcing2features.PARA.variables = []; 
            forcing2features.PARA.start_period = [];
            forcing2features.PARA.end_period = [];
            forcing2features.PARA.number_of_periods = [];
        end


        function forcing2features = provide_CONST(forcing2features)

        end

        function forcing2features = provide_STATVAR(forcing2features)

        end

        function forcing2features = finalize_init(forcing2features, tile)
            test_year=2000;
            start_period = datenum(test_year, forcing2features.PARA.start_period(1), forcing2features.PARA.start_period(2)); 
            end_period = datenum(test_year, forcing2features.PARA.end_period(1), forcing2features.PARA.end_period(2)); 
            if start_period > end_period
                forcing2features.TEMP.offset_years = 1;
            else
                forcing2features.TEMP.offset_years = 0;
            end
        end

        function [out, data_groups] = features_from_forcing(forcing2features, features, forcing)

            test_year=2000;
            end_period1 = datenum(test_year, month(features.DATA.timestamp(1,1)), day(features.DATA.timestamp(1,1))); 
            end_period2 = datenum(test_year, forcing2features.PARA.end_period(1), forcing2features.PARA.end_period(2)); 
            if end_period2 > end_period1
                offset_years2 = 1;
            else
                offset_years2 = 0;
            end

            out = [];
            data_groups = [];
            data_group_index = 1;
            for j=1:size(features.DATA.timestamp,1)
                out_i = [];
                for v=1:size(forcing2features.PARA.variables,1)
                    for i=forcing2features.PARA.number_of_periods-1:-1:0
                        start_date = datenum(year(features.DATA.timestamp(j,1))-i-forcing2features.TEMP.offset_years-offset_years2, forcing2features.PARA.start_period(1), forcing2features.PARA.start_period(2));
                        end_date = datenum(year(features.DATA.timestamp(j,1))-i-offset_years2, forcing2features.PARA.end_period(1), forcing2features.PARA.end_period(2));
                        range = find(forcing.DATA.timeForcing>=start_date & forcing.DATA.timeForcing<end_date);
                        if strcmp(forcing2features.PARA.variables{v,1}, 'FDD')
                            out_i = [out_i sum(forcing.DATA.Tair(range,1) .* double(forcing.DATA.Tair(range,1)<0))];
                        elseif strcmp(forcing2features.PARA.variables{v,1}, 'TDD')
                            out_i = [out_i sum(forcing.DATA.Tair(range,1) .* double(forcing.DATA.Tair(range,1)>0))];
                        else
                            out_i = [out_i mean(forcing.DATA.(forcing2features.PARA.variables{v,1})(range,1))];
                        end
                    end
                    if j==1
                        %data_groups = [data_groups repmat(data_group_index,1,size(out_i,2))];
                        data_groups = [data_groups repmat(data_group_index,1,forcing2features.PARA.number_of_periods)];
                        data_group_index = data_group_index + 1;
                    end
                end
                out = [out; out_i];

            end
        end

    end
end

