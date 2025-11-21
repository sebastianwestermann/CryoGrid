%========================================================================
% CryoGrid DATA_PROVIDER class GET_VARIABLES_DEM
% Loads/calculates spatial variables similar to dem_DEM
% R.B. Zweigel, November 2025
%========================================================================

classdef GET_VARIABLES_DEM < DEM_BASE

    properties
        PARENT
    end

    methods
        function dem = provide_PARA(dem)
            dem.PARA.variables = []; % altitude or altitude, slope_angle, aspect or altitude, slope_angle, aspect, horizon_angles

            dem.PARA.number_of_horizon_bins = []; %multiples of 4!
            dem.PARA.DEM_folder = [];
            dem.PARA.DEM_filename = [];
            dem.PARA.reproject2utm = 1; %select 1 when using a geographic coordinate system and computing more than just altitude; select 0 to speed up altitde computation in big DEMs

        end

        function dem = provide_STATVAR(dem)

        end

        function dem = provide_CONST(dem)

        end

        function dem = finalize_init(dem)

        end

        function dem = load_data(dem)
            % Same as in dem_DEM
            % dem.STATVAR.latitude = dem.PARA.latitude;
            % dem.STATVAR.longitude = dem.PARA.longitude;
            % dem.STATVAR.area = dem.PARA.area;
            % dem.STATVAR.altitude = 0;
            % dem.STATVAR.slope_angle = 0;
            % dem.STATVAR.aspect = 0;
            % dem.STATVAR.skyview_factor = 1;
            % dem.STATVAR.horizon_bins = 0;
            % dem.STATVAR.horizon_angles = 0;

            dem = read_DEM_raster(dem);
            temp = dem.PARENT.RUN_INFO.SPATIAL.STATVAR;
            temp.horizon_angles = nan(length(temp.latitude),dem.PARA.number_of_horizon_bins+1); % These are not provided by ENSEMBLE_OF_POINTS
            temp.horizon_bins = temp.horizon_angles;

            [C, ia, ic] = unique([temp.latitude temp.longitude], 'rows');

            for i = 1:length(C)
                dem.PARA.latitude = temp.latitude(ia(i));
                dem.PARA.longitude = temp.longitude(ia(i));

                dem = project_target_coordinates(dem);
                dem = compute_global_offset_from_north(dem);

                for j=1:size(dem.PARA.variables,1)
                    a = str2func(['get_' dem.PARA.variables{j,1}]);
                    dem = a(dem);
                    if strcmp(dem.PARA.variables{j,1},'horizon_angles')
                        temp.horizon_angles(ic==i,:) = repmat(dem.STATVAR.horizon_angles', sum(ic==i),1);
                        temp.horizon_bins(ic==i,:) = repmat(dem.STATVAR.horizon_bins', sum(ic==i),1);
                    else
                        temp.(dem.PARA.variables{j,1})(ic==i,:) = dem.STATVAR.(dem.PARA.variables{j,1})';
                    end
                end
            end
            dem.PARENT.RUN_INFO.SPATIAL.STATVAR = temp;
        end

        function dem = get_altitude(dem)
            dem = get_altitude_base(dem);
        end

        function dem = get_slope_angle(dem)
            dem = get_slope_angle_base(dem);
        end

        function dem = get_aspect(dem)
            dem = get_aspect_base(dem);
        end

        function dem = get_horizon_angles(dem)
            dem = get_horizon_angles_single_point(dem);
        end


    end
end

