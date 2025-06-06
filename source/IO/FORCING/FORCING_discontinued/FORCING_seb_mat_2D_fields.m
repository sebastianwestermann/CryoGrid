%========================================================================
% CryoGrid FORCING class FORCING_seb_mat
%
% simple model forcing for GROUND classes computing the surface energy 
% balance (keyword "seb"). 
% 
% The data is obtained using the READ_FORCING_mat class. See this class for
% instructions about mat-file data format.
% 
% The mandatory forcing variables are:
%
% Tair:      Air temperature (in degree Celsius)
% Lin:       incoming long-wave radiation (in W/m2)
% Sin:       incoming short-wave radiation (in W/m2)
% rainfall:  Rainfall (in mm/day)
% snowfall:  Snowfall (in mm/day)
% q:         absolute humidity (in kg water vapor / kg air)
% p:         air pressure (OPTIONAL, in Pa)
% wind:       wind speed (in m/sec)
% 
% All forcing variables must be discretized identically, and one array of
% timestamps must be provided (t_span or timeForcing, in Matlab time / increment 1 
% corresponds to one day). 
%
% IMPORTANT POINT: the time series must be equally spaced in time, and this 
% must be really exact. When reading the timestamps from an existing data 
% set (e.g. an Excel file), rounding errors can result in small differences 
% in the forcing timestep, often less than a second off. In this case, it 
% is better to manually compile a new, equally spaced timestep in Matlab.
%
% Authors:
% S. Westermann, T. Ingeman-Nielsen, J. Scheer, October 2020
% T. Ingeman-Nielsen, October 2022
%
%========================================================================

classdef FORCING_seb_mat_2D_fields < FORCING_base
    
    methods
        
        function forcing = provide_PARA(forcing)         
            % INITIALIZE_PARA  Initializes PARA structure, setting the variables in PARA.  

            forcing.PARA.filename = [];   %filename of Matlab file containing forcing data
            forcing.PARA.forcing_path = [];
            forcing.PARA.start_time = []; % start time of the simulations (must be within the range of data in forcing file)
            forcing.PARA.end_time = [];   % end time of the simulations (must be within the range of data in forcing file)
            
            forcing.PARA.post_proc_class = [];  %optional post-processing classes
            forcing.PARA.post_proc_class_index = [];
        end
        
        
        function forcing = provide_CONST(forcing)
            forcing.CONST.Tmfw = [];
            forcing.CONST.sigma = [];            
        end
        

        function forcing = provide_STATVAR(forcing)
            
        end
                
        
        function forcing = finalize_init(forcing, tile)
            
            variables = {'precip'; 'Tair'; 'Lin'; 'Sin'; 'q'; 'wind'; 'latitude'; 'longitude'};
            [data, timestamp] = read_mat(forcing, [forcing.PARA.forcing_path forcing.PARA.filename], variables);
            
            dist_lat = abs(tile.PARA.latitude - data.latitude);
            dist_lon=abs(tile.PARA.longitude - data.longitude);
            [dist_lat, ind_lat] = sort(dist_lat);
            [dist_lon, ind_lon] = sort(dist_lon);
            dist_lat=dist_lat(1:2);
            dist_lon=dist_lon(1:2);
            ind_lat = ind_lat(1:2);
            ind_lon = ind_lon(1:2);
            weights_lat = 1 - dist_lat./sum(dist_lat);
            weights_lon = 1 - dist_lon./sum(dist_lon);
            
            variables = {'precip'; 'Tair'; 'Lin'; 'Sin'; 'q'; 'wind'};
            for i=1:size(variables,1)
                if isfield(data, variables{i,1})
                    forcing.DATA.(variables{i,1}) = weights_lon(1) .* (weights_lat(1) .* squeeze(data.(variables{i,1})(ind_lon(1),ind_lat(1),:)) + weights_lat(2) .* squeeze(data.(variables{i,1})(ind_lon(1), ind_lat(2),:))) + ...
                        weights_lon(2) .* (weights_lat(1) .* squeeze(data.(variables{i,1})(ind_lon(2), ind_lat(1), :)) + weights_lat(2) .* squeeze(data.(variables{i,1})(ind_lon(2), ind_lat(2),:)));
                end
            end
            

%             forcing.DATA.rainfall = data.rainfall.*forcing.PARA.rain_fraction;
%             forcing.DATA.snowfall = data.snowfall.*forcing.PARA.snow_fraction;
            forcing.DATA.timeForcing = timestamp;
            
           % forcing = check_and_correct(forcing); % Remove known errors
            forcing = set_start_and_end_time(forcing); % assign start/end time
            forcing = initialize_TEMP(forcing);

            %set pressure to mean pressure at corresponding altitude (international
            %altitude formula) if not provided by the forcing time series
            if ~isfield(forcing.DATA, 'p')
                altitude = tile.PARA.altitude;
                forcing.DATA.p=forcing.DATA.Tair.*0 + 1013.25.*100.*(1-0.0065./288.15.*altitude).^5.255;
            end
            
            %optional post-processing with dedicated classes
            if ~isempty(forcing.PARA.post_proc_class) && sum(isnan(forcing.PARA.post_proc_class_index)) == 0
                for i=1:size(forcing.PARA.post_proc_class,1)
                    post_proc_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(forcing.PARA.post_proc_class{i,1}){forcing.PARA.post_proc_class_index(i,1),1});
                    post_proc_class = finalize_init(post_proc_class, tile);
                    forcing = post_process(post_proc_class, forcing, tile);
                end
            end
        end
        
        
        function forcing = interpolate_forcing(forcing, tile)
            forcing = interpolate_forcing@FORCING_base(forcing, tile);
                        
            forcing.TEMP.rainfall = forcing.TEMP.rainfall + double(forcing.TEMP.Tair > 2) .* forcing.TEMP.snowfall;  %reassign unphysical snowfall
            forcing.TEMP.snowfall = double(forcing.TEMP.Tair <= 2) .* forcing.TEMP.snowfall;
        end



        %-------------param file generation-----
        function forcing = param_file_info(forcing)
            forcing = provide_PARA(forcing);

            forcing.PARA.STATVAR = [];
            forcing.PARA.class_category = 'FORCING';
            
            forcing.PARA.comment.filename = {'filename of Matlab file containing forcing data'};
            
            forcing.PARA.default_value.forcing_path = {'../CryoGridCommunity_forcing/'};
            forcing.PARA.comment.forcing_path = {'path where forcing data file is located'};
            
            forcing.PARA.comment.start_time = {'start time of the simulations (must be within the range of data in forcing file) - year month day'};
            forcing.PARA.options.start_time.name =  'H_LIST';
            forcing.PARA.options.start_time.entries_x = {'year' 'month' 'day'};
            
            forcing.PARA.comment.end_time = {'end_time time of the simulations (must be within the range of data in forcing file) - year month day'};
            forcing.PARA.options.end_time.name =  'H_LIST'; % 
            forcing.PARA.options.end_time.entries_x = {'year' 'month' 'day'};
            
            forcing.PARA.default_value.rain_fraction = {1};  
            forcing.PARA.comment.rain_fraction = {'rainfall fraction assumed in simulations (rainfall from the forcing data file is multiplied by this parameter)'};
            
            forcing.PARA.default_value.snow_fraction = {1};  
            forcing.PARA.comment.snow_fraction = {'snowfall fraction assumed in simulations (rainfall from the forcing data file is multiplied by this parameter)'};

            forcing.PARA.default_value.heatFlux_lb = {0.05};
            forcing.PARA.comment.heatFlux_lb = {'heat flux at the lower boundary [W/m2] - positive values correspond to energy gain'};
            
            forcing.PARA.default_value.airT_height = {2};  
            forcing.PARA.comment.airT_height = {'height above ground surface where air temperature from forcing data is applied'};
            
            forcing.PARA.comment.post_proc_class = {'list of postprocessing classes to modify forcing data in user-defined ways; no post-processing applied when empty'};
            forcing.PARA.options.post_proc_class.name = 'H_LIST';
            
            forcing.PARA.comment.post_proc_class_index = {''};
            forcing.PARA.options.post_proc_class_index.name = 'H_LIST';
        end
    end
end