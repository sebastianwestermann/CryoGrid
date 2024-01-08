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

classdef FORCING_seb_mat < FORCING_base
    
    methods
        
        function forcing = provide_PARA(forcing)         
            % INITIALIZE_PARA  Initializes PARA structure, setting the variables in PARA.  

            forcing.PARA.filename = [];   %filename of Matlab file containing forcing data
            forcing.PARA.forcing_path = [];
            forcing.PARA.start_time = []; % start time of the simulations (must be within the range of data in forcing file)
            forcing.PARA.end_time = [];   % end time of the simulations (must be within the range of data in forcing file)
            forcing.PARA.rain_fraction = [];  %rainfall fraction assumed in sumulations (rainfall from the forcing data file is multiplied by this parameter)
            forcing.PARA.snow_fraction = [];  %snowfall fraction assumed in sumulations (snowfall from the forcing data file is multiplied by this parameter)
            forcing.PARA.heatFlux_lb = [];  % heat flux at the lower boundary [W/m2] - positive values correspond to energy gain
            forcing.PARA.airT_height = [];  % height above ground at which air temperature (and wind speed!) from the forcing data are applied.
            
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
            forcing.PARA.is_slope = 0;
            forcing.PARA.overwrite = 0;
            forcing.PARA.spatial_class = [];
            forcing.PARA.spatial_class_index = [];
            forcing.PARA.proc_class = [];
            forcing.PARA.proc_class_index = [];
            forcing = finalize_init@FORCING_base(forcing,tile);
            
            proc_functions = {'read_mat'; 'set_start_end_time'; 'assign_lb_heatflux_airT_height';	'scale_precip';	'check_and_correct'; 'check_air_pressure'};
            for i=1:size(proc_functions,1)
                proc = str2func(proc_functions{i,1});
                proc = proc();
                proc.PARA = forcing.PARA;
                proc.CONST = forcing.CONST;
                forcing = process(proc, forcing, tile);
            end
            
            %optional post-processing with dedicated classes
            if ~isempty(forcing.PARA.post_proc_class) && sum(isnan(forcing.PARA.post_proc_class_index)) == 0
                for i=1:size(forcing.PARA.post_proc_class,1)
                    post_proc_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(forcing.PARA.post_proc_class{i,1}){forcing.PARA.post_proc_class_index(i,1),1});
                    post_proc_class = finalize_init(post_proc_class, tile);
                    forcing = process(post_proc_class, forcing, tile);
                end
            end
            
            proc_functions = {'initialize_TEMP'};
            for i=1:size(proc_functions,1)
                proc = str2func(proc_functions{i,1});
                proc = proc();
                proc.PARA = forcing.PARA;
                proc.CONST = forcing.CONST;
                forcing = process(proc, forcing, tile);
            end
        end
        
        
        function forcing = interpolate_forcing(forcing, tile)
            forcing = interpolate_forcing@FORCING_base(forcing, tile);
                        
%             forcing.TEMP.rainfall = forcing.TEMP.rainfall + double(forcing.TEMP.Tair > 2) .* forcing.TEMP.snowfall;  %reassign unphysical snowfall
%             forcing.TEMP.snowfall = double(forcing.TEMP.Tair <= 2) .* forcing.TEMP.snowfall;
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