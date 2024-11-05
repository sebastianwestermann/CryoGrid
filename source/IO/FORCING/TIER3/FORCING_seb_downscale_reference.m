%========================================================================
% CryoGrid FORCING class FORCING_seb_downscale_CMIP
%
% forcing for GROUND classes computing the surface energy 
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
% Authors:
% S. Westermann, January 2024
% replaces FORCING_downscale_slope_seb
%
%========================================================================

classdef FORCING_seb_downscale_reference < FORCING_base_carrier_reference
    
    methods
        
        function forcing = provide_PARA(forcing)         

            forcing.PARA.carrier_forcing_class = [];
            forcing.PARA.carrier_forcing_class_index = [];
            forcing.PARA.offset_from_GMT_carrier = []; %in hours
            forcing.PARA.reference_forcing_class = [];   
            forcing.PARA.reference_forcing_class_index = [];
            forcing.PARA.offset_from_GMT_reference = []; %in hours
            
            forcing.PARA.detrend = [];
            
            forcing.PARA.overlap_interval = []; %interval for which the list of overlap-pairs is created
            forcing.PARA.selection_period_before = [];
            forcing.PARA.selection_period_after = [];
 
            forcing.PARA.start_time = []; % start time of the simulations (must be within the range of data in forcing file)
            forcing.PARA.end_time = [];   % end time of the simulations (must be within the range of data in forcing file)
            forcing.PARA.rain_fraction = [];  %rainfall fraction assumed in sumulations (rainfall from the forcing data file is multiplied by this parameter)
            forcing.PARA.snow_fraction = [];  %snowfall fraction assumed in sumulations (snowfall from the forcing data file is multiplied by this parameter)
            forcing.PARA.all_rain_T = [];
            forcing.PARA.all_snow_T = [];
            
            forcing.PARA.heatFlux_lb = [];  % heat flux at the lower boundary [W/m2] - positive values correspond to energy gain
            forcing.PARA.airT_height = [];  % height above ground at which air temperature (and wind speed!) from the forcing data are applied

            forcing.PARA.post_proc_class = [];  %optional post-processing classes
            forcing.PARA.post_proc_class_index = [];
            
            %             forcing.PARA.variables = [];
%             forcing.PARA.relative_correction = [];
        end
        
        
        function forcing = provide_CONST(forcing)
         
        end
        

        function forcing = provide_STATVAR(forcing)
            
        end
                
        
        function forcing = finalize_init(forcing, tile)
            
            monthly_offset_class = month_offset();
            monthly_offset_class.PARA.relative_correction = 0;
            monthly_offset_class.PARA.trend_reference_to_carrier = 1;
            tile.RUN_INFO.PPROVIDER.CLASSES.month_offset{1,1} = monthly_offset_class;

            monthly_offset_class = month_offset();
            monthly_offset_class.PARA.relative_correction = 1;
            monthly_offset_class.PARA.trend_reference_to_carrier = 1;
            tile.RUN_INFO.PPROVIDER.CLASSES.month_offset{2,1} = monthly_offset_class;
            
            forcing.PARA.start_overlap = [forcing.PARA.overlap_interval(1); 1; 1];
            forcing.PARA.end_overlap = [forcing.PARA.overlap_interval(2); 12; 31];
            forcing.PARA.is_slope = 0;
            forcing.PARA.overwrite = 0;
            forcing.PARA.save_transform2file = 0;
            forcing.PARA.number_of_quantiles = 20;
            forcing.PARA.spatial_class = [];
            forcing.PARA.spatial_class_index = [];
            forcing.PARA.proc_class = [];
            forcing.PARA.proc_class_index = [];

            carrier_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(forcing.PARA.carrier_forcing_class){forcing.PARA.carrier_forcing_class_index,1});
            carrier_class = finalize_init(carrier_class, tile);
            carrier_class.DATA.timeForcing = carrier_class.DATA.timeForcing - forcing.PARA.offset_from_GMT_carrier ./ 24;
%             carrier_class.PARA.start_time = datenum(forcing.PARA.overlap_interval(1), 1,1);
%             carrier_class.PARA.end_time = datenum(forcing.PARA.overlap_interval(2)+1, 1,1)-0.01;
%             proc = clip2start_end_time(); %clip carrier to overlap interval
%             carrier_class = process(proc, carrier_class, tile);
            forcing.CARRIER = carrier_class;
            
            reference_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(forcing.PARA.reference_forcing_class){forcing.PARA.reference_forcing_class_index,1});
            reference_class.PARA.start_time = datenum(forcing.PARA.overlap_interval(1,1),1,1);
            reference_class.PARA.end_time = datenum(forcing.PARA.overlap_interval(2,1),12,31);
            reference_class = finalize_init(reference_class, tile);
            reference_class.DATA.timeForcing = reference_class.DATA.timeForcing - forcing.PARA.offset_from_GMT_reference ./ 24;
         
            forcing.REFERENCE = reference_class;
            
            reload_REFERENCE = 0;
            
            forcing = finalize_init@FORCING_base(forcing, tile);

            
            
            if forcing.PARA.detrend
                proc_functions = {'set_start_end_time'; 'detrend_reference_and_carrier_monthly';  'extend_carrier_random_months'; 'fit_TRANSFORM'; 'apply_TRANSFORM'; ...
                    'retrend_forcing_monthly'; 'correct_Sin_time_series'; 'set_min_max'; 'calculate_q_from_RH'; 'split_precip_snow_rain'; 'scale_precip';  'assign_lb_heatflux_airT_height'; 'check_air_pressure'};
            else
                    proc_functions = {'set_start_end_time'; 'fit_TRANSFORM'; 'extend_carrier_random_months'; 'apply_TRANSFORM'; ...
                       'correct_Sin_time_series'; 'set_min_max'; 'calculate_q_from_RH'; 'split_precip_snow_rain'; 'scale_precip'; 'assign_lb_heatflux_airT_height'; 'check_air_pressure'};
            end

            for i=1:size(proc_functions,1)
                proc = str2func(proc_functions{i,1});
                proc = proc();
                proc.PARA = forcing.PARA;
                proc.CONST = forcing.CONST;
                if strcmp(proc_functions{i,1}, 'detrend_reference_and_carrier_monthly') || strcmp(proc_functions{i,1}, 'retrend_forcing_monthly')
                    proc.PARA.variables = {'Tair'; 'precip'; 'Lin'; 'wind'; 'Sin'; 'RH'};
                    proc.PARA.relative_correction =[0; 1; 0; 1;	1; 1];
                end
%                  if strcmp(proc_functions{i,1}, 'detrend_reference_and_carrier_monthly') 
%                     proc.PARA.detrend_interval = forcing.PARA.overlap_interval;
%                 end
                if strcmp(proc_functions{i,1}, 'fit_TRANSFORM')
                    proc.PARA.variables	={'Tair'; 'precip'; 'Lin'; 'RH'; 'Sin'; 'wind'};
                    proc.PARA.overlap_target_interval = {'month';'month'; 'month'; 'month'; 'month'; 'month'};
                    proc.PARA.transform_class = {'month_offset'; 'month_offset'; 'month_offset'; 'month_offset'; 'month_offset'; 'month_offset'};
                    proc.PARA.transform_class_index	= [1; 2; 1;	2; 2; 2];
                end
                if strcmp(proc_functions{i,1},  'set_min_max')
                    proc.PARA.variables	= {'wind'; 'RH'; 'Sin'; 'precip'; 'Lin'};
                    proc.PARA.minimum = [0.5; 0.1; 0; 0; 25];
                    proc.PARA.maximum =[NaN; 1; NaN; NaN; NaN];
                end
                %reinitialze reference in case
                if  strcmp(proc_functions{i,1},  'apply_TRANSFORM') && (forcing.PARA.end_time < forcing.REFERENCE.DATA.timeForcing(1,1) || forcing.PARA.start_time > forcing.REFERENCE.DATA.timeForcing(end,1)) 
                    forcing.REFERENCE.PARA.start_time = forcing.PARA.start_time;
                    forcing.REFERENCE.PARA.end_time = forcing.PARA.end_time;
                    forcing.REFERENCE = finalize_init(forcing.REFERENCE, tile);
                end
                
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