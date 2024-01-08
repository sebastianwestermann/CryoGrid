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

classdef FORCING_base <  matlab.mixin.Copyable
    
    properties
        PARA
        CONST
        STATVAR
        TEMP
        SPATIAL
        TIME
        DATA
    end
    
    methods
        
        function forcing = provide_PARA(forcing)         
            forcing.PARA.proc_class = [];  
            forcing.PARA.proc_class_index = [];
            forcing.PARA.spatial_class = [];
            forcing.PARA.spatial_class_index = [];
        end
        
        
        function forcing = provide_CONST(forcing)
         
        end
        

        function forcing = provide_STATVAR(forcing)
            
        end
                
        
        function forcing = finalize_init(forcing, tile)
            if ~isempty(forcing.PARA.spatial_class) && sum(isnan(forcing.PARA.spatial_class)) == 0
                forcing.SPATIAL = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(forcing.PARA.spatial_class){forcing.PARA.spatial_class_index,1});
                forcing.SPATIAL = finalize_init(forcing.SPATIAL);
            else
                forcing.SPATIAL.STATVAR.altitude = tile.PARA.altitude;
                forcing.SPATIAL.STATVAR.slope_angle = tile.PARA.slope_angle;
                forcing.SPATIAL.STATVAR.aspect = tile.PARA.aspect;
                forcing.SPATIAL.STATVAR.skyview_factor = tile.PARA.skyview_factor;                
                forcing.SPATIAL.STATVAR.horizon_bins = tile.PARA.horizon_bins;
                forcing.SPATIAL.STATVAR.horizon_angles = tile.PARA.horizon_angles;
                forcing.SPATIAL.STATVAR.latitude = tile.PARA.latitude;
                forcing.SPATIAL.STATVAR.longitude = tile.PARA.longitude;
            end
            
            if ~isempty(forcing.PARA.proc_class) && sum(isnan(forcing.PARA.proc_class_index)) == 0
                for i=1:size(forcing.PARA.proc_class,1)
                    proc_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(forcing.PARA.proc_class{i,1}){forcing.PARA.proc_class_index(i,1),1});
                    proc_class = finalize_init(proc_class, tile);
                    forcing = process(proc_class, forcing, tile);
                end
            end
        end
        
        
        function forcing = interpolate_forcing(forcing, tile)
            t = tile.t;
            
            posit = floor((t-forcing.DATA.timeForcing(1,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1)))+1;
            
            variables = fieldnames(forcing.TEMP);
            t_weight = (t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1)); % current distance from last timestep (0,1)
            
            for i = 1:length(variables)
                if ~strcmp(variables{i},'t') && isfield(forcing.DATA, variables{i})
                    forcing.TEMP.(variables{i}) = forcing.DATA.(variables{i})(posit,:)+(forcing.DATA.(variables{i})(posit+1,:)-forcing.DATA.(variables{i})(posit,:)).*t_weight;
                end
            end
            forcing.TEMP.t = t;
        end

        function forcing = adjust_forcing(forcing, tile)
            
            %optional post-processing with dedicated classes
            if ~isempty(forcing.PARA.post_proc_class) && sum(isnan(forcing.PARA.post_proc_class_index)) == 0
                for i=1:size(forcing.PARA.post_proc_class,1)
                    post_proc_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(forcing.PARA.post_proc_class{i,1}){forcing.PARA.post_proc_class_index(i,1),1});
                    forcing = process(post_proc_class, forcing, tile);
                end
            end
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