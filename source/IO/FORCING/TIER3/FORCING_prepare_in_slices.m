%========================================================================
% CryoGrid FORCING class FORCING_prepare_in_slices
%
%
% Authors:
% S. Westermann, December 2022, January 2024
%
%========================================================================

classdef FORCING_prepare_in_slices < FORCING_base
    
    properties
        
    end
    
    methods
        function forcing = provide_PARA(forcing)
            
            forcing.PARA.forcing_class = [];   %filename of Matlab file containing forcing data
            forcing.PARA.forcing_class_index = [];
            forcing.PARA.start_time = [];
            forcing.PARA.end_time = [];
            forcing.PARA.years_per_slice = [];
        end
        
        
        function forcing = provide_CONST(forcing)

        end
        
        
        function forcing = provide_STATVAR(forcing)
            
        end
        
        
        function forcing = finalize_init(forcing, tile)
            
            forcing.PARA.start_time = datenum(forcing.PARA.start_time(1,1), forcing.PARA.start_time(2,1), forcing.PARA.start_time(3,1));
            forcing.PARA.end_time = datenum(forcing.PARA.end_time(1,1), forcing.PARA.end_time(2,1),forcing.PARA.end_time(3,1));
            
            forcing.TEMP.time_next_slice = datenum(year(forcing.PARA.start_time)+forcing.PARA.years_per_slice, 1, 1);
            
            %load first slice
            forcing_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(forcing.PARA.forcing_class){forcing.PARA.forcing_class_index,1});
            %forcing_class.PARA.start_time = [year(forcing.PARA.start_time); month(forcing.PARA.start_time); day(forcing.PARA.start_time)];
            forcing_class.PARA.start_time = forcing.PARA.start_time;
            forcing_class.PARA.end_time = min(forcing.TEMP.time_next_slice, forcing.PARA.end_time);
            %forcing_class.PARA.end_time = [year(forcing_class.PARA.end_time); month(forcing_class.PARA.end_time); day(forcing_class.PARA.end_time)];
            forcing_class = finalize_init(forcing_class, tile);
            if forcing.TEMP.time_next_slice < forcing.PARA.end_time
                range = find(forcing_class.DATA.timeForcing>=forcing.PARA.start_time & forcing_class.DATA.timeForcing < min(forcing.TEMP.time_next_slice, forcing.PARA.end_time));
            else
                range = find(forcing_class.DATA.timeForcing>=forcing.PARA.start_time & forcing_class.DATA.timeForcing <= min(forcing.TEMP.time_next_slice, forcing.PARA.end_time));
            end
            variable_list = fieldnames(forcing_class.DATA);
            for i=1:size(variable_list,1)
                forcing.DATA.(variable_list{i,1}) = forcing_class.DATA.(variable_list{i,1})(range,1);
            end
            forcing_class.TEMP.time_next_slice = forcing.TEMP.time_next_slice;
            forcing.TEMP = forcing_class.TEMP;
            
            forcing.PARA.heatFlux_lb = forcing_class.PARA.heatFlux_lb;
            forcing.PARA.airT_height = forcing_class.PARA.airT_height;
        end
        
        
        function forcing = interpolate_forcing(forcing, tile)
            forcing = interpolate_forcing@FORCING_base(forcing, tile);
            
            forcing.TEMP.rainfall = forcing.TEMP.rainfall + double(forcing.TEMP.Tair > 2) .* forcing.TEMP.snowfall;  %reassign unphysical snowfall
            forcing.TEMP.snowfall = double(forcing.TEMP.Tair <= 2) .* forcing.TEMP.snowfall;
            
            %load new chunk
            if tile.t > forcing.DATA.timeForcing(end,1)-2 && forcing.DATA.timeForcing(end,1) < forcing.PARA.end_time
                
                old_slice = forcing.DATA;
                
                start_time = forcing.TEMP.time_next_slice;
                forcing.TEMP.time_next_slice = datenum(year(start_time)+forcing.PARA.years_per_slice, 1, 1);
                
                %load next slice
                forcing_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(forcing.PARA.forcing_class){forcing.PARA.forcing_class_index,1});
                %forcing_class.PARA.start_time = [year(start_time); month(start_time); day(start_time)];
                forcing_class.PARA.start_time = start_time;
                forcing_class.PARA.end_time = min(forcing.TEMP.time_next_slice, forcing.PARA.end_time);
                %forcing_class.PARA.end_time = [year(forcing_class.PARA.end_time); month(forcing_class.PARA.end_time); day(forcing_class.PARA.end_time)];
                forcing_class = finalize_init(forcing_class, tile);

                if forcing.TEMP.time_next_slice < forcing.PARA.end_time
                    range = find(forcing_class.DATA.timeForcing>=start_time & forcing_class.DATA.timeForcing < forcing.TEMP.time_next_slice);
                else
                    range = find(forcing_class.DATA.timeForcing>=start_time & forcing_class.DATA.timeForcing <= forcing.PARA.end_time);
                end
                variable_list = fieldnames(forcing_class.DATA);
                for i=1:size(variable_list,1)
                    forcing.DATA.(variable_list{i,1}) = forcing_class.DATA.(variable_list{i,1})(range,1);
                end                
                
                %overwrite with new data completed
                
                %cut away the old chunk except for the last 4 days (2 days
                %buffer)
                variables = fieldnames(old_slice);
                offset_index = min(size(old_slice.timeForcing,1)+1, floor(4./(old_slice.timeForcing(2)-old_slice.timeForcing(1))));
                for i=1:size(variables,1)
                    data_chunk = old_slice.(variables{i,1});
                    forcing.DATA.(variables{i,1}) = [data_chunk(end-offset_index:end,1); forcing.DATA.(variables{i,1})];
                end
            end
            
        end
        
        
        %-------------param file generation-----
%         function forcing = param_file_info(forcing)
%             forcing = provide_PARA(forcing);
% 
%             forcing.PARA.STATVAR = [];
%             forcing.PARA.class_category = 'FORCING';
%             
%             forcing.PARA.comment.nc_folder = {'folder containing nc-files with forcing data'};
%             
%             forcing.PARA.default_value.forcing_path = {'../CryoGridCommunity_forcing/'};
%             forcing.PARA.comment.forcing_path = {'path where forcing data file is located'};
%             
%             forcing.PARA.default_value.time_resolution_input = {'month'};
%             forcing.PARA.comment.time_resolution_input = {'time resolution of nc-files; three options: "month", "quarter" or "year"'};
%             
%             forcing.PARA.comment.start_time = {'start time of the simulations (must be within the range of data in forcing file) - year month day'};
%             forcing.PARA.options.start_time.name =  'H_LIST';
%             forcing.PARA.options.start_time.entries_x = {'year' 'month' 'day'};
%             
%             forcing.PARA.comment.end_time = {'end_time time of the simulations (must be within the range of data in forcing file) - year month day'};
%             forcing.PARA.options.end_time.name =  'H_LIST'; % 
%             forcing.PARA.options.end_time.entries_x = {'year' 'month' 'day'};
%             
%             forcing.PARA.default_value.rain_fraction = {1};  
%             forcing.PARA.comment.rain_fraction = {'rainfall fraction assumed in simulations (rainfall from the forcing data file is multiplied by this parameter)'};
%             
%             forcing.PARA.default_value.snow_fraction = {1};  
%             forcing.PARA.comment.snow_fraction = {'snowfall fraction assumed in simulations (rainfall from the forcing data file is multiplied by this parameter)'};
% 
%             forcing.PARA.default_value.all_rain_T = {0.5};    
%             forcing.PARA.comment.all_rain_T = {'Temperature above which all precipitation is considered as rain'};
% 
%             forcing.PARA.default_value.all_snow_T =  {-0.5};    
%             forcing.PARA.comment.all_snow_T = {'Temperature below which all precipitation is considered as snow'};
%              
%             forcing.PARA.default_value.albedo_surrounding_terrain = {0.2};
%             forcing.PARA.comment.albedo_surrounding_terrain = {'albedo of terrain in the field of view of the target location, reflecting short-wave radiation; considered static throughout the year'}; 
%             
%             forcing.PARA.default_value.heatFlux_lb = {0.05};
%             forcing.PARA.comment.heatFlux_lb = {'heat flux at the lower boundary [W/m2] - positive values correspond to energy gain'};
%             
%             forcing.PARA.default_value.airT_height = {2};  
%             forcing.PARA.comment.airT_height = {'height above ground surface where air temperature from forcing data is applied'};
%             
%             forcing.PARA.comment.post_proc_class = {'list of postprocessing classes to modify forcing data in user-defined ways; no post-processing applied when empty'};
%             forcing.PARA.options.post_proc_class.name = 'H_LIST';
%             
%             forcing.PARA.comment.post_proc_class_index = {''};
%             forcing.PARA.options.post_proc_class_index.name = 'H_LIST';
%         end
        
    end
    
end