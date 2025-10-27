%========================================================================
% CryoGrid FORCING class FORCING_bias_correct_with_measurements
%
% Authors:
% S. Westermann, T. Ingeman-Nielsen, J. Scheer, October 2020
% T. Ingeman-Nielsen, October 2022
%
%========================================================================

classdef FORCING_base_carrier_reference <  FORCING_base
    
    properties
       CARRIER
       REFERENCE
       TRANSFORM %must be a cell array containing different classes
    end
    
    
    methods
        
        function forcing = provide_PARA(forcing)         

            forcing.PARA.carrier_forcing_class = [];  %must be filled
            forcing.PARA.carrier_forcing_class_index = [];
            forcing.PARA.offset_from_GMT_carrier = []; %in hours
            forcing.PARA.reference_forcing_class = [];    %can be empty if no modeification is applied
            forcing.PARA.reference_forcing_class_index = [];
            forcing.PARA.offset_from_GMT_reference = []; %in hours
            
            forcing.PARA.proc_class = [];
            forcing.PARA.proc_class_index = [];
            forcing.PARA.spatial_class = [];
            forcing.PARA.spatial_class_index = [];
            
            forcing.PARA.save_transform2file = 0;
            forcing.PARA.transform_name = [];
          
        end
        
        
        function forcing = provide_CONST(forcing)
      
        end
        

        function forcing = provide_STATVAR(forcing)
            
        end
                
        
        function forcing = finalize_init(forcing, tile)
                       
            %load carrier and reference classes
            if ~isempty(forcing.PARA.carrier_forcing_class) && sum(isnan(forcing.PARA.carrier_forcing_class))==0
                carrier_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(forcing.PARA.carrier_forcing_class){forcing.PARA.carrier_forcing_class_index,1});
                carrier_class = finalize_init(carrier_class, tile);
                carrier_class.DATA.timeForcing = carrier_class.DATA.timeForcing - forcing.PARA.offset_from_GMT_carrier ./ 24;
                forcing.CARRIER = carrier_class;
            end
            if ~isempty(forcing.PARA.reference_forcing_class) && sum(isnan(forcing.PARA.reference_forcing_class))==0
                reference_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(forcing.PARA.reference_forcing_class){forcing.PARA.reference_forcing_class_index,1});
                reference_class = finalize_init(reference_class, tile);
                reference_class.DATA.timeForcing = reference_class.DATA.timeForcing - forcing.PARA.offset_from_GMT_reference ./ 24;
                forcing.REFERENCE = reference_class;
            end
            
            forcing = finalize_init@FORCING_base(forcing, tile);

            if forcing.PARA.save_transform2file
                transform = forcing.TRANSFORM;
                spatial = forcing.SPATIAL;
                save([tile.PARA.result_path tile.PARA.run_name '/' forcing.PARA.transform_name '.mat'], 'transform', 'spatial')
            end

            forcing.CARRIER = [];
            forcing.REFERENCE = [];
            %only TRANSFORM is filled at the end, unless data, etc is
            %populated by a PROCES class -> TRANSFORM can be passsed on to
            %another class which uses it
        end
        
        
        function forcing = interpolate_forcing(forcing, tile)
            
            forcing = interpolate_forcing@FORCING_base(forcing, tile);
                        
        end

        
        

        
        


        %-------------param file generation-----
        function forcing = param_file_info(forcing)
            forcing = provide_PARA(forcing);

%             forcing.PARA.STATVAR = [];
%             forcing.PARA.class_category = 'FORCING';
%             
%             forcing.PARA.comment.filename = {'filename of Matlab file containing forcing data'};
%             
%             forcing.PARA.default_value.forcing_path = {'../CryoGridCommunity_forcing/'};
%             forcing.PARA.comment.forcing_path = {'path where forcing data file is located'};
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
%             forcing.PARA.comment.rain_fraction = {'rainfall fraction assumed in sumulations (rainfall from the forcing data file is multiplied by this parameter)'};
%             
%             forcing.PARA.default_value.snow_fraction = {1};  
%             forcing.PARA.comment.snow_fraction = {'snowfall fraction assumed in sumulations (rainfall from the forcing data file is multiplied by this parameter)'};
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
        end
    end
end