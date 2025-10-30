%========================================================================
% CryoGrid DATA_PROVIDER class SET_FLAT_HORIZON
% S. Westermann, Oct 2025
%========================================================================

classdef SET_FLAT_HORIZON < matlab.mixin.Copyable

    properties
        PARA
        CONST
        STATVAR
        TEMP
        PARENT
    end
    
    methods
        function dem = provide_PARA(dem)            
        end
        
        function dem = provide_STATVAR(dem)

        end
        
        function dem = provide_CONST(dem)
            
        end
        
        function dem = finalize_init(dem)
            
              
        end
        
        function dem = load_data(dem)
            
            dem.PARENT.STATVAR.horizon_angles = repmat([0 0], size(dem.PARENT.STATVAR.key,1),1);
            dem.PARENT.STATVAR.horizon_bins = repmat([0 360], size(dem.PARENT.STATVAR.key,1),1);

        end
        

        

        
        %-------------param file generation-----
%         function point = param_file_info(point)
%             point = provide_PARA(point);
%             
%             point.PARA.STATVAR = [];
%             point.PARA.class_category = 'DATA_PROVIDER';
%             
%             point.PARA.comment.variables = {'properties calculated from DEM: altitude OR altitude, slope_angle, aspect'};
%             point.PARA.options.variables.name = 'H_LIST';
%             point.PARA.options.variables.entries_x = {'altitude' 'slope_angle' 'aspect'};
%                         
% %             point.PARA.comment.number_of_horizon_bins = {'number of angular points for which horizon is calculated; must be multiple of 4'};
% %             point.PARA.default_value.number_of_horizon_bins = {24};
%             
%             point.PARA.comment.DEM_folder = {'folder in which DEM file is located'}; 
%             
%             point.PARA.comment.DEM_filename = {'name of DEM file'}; 
%             
%             point.PARA.comment.reproject2utm = {'select 1 when using a DEM in geographic coordinates (or similar) and computing more than just altitude; select 0 to speed up altitde computation in big DEMs'};
%             point.PARA.default_value.reproject2utm = {1};
%         end
        
    end
end

