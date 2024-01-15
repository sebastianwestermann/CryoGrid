%========================================================================
% CryoGrid DATA_PROVIDER class DEM
% DEM class deriving information from a digital 
% elevation model (DEM),including slope, aspect (terrain shading still 
% lacking). The class can ingest Copernicus 30m DEMs downloaded from 
% https://portal.opentopography.org/   
%
% S. Westermann, Dec 2022
%========================================================================

classdef READ_DATASET < matlab.mixin.Copyable

    properties
        PARA
        CONST
        STATVAR
        TEMP
        PARENT
    end
    
    methods
        function dem = provide_PARA(dem)            
            dem.PARA.variable_name = []; % name that the variable will be given
            dem.PARA.data_folder = [];
            dem.PARA.data_filename = [];
            
            dem.PARA.reproject2utm = [];

        end
        
        function dem = provide_STATVAR(dem)

        end
        
        function dem = provide_CONST(dem)
            
        end
        
        function dem = finalize_init(dem)
            
              
        end
        
        function dem = load_data(dem)
            
            data = readgeoraster([dem.PARA.data_folder dem.PARA.data_filename]);
            
            data = data(dem.PARENT.STATVAR.key);
            dem.PARENT.STATVAR.(dem.PARA.variable_name) = data;
            
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

