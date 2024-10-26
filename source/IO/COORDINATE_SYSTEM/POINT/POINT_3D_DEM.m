%========================================================================
% CryoGrid SPATIAL_REFERENCE class POINT_3D_DEM
% 3D POINT class designed to provide information for several 
% coupled TILE classes that are run in parallel with the RUN_INFO class 
% RUN_3D_POINT. The topological relationships between the TILE classes are
% provided in this class
%
% S. Westermann, Dec 2022
%========================================================================

classdef POINT_3D_DEM < POINT_DEM

    
    methods
        function point = provide_PARA(point)
            point.PARA.number_of_tiles = [];     
            point.PARA.latitude = [];
            point.PARA.longitude = [];
            point.PARA.area = [];
            point.PARA.param_file_number = []; %[1;2;3];
            
            point.PARA.connected = [];
            point.PARA.contact_length = [];
            point.PARA.distance = [];
            
            point.PARA.variables = []; % altitude or altitude, slope_angle, aspect or altitude, slope_angle, aspect, horizon_angles
            
            point.PARA.number_of_horizon_bins = []; %multiples of 4!
            point.PARA.DEM_folder = [];
            point.PARA.DEM_filename = [];
            point.PARA.reproject2utm = 1; %select 1 when using a geographic coordinate system and computing more than just altitude; select 0 to speed up altitde computation in big DEMs

        end
        
        function point = provide_STATVAR(point)

        end
        
        function point = provide_CONST(point)
            
        end
        
        function point = finalize_init(point)
            %provide default values
            point.PARA.slope_angle = 0;     
            point.PARA.aspect = 0;     
            point.PARA.skyview_factor = 0;     
            point.PARA.horizon_bins = 0;
            point.PARA.horizon_angles = 0;
            
            point = read_DEM_raster(point);
            
            if size(point.PARA.latitude,1) == 1 && size(point.PARA.longitude,1) == 1
                point.STATVAR.latitude = repmat(point.PARA.latitude, 1, point.PARA.number_of_tiles);
                point.STATVAR.longitude = repmat(point.PARA.longitude, 1, point.PARA.number_of_tiles);
                
                %do the DEM analysis once and assign the values to all
                %tiles
                point = project_target_coordinates(point);
                point = compute_global_offset_from_north(point);
                
                for i=1:size(point.PARA.variables,1)
                    a = str2func(['get_' point.PARA.variables{i,1}]);
                    point = a(point);
                end
                %then duplicate as normal
                point.STATVAR.slope_angle = repmat(point.STATVAR.slope_angle, 1, point.PARA.number_of_tiles);
                point.STATVAR.aspect = repmat(point.STATVAR.aspect, 1, point.PARA.number_of_tiles);
                point.STATVAR.skyview_factor = repmat(point.STATVAR.skyview_factor, 1, point.PARA.number_of_tiles);
                point.STATVAR.horizon_angles = repmat(point.STATVAR.horizon_angles, 1, point.PARA.number_of_tiles);
                point.STATVAR.horizon_bins = repmat(point.STATVAR.horizon_bins, 1, point.PARA.number_of_tiles);
            else
                point.STATVAR.latitude = point.PARA.latitude';
                point.STATVAR.longitude = point.PARA.longitude';
                point.STATVAR.altitude = [];
                point.STATVAR.slope_angle = [];
                point.STATVAR.aspect = [];
                point.STATVAR.skyview_factor = [];
                point.STATVAR.horizon_angles = [];
                point.STATVAR.horizon_bins = [];
                
                %do the DEM analysis for each worker
                for j=1:size(point.STATVAR.latitude,2)
                    point2 = copy(point);
                    point2.PARA.latitude = point.PARA.latitude(j,1);
                    point2.PARA.longitude = point.PARA.longitude(j,1);
                    point2 = project_target_coordinates(point2);
                    point2 = compute_global_offset_from_north(point2);
                    for i=1:size(point2.PARA.variables,1)
                        a = str2func(['get_' point2.PARA.variables{i,1}]);
                        point2 = a(point2);
                    end
                    %assign from poin to point2
                    point.STATVAR.altitude = [point.STATVAR.altitude point2.STATVAR.altitude];
                    point.STATVAR.slope_angle = [point.STATVAR.slope_angle point2.STATVAR.slope_angle];
                    point.STATVAR.aspect = [point.STATVAR.aspect point2.STATVAR.aspect];
                    point.STATVAR.skyview_factor = [point.STATVAR.skyview_factor point2.STATVAR.skyview_factor];
                    point.STATVAR.horizon_angles = [point.STATVAR.horizon_angles point2.STATVAR.horizon_angles];
                    point.STATVAR.horizon_bins = [point.STATVAR.horizon_bins point2.STATVAR.horizon_bins];
                end
            end

            if size(point.PARA.area,1) == 1
                point.STATVAR.area = repmat(point.PARA.area, 1, point.PARA.number_of_tiles);
            else
                point.STATVAR.area = point.PARA.area';
            end
            
        end
        
        
        
        %-------------param file generation-----
        function point = param_file_info(point)
            point = provide_PARA(point);
            
            point.PARA.STATVAR = [];
            point.PARA.class_category = 'SPATIAL_REFERENCE';
            point.PARA.default_value = [];
            
            point.PARA.default_value.number_of_tiles = {3};
            point.PARA.comment.number_of_tiles = {'number of tiles/cores'};
            
            point.PARA.comment.latitude = {'latitude in decimal degrees, either single value assumed for all tiles, or list for all TILEs'};
            
            point.PARA.comment.longitude = {'longitude in decimal degrees, either single value assumed for all tiles, or list for all TILEs'};
            
            point.PARA.comment.altitude = {'altitude in m, either single value assumed for all tiles, or list for all TILEs'};
            point.PARA.options.altitude.name = 'H_LIST';
            point.PARA.options.altitude.entries_x = {31 32.4 33};
            
            point.PARA.comment.area = {'area in m2, either single value assumed for all tiles, or list for all TILEs'};
            point.PARA.default_value.area = {1};
            
            point.PARA.comment.param_file_number = {'index of parameter file used for each TILE'};
            point.PARA.options.param_file_number.name = 'H_LIST';
            point.PARA.options.param_file_number.entries_x = {'1' '2' '3'};

            point.PARA.comment.connected = {'matrix of conectivity between TILES, 1: connected; 0: not connected'};
            point.PARA.options.connected.name = 'MATRIX';
            point.PARA.options.connected.entries_matrix = {'0' '1' '0'; '1' '0' '1'; '0' '1' '0'};

            point.PARA.comment.contact_length = {'contact length in m between TILES'};
            point.PARA.options.contact_length.name = 'MATRIX';
            point.PARA.options.contact_length.entries_matrix = {'0' '1' '0'; '1' '0' '1'; '0' '1' '0'};

            point.PARA.comment.distance = {'distance in m between TILES'};
            point.PARA.options.distance.name = 'MATRIX';
            point.PARA.options.distance.entries_matrix = {'0' '1' '0'; '1' '0' '1'; '0' '1' '0'};
        end
        
    end
end

