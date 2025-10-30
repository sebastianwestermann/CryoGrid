%========================================================================
% CryoGrid DATA_PROVIDER class CALCULATE_SKYVIEW_FACTOR
% S. Westermann, Oct 2025
%========================================================================

classdef CALCULATE_SKYVIEW_FACTOR < matlab.mixin.Copyable

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

            dem.PARENT.STATVAR.skyview_factor = [];

            for i=1:size(dem.PARENT.STATVAR.key,1)

                horizon_bins = [0:360]';
                horizon_angles = interp1(dem.PARENT.STATVAR.horizon_bins(i,:)', dem.PARENT.STATVAR.horizon_angles(i,:)', horizon_bins);

                azmRadian = (pi/180).*horizon_bins;

                % convert output from horizon program to radians and translate to angle
                % from zenith
                H = (pi/180).*(90 - horizon_angles);

                aspectRadian = (pi/180)*(dem.PARENT.STATVAR.aspect(i,1));
                % modify limits of integration for slopes facing away from horizons
                t = cosd(dem.PARENT.STATVAR.aspect(i,1) - horizon_bins)<0;
                %Simplified trig, the original was H(t) = min(H(t),...
                %  acos(-cos(azmRadian(t)-aspectRadian)*sind(slopeDegrees)./...
                %  sqrt(cosd(slopeDegrees)^2+sind(slopeDegrees)^2*cos(azmRadian(t)-aspectRadian).^2)));
                % but same as
                H(t) = min(H(t), acos(sqrt(1-1./(1+tand(dem.PARENT.STATVAR.slope_angle(i,1))^2*cos(azmRadian(t)-aspectRadian).^2))));
                qIntegrand = (cosd(dem.PARENT.STATVAR.slope_angle(i,1))*sin(H).^2 + sind(dem.PARENT.STATVAR.slope_angle(i,1))*cos(aspectRadian-azmRadian).*(H-cos(H).*sin(H)))/2;

                % shouldn't be any negative, except perhaps rounding error, so just in case
                qIntegrand(qIntegrand<0) = 0;

                % integrate
                dem.PARENT.STATVAR.skyview_factor = [dem.PARENT.STATVAR.skyview_factor; trapz(azmRadian,qIntegrand)./pi];
            end
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

