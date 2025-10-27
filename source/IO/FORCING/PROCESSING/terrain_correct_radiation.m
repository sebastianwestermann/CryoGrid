%========================================================================
% CryoGrid FORCING processing class
%
%
% Authors:
% S. Westermann, December 2022
%
%========================================================================

classdef terrain_correct_radiation < process_BASE
    
    
    methods
        function proc = provide_PARA(proc)
            proc.PARA.albedo_surrounding_terrain = [];
        end
        
        
        function proc = provide_CONST(proc)
            forcing.CONST.sigma = []; %Stefan-Boltzman constant
            forcing.CONST.Tmfw = [];
        end
        
        
        function proc = provide_STATVAR(proc)
            
        end
        
        
        function proc = finalize_init(proc, tile)

        end
        
        
        function forcing = process(proc, forcing, tile)
            
            forcing.PARA.is_slope = 1;
                        
            forcing = SolarAzEl(proc, forcing, tile);
            forcing = check_S_TOA(proc, forcing, tile);
            if ~isfield(forcing.DATA, 'Sin_dir')
                forcing = split_Sin(proc, forcing, tile); % split Sin in dir and dif
            end
            forcing = terrain_corr_Sin_dif(proc, forcing, tile);
            forcing = reproject_Sin_dir(proc, forcing, tile);
            forcing = terrain_shade(proc, forcing, tile);
            forcing.DATA.Sin = forcing.DATA.Sin_dir + forcing.DATA.Sin_dif;
            
        end
        
        function forcing = check_S_TOA(proc, forcing, tile)
            if ~isfield(forcing.DATA, 'S_TOA')
                mu0=max(sind(forcing.DATA.sunElevation),0); % Trunacte negative values.
                sunset=mu0<cosd(89);%(mu0==0); % Sunset switch.
                forcing.DATA.S_TOA = 1370.*mu0;
            end
        end
        

        
        function forcing = split_Sin(proc, forcing, tile) 
            % Split Sin into direct and diffuse parts, see Fiddes(2014)
            Sin = forcing.DATA.Sin;
            Sin_TOA = forcing.DATA.S_TOA;
            
            kt = Sin ./ Sin_TOA; %clearness index
            kt(isnan(kt))=0;
            kd = 0.952 - 1.041.*exp(-exp(2.300 - 4.702 .* kt)); %diffuse fraction
            kd(isnan(kd))=0;
            kd = max(0,kd);
            
            forcing.DATA.Sin_dif = kd .* forcing.DATA.Sin; % full hemisphere diffuse Sin
            forcing.DATA.Sin_dir = (1-kd) .* forcing.DATA.Sin; % direct Sin on horizontal surface
            
        end
        
        function forcing = terrain_corr_Sin_dif(proc, forcing, tile)
            % include effect of terrain on diffuse by removin the fraction
            % of the hemisphere covered by the horizon, and adding
            % reflected Sin from surrounding terrain
            
            Sin = forcing.DATA.Sin; % Total Sin (horizontal)
            Sin_dif = forcing.DATA.Sin_dif;% Diffuse Sin (horizontal)
            alpha = proc.PARA.albedo_surrounding_terrain; %Albedo at the foot of the slope
            svf = forcing.SPATIAL.STATVAR.skyview_factor; % hemispheric fraction of sky not occluded by terrain
            
            forcing.DATA.Sin_dif = Sin_dif.*svf + Sin.*alpha.*(1-svf);
        end
        
        function forcing = reproject_Sin_dir(proc,forcing, tile)
            Sin_dir = forcing.DATA.Sin_dir;
            aspect = forcing.SPATIAL.STATVAR.aspect;
            slope = forcing.SPATIAL.STATVAR.slope_angle;
            
            surf_norm_vec = repmat([0.0,0.0,1.0], size(forcing.DATA.Sin,1), 1); %Unit vector on the horizontal
            
            alpha = aspect.*pi./180; %Degree to radians of exposition of the slope
            beta = (90-slope).*pi./180; %Degree to radians of inclination of the slope
            face_vec = repmat([sin(alpha).*cos(beta) cos(alpha).*cos(beta) sin(beta)], size(forcing.DATA.Sin,1), 1); %Unit vector of the slope
            
            % Calculation the solar azimuth and elevation angle relative to
            % the coordinates of the site (revised after Darin C. Koblick)
            % ->CHANGED to Kris script
%             forcing = SolarAzEl(forcing, tile);
            
            alpha = forcing.DATA.azimuth.*pi/180; %Degree to radians of the azimuth
            beta = forcing.DATA.sunElevation.*pi/180; %Degree to radians of elevation
            sun_vec = [sin(alpha).*cos(beta) cos(alpha).*cos(beta) sin(beta)]; %Unit vector of the radiation
            
            delta_angle_surf_norm = acos(dot(surf_norm_vec' ,sun_vec')').*180./pi; %angle between the radiation and the normal on the horizontal in degrees
            delta_angle_face=acos(dot(face_vec', sun_vec')').*180./pi; %angle between the radiation and the normal on the slope in degrees
            
            Sin_sun_direction = Sin_dir./cos(delta_angle_surf_norm.*pi./180); %direct Sin in direction of the sun
            Sin_face_direction = double(delta_angle_surf_norm < 89.5 & delta_angle_face < 89.5) .* Sin_sun_direction.*cos(delta_angle_face.*pi/180); % Direct Sin in direction of surface face
            
            forcing.DATA.Sin_dir = Sin_face_direction;
        end
        
        function forcing = terrain_corr_Lin(proc, forcing, tile)
            % Reduce Lin from atmosphere to sky view fraction, add Lin from
            % surroundings assuming Tair everywhere
            sigma = proc.CONST.sigma; %Stefan-Boltzman constant
            Lin = forcing.DATA.Lin;
            Tair = forcing.DATA.Tair;
            svf = forcing.SPATIAL.STATVAR.skyview_factor;
            Tmfw = proc.CONST.Tmfw;
            
            forcing.DATA.Lin = svf.*Lin + (1-svf).*sigma.*(Tair+Tmfw).^4;
        end
        
        function forcing = terrain_shade(proc, forcing, tile)
            az = forcing.DATA.azimuth;
            el = forcing.DATA.sunElevation;
            hbins = forcing.SPATIAL.STATVAR.horizon_bins;
            h = forcing.SPATIAL.STATVAR.horizon_angles;
            
            I = knnsearch(hbins,az); % hbin containing current solar azimuth -> INTERPOLATE INSTEAD??? in the DEM analysis, this is points along lines, not bins!!
            forcing.DATA.Sin_dir(h(I)>el) = 0; % remove direct Sin if sun is below horizon

        end
        
          function forcing = SolarAzEl(proc, forcing, tile)
           [forcing.DATA.azimuth,forcing.DATA.sunElevation] = solargeom(proc, forcing.DATA.timeForcing, forcing.SPATIAL.STATVAR.latitude, forcing.SPATIAL.STATVAR.longitude);
           %-----
          % [forcing.DATA.azimuth,forcing.DATA.sunElevation, forcing.DATA.hour_angle, forcing.DATA.DecR] = solargeom(forcing, forcing.DATA.timeForcing ,tile.PARA.latitude,tile.PARA.longitude); %CHANGE SW
%            forcing.DATA.hour_angle = rad2deg(forcing.DATA.hour_angle);
%            forcing.DATA.hour_angle(forcing.DATA.hour_angle<0) = forcing.DATA.hour_angle(forcing.DATA.hour_angle<0)+360;
%            forcing.DATA.hour_angle(forcing.DATA.hour_angle>=360) = forcing.DATA.hour_angle(forcing.DATA.hour_angle>=360)-360;
           %----------------

           forcing.DATA.azimuth = rad2deg(forcing.DATA.azimuth);
           forcing.DATA.azimuth(forcing.DATA.azimuth<0) = forcing.DATA.azimuth(forcing.DATA.azimuth<0) + 360;
           forcing.DATA.sunElevation = 90-rad2deg(forcing.DATA.sunElevation);
        end
        
        
                %-------------param file generation-----
%         function post_proc = param_file_info(post_proc)
%             post_proc = provide_PARA(post_proc);
% 
%             post_proc.PARA.STATVAR = [];
%             post_proc.PARA.class_category = 'FORCING POST_PROCESSING';
%             post_proc.PARA.options = [];
%             
%             post_proc.PARA.eliminate_fraction = [];
%             post_proc.PARA.survive_fraction = [];
%                         
%             post_proc.PARA.default_value.window_size = {7};
%             post_proc.PARA.comment.window_size = {'window size in days within which precipitation is reallocated'};
%             
%             post_proc.PARA.default_value.eliminate_fraction = {0.5};
%             post_proc.PARA.comment.eliminate_fraction = {'fraction of smallest precipitation events (= timestamps with precipitation) that is reallocated to larger events'};
%             
%             post_proc.PARA.default_value.survive_fraction = {0.5};  
%             post_proc.PARA.comment.survive_fraction = {'fraction of largest precipitation events (= timestamps with precipitation) that the small events are reallocated to'};
%             
%         end
        
    end
    
end

