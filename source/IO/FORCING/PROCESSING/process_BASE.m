%========================================================================
% CryoGrid FORCING processing class
%
%
% Authors:
% S. Westermann, December 2022
%
%========================================================================

classdef process_BASE < matlab.mixin.Copyable 
    
    properties
        CONST
        PARA
        STATVAR
        TEMP
    end
    
    methods
        function proc = provide_PARA(proc)

        end
        
        
        function proc = provide_CONST(proc)

            
        end
        
        
        function proc = provide_STATVAR(proc)
            
        end
        
        
        function proc = finalize_init(proc, tile)

        end
        
        
        function forcing = process(proc, forcing, tile)
            

        end
        
        function p = satPresIce(proc, T)
            Tmfw = proc.CONST.Tmfw;
            p = 6.112.* 100.* exp(22.46.*(T-Tmfw)./(272.61+T-Tmfw));
        end
        
        function p = satPresWater(proc, T)
            Tmfw = proc.CONST.Tmfw;
            p = 6.112 .* 100 .* exp(17.62.*(T-Tmfw)./(243.12+T-Tmfw));
        end
        
%         function q = convertRelative2absoluteHumidity(forcing)
%             p=1005*100;
%             
%             q=(double((forcing.DATA.Tair)<=0).*(forcing.DATA.RH.*satPresIce(forcing, forcing.DATA.Tair+273.15)) + double(forcing.DATA.Tair>0).*(forcing.DATA.RH.*satPresWater(forcing, forcing.DATA.Tair+273.15)))./p;
%         end

        
        %Kris script, used now to be consistent with TopoScale
        function [solar_azimuth,solar_zenith]=solargeom(proc, Time, Latitude,Longitude)  %added two arguments, SW
      %              function [solar_azimuth,solar_zenith, HrAngleR, DecR]=solargeom(forcing, Time, Latitude,Longitude)  %added two arguments, SW
            % [saz,szen]=solargeom(time,latitude,longitude)
            % Adopted from the Sandia National Labs PVL Toolbox ephemeris routine.
            % Inputs:
            %   time = Time stamp vector (matlab datenum format) assumed to be in UTC
            %   latitude = Latitude
            %   longitude = Longitude
            % Outputs:
            %   saz = Solar azimuth angle [radians, anticlockwise from south]
            %   szen = Solar zentih angle [radians].
            % Link to the original toolbox:
            % https://pvpmc.sandia.gov/applications/pv_lib-toolbox/
            % References:
            % Stein et al. (2012), doi:10.1109/PVSC.2012.6318225 [MATLAB version]
            % Holmgren et al. (2018), doi:10.21105/joss.00884 [Python version]
            
            
            TZone=0;
            Longitude=-Longitude;
            % tv=datevec(Time);
            Year=year(Time);
            % v0=zeros(size(Year)); v1=ones(size(Year));
            DayOfYear=floor(Time-datenum(Year,1, 1))+1;
            DecHours=(Time - floor(Time)) .* 24;
            RadtoDeg=180/pi;
            DegtoRad=pi/180;
            Abber = 20/3600;
            LatR = Latitude * DegtoRad;
            UnivDate = DayOfYear + floor((DecHours + TZone)/24);
            UnivHr = mod((DecHours + TZone), 24);
            Yr = Year-1900;
            YrBegin = 365 * Yr + floor((Yr-1)/4)-0.5;
            Ezero = YrBegin + UnivDate;
            T = Ezero / 36525;
            GMST0 = 6/24 +38/1440 + (45.836 + 8640184.542 * T + 0.0929 * T.^2)/86400;
            GMST0 = 360 * (GMST0 - floor(GMST0));
            GMSTi = mod(GMST0 + 360*(1.0027379093 * UnivHr / 24),360);
            LocAST = mod((360 + GMSTi - Longitude), 360);
            EpochDate = Ezero + UnivHr / 24;
            T1 = EpochDate / 36525;
            ObliquityR = DegtoRad * (23.452294 - 0.0130125 * T1 - 0.00000164 * T1.^2 ...
                + 0.000000503 * T1.^3);
            MlPerigee = 281.22083 + 0.0000470684 * EpochDate + 0.000453 * T1 .^ 2 + ...
                0.000003 * T1 .^ 3;
            MeanAnom = mod((358.47583 + 0.985600267 * EpochDate - 0.00015 * T1 .^ 2 - ...
                0.000003 * T1 .^ 3), 360);
            Eccen = 0.01675104 - 0.0000418 * T1 - 0.000000126 * T1 .^ 2;
            EccenAnom = MeanAnom;
            E=0;
            while max(abs(EccenAnom - E)) > 0.0001
                E = EccenAnom;
                EccenAnom = MeanAnom + RadtoDeg .* Eccen .* sin(DegtoRad .* E);
            end
            TrueAnom = 2 * mod(RadtoDeg * atan2(((1 + Eccen) ./ (1 - Eccen)).^ 0.5 .* tan(DegtoRad * EccenAnom / 2), 1), 360) ;
            EcLon = mod(MlPerigee + TrueAnom, 360) - Abber ;
            EcLonR = DegtoRad * EcLon;
            DecR = asin(sin(ObliquityR) .* sin(EcLonR));
            %Dec = RadtoDeg * DecR;
            RtAscen = RadtoDeg * atan2(cos(ObliquityR).*(sin(EcLonR)),cos(EcLonR));
            HrAngle = LocAST - RtAscen ;
            HrAngleR = DegtoRad .* HrAngle ;
            %HrAngle = HrAngle - (360 .* sign(HrAngle) .* (abs(HrAngle) > 180));
            SunAz = RadtoDeg .* atan2(-1 * sin(HrAngleR), cos(LatR) .* tan(DecR) - sin(LatR) .* cos(HrAngleR));
            SunAz = SunAz + (SunAz < 0) * 360; %shift from range of [-180,180] to [0,360]
            SunEl = asind(cos(LatR) .* cos(DecR) .* cos(HrAngleR) + sin(LatR) .* sin(DecR));
            
            % Convert solar azimuth angle from [N,E,S,W]=[0,90,180,270] to [180, 90, 0
            % -90], i.e. the same as the aspect and horizon angle system.
            solar_azimuth=deg2rad(SunAz);
            solar_azimuth=(5*pi/2)-solar_azimuth;
            solar_azimuth=solar_azimuth-2.*pi.*(solar_azimuth>2.*pi);
            solar_azimuth=solar_azimuth+pi/2;
            solar_azimuth=solar_azimuth-2.*pi.*(solar_azimuth>pi);
            
            % Calculate solar zenith angle from solar elevation angle
            SunEl=deg2rad(SunEl);
            solar_zenith=(pi/2)-SunEl;
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

