%========================================================================
% CryoGrid FORCING processing class interpolate_sl_ERA2
% uses interpolation of Sin relying on S_TOA
%
% Authors:
% S. Westermann, December 2022
%
%========================================================================

classdef interpolate_sl_MET < process_BASE
    
    
    methods
        function proc = provide_PARA(proc)
            
        end
        
        
        function proc = provide_CONST(proc)
            proc.CONST.Tmfw = [];
        end
        
        
        function proc = provide_STATVAR(proc)
            
        end
        
        
        function proc = finalize_init(proc, tile)
            
        end
        
        
        function forcing = process(proc, forcing, tile)
            
            disp('interpolating surface level data')
            era = forcing.TEMP.era;
            
            [dist,~] = distance(forcing.SPATIAL.STATVAR.latitude,forcing.SPATIAL.STATVAR.longitude,era.latitude,era.longitude);
            
            [dist, ind] = sort(dist(:));
            
            dist=dist(1:4);
            ind = ind(1:4);
            [row, col] = ind2sub(size(era.latitude), ind);
            
            weights = 1./dist ./ sum(1./dist);
            
            altitude  = era.altitude;
            altitude = altitude(row(1), col(1)) .* weights(1) + altitude(row(2), col(2)) .* weights(2) + altitude(row(3), col(3)) .* weights(3) + altitude(row(4), col(4)) .* weights(4);
            
            altitude_4cells = diag(era.altitude(row, col));
            
            Tair=[];
            Tair2=[];
            lr=[];
            for i=1:size(era.T2,3)
                T  = double(era.T2(:,:,i)) .* era.T_sf;
                T_4cells = diag(T(row, col));
                T = T(row(1), col(1)) .* weights(1) + T(row(2), col(2)) .* weights(2) + T(row(3), col(3)) .* weights(3) + T(row(4), col(4)) .* weights(4);
                %     lapse_rate = mean(T_4cells (1:2) - T_4cells (3:4)) ./ mean(altitude_4cells(1:2) - altitude_4cells(3:4));
                [~,ind_max] = max(altitude_4cells);
                [~,ind_min] = min(altitude_4cells);
                lapse_rate = (T_4cells(ind_max) - T_4cells(ind_min)) ./ (altitude_4cells(ind_max) - altitude_4cells(ind_min));
                lr = [lr; lapse_rate];
                Tair = [Tair; T + lapse_rate .* (forcing.SPATIAL.STATVAR.altitude - altitude)];
                Tair2 = [Tair2; T];
            end
            
            
            wind = double(era.wind) .* era.wind_sf;
            wind = squeeze(wind(row(1), col(1), :) .* weights(1) + wind(row(2), col(2),:) .* weights(2) + wind(row(3), col(3),:) .* weights(3) + wind(row(4), col(4), :) .* weights(4));
            
            Sin = double(era.SW) .* era.rad_sf;
            Sin = squeeze(Sin(row(1), col(1), :) .* weights(1) + Sin(row(2), col(2),:) .* weights(2) + Sin(row(3), col(3),:) .* weights(3) + Sin(row(4), col(4), :) .* weights(4));
            
            Lin = double(era.LW) .* era.rad_sf;
            Lin = squeeze(Lin(row(1), col(1), :) .* weights(1) + Lin(row(2), col(2),:) .* weights(2) + Lin(row(3), col(3),:) .* weights(3) + Lin(row(4), col(4), :) .* weights(4));
            
            RH = double(era.RH) .* era.RH_sf;
            RH = squeeze(RH(row(1), col(1), :) .* weights(1) + RH(row(2), col(2),:) .* weights(2) + RH(row(3), col(3),:) .* weights(3) + RH(row(4), col(4), :) .* weights(4));
            
            precip = double(era.P) .* era.P_sf;
            precip = squeeze(precip(row(1), col(1), :) .* weights(1) + precip(row(2), col(2),:) .* weights(2) + precip(row(3), col(3),:) .* weights(3) + precip(row(4), col(4), :) .* weights(4));
            
            pressure = double(era.ps) .* era.ps_sf;
            pressure = squeeze(pressure(row(1), col(1), :) .* weights(1) + pressure(row(2), col(2),:) .* weights(2) + pressure(row(3), col(3),:) .* weights(3) + pressure(row(4), col(4), :) .* weights(4));
            pressure = pressure .* (1 - 0.0065./288.15.*forcing.SPATIAL.STATVAR.altitude).^5.255;
            
            %             [~, sunElevation] = solargeom(proc, era.t, forcing.SPATIAL.STATVAR.latitude, forcing.SPATIAL.STATVAR.longitude);
            %             sunElevation = 90-rad2deg(sunElevation);
            %             mu0=max(sind(sunElevation),0); % Trunacte negative values.
            %             S_TOA = 1370.*mu0';
            
            forcing.DATA.Tair = Tair2;
            forcing.DATA.RH = RH;
            forcing.DATA.wind = wind;
            forcing.DATA.Sin =  Sin;
            forcing.DATA.Lin = Lin;
            forcing.DATA.p = pressure;
            forcing.DATA.precip = precip;  %mm/hour to mm/day
            forcing.DATA.timeForcing = era.t;
            
            forcing.TEMP.era = [];
            
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

