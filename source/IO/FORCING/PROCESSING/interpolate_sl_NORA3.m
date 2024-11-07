%========================================================================
% CryoGrid FORCING processing class interpolate_sl_NORA3
% uses interpolation of Sin relying on S_TOA
%
% Authors:
% S. Westermann, December 2022
%
%========================================================================

classdef interpolate_sl_NORA3 < process_BASE
    
    
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
            
            wind = sqrt((double(era.u10) .* era.wind_sf).^2 + (double(era.v10) .* era.wind_sf).^2);
            wind = squeeze(wind(row(1), col(1), :) .* weights(1) + wind(row(2), col(2),:) .* weights(2) + wind(row(3), col(3),:) .* weights(3) + wind(row(4), col(4), :) .* weights(4));
            
            Sin = double(era.SW) .* era.rad_sf;
            Sin = squeeze(Sin(row(1), col(1), :) .* weights(1) + Sin(row(2), col(2),:) .* weights(2) + Sin(row(3), col(3),:) .* weights(3) + Sin(row(4), col(4), :) .* weights(4));
            
            Lin = double(era.LW) .* era.rad_sf;
            Lin = squeeze(Lin(row(1), col(1), :) .* weights(1) + Lin(row(2), col(2),:) .* weights(2) + Lin(row(3), col(3),:) .* weights(3) + Lin(row(4), col(4), :) .* weights(4));
            
            RH = double(era.RH) .* era.T_sf;
            RH = squeeze(RH(row(1), col(1), :) .* weights(1) + RH(row(2), col(2),:) .* weights(2) + RH(row(3), col(3),:) .* weights(3) + RH(row(4), col(4), :) .* weights(4));
            
            precip = double(era.P) .* era.P_sf;
            precip = squeeze(precip(row(1), col(1), :) .* weights(1) + precip(row(2), col(2),:) .* weights(2) + precip(row(3), col(3),:) .* weights(3) + precip(row(4), col(4), :) .* weights(4));
            
            pressure = double(era.ps) .* era.ps_sf;
            pressure = squeeze(pressure(row(1), col(1), :) .* weights(1) + pressure(row(2), col(2),:) .* weights(2) + pressure(row(3), col(3),:) .* weights(3) + pressure(row(4), col(4), :) .* weights(4));

            Tair = double(era.T2) .* era.T_sf;
            Tair = squeeze(Tair(row(1), col(1), :) .* weights(1) + Tair(row(2), col(2),:) .* weights(2) + Tair(row(3), col(3),:) .* weights(3) + Tair(row(4), col(4), :) .* weights(4));

            
            forcing.DATA.Tair = Tair;
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

