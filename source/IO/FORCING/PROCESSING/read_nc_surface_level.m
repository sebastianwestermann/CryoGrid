%========================================================================
% CryoGrid FORCING processing class
%
%
% Authors:
% S. Westermann, December 2022
%
%========================================================================

classdef read_nc_surface_level < matlab.mixin.Copyable 
    
    properties
        CONST
        PARA
        STATVAR
        TEMP
    end
    
    methods
        function proc = provide_PARA(proc)
            proc.PARA.nc_folder = [];  
            proc.PARA.forcing_path = [];
            proc.PARA.time_resolution_input = [];
        end
        
        
        function proc = provide_CONST(proc)

        end
        
        
        function proc = provide_STATVAR(proc)
            
        end
        
        
        function proc = finalize_init(proc, tile)

        end
        
        
        function forcing = process(proc, forcing, tile)
            
            disp('load nc-files')
            
            if strcmp(proc.PARA.time_resolution_input, 'quarter')
                quarter = forcing.TEMP.current_quarter;
            elseif strcmp(proc.PARA.time_resolution_input, 'month')
                month = forcing.TEMP.current_month;
            end
            year = forcing.TEMP.current_year;
                        
            if strcmp(proc.PARA.time_resolution_input, 'quarter')
                ncf=[proc.PARA.forcing_path proc.PARA.nc_folder 'surf_q' num2str(quarter) '_' num2str(year) '.nc'];
            elseif strcmp(proc.PARA.time_resolution_input, 'month')
                ncf=[proc.PARA.forcing_path proc.PARA.nc_folder 'surf_m' num2str(month) '_' num2str(year) '.nc'];
            elseif strcmp(proc.PARA.time_resolution_input, 'year')
                ncf=[proc.PARA.forcing_path proc.PARA.nc_folder 'surf_' num2str(year) '.nc'];                
            end
                        
            t = ncread(ncf,'time');
            target_timestep = (double(t(2))-double(t(1)))./24;
            
            if strcmp(proc.PARA.time_resolution_input, 'quarter')
                era.t=datenum(year,(quarter-1).*3+1, 1):target_timestep:datenum(year, quarter.*3+1,1)-target_timestep; % Target time steps
            elseif strcmp(proc.PARA.time_resolution_input, 'month')
                era.t=datenum(year,month, 1):target_timestep:datenum(year, month+1,1)-target_timestep; % Target time steps
            elseif strcmp(proc.PARA.time_resolution_input, 'year')
                era.t=datenum(year,1, 1):target_timestep:datenum(year+1,1,1)-target_timestep; % Target time steps
            end

            %only needed during initialization
            era.lon=ncread(ncf,'longitude');
            era.lat=ncread(ncf,'latitude');
            
            % [era.P,era.SW,era.LW,era.ps,era.T2,era.Td2,era.u10,era.v10]=...
            %     deal(nan(numel(era.lon),numel(era.lat),numel(era.t),'single'));
            
            ncfz=[proc.PARA.forcing_path proc.PARA.nc_folder 'gp_surf.nc'];
            tmp=ncread(ncfz,'z');
            era.Zs=tmp(:,:,1)./9.81;

            u10=ncread(ncf,'u10');
            era.u10=u10;
            v10=ncread(ncf,'v10');
            era.v10=v10;
            Td2=ncread(ncf,'d2m');
            era.Td2=Td2;
            T2=ncread(ncf,'t2m');
            era.T2=T2;
            ps=ncread(ncf,'sp');
            era.ps=ps;
            SW=ncread(ncf,'ssrd');
            era.SW=SW./(60.^2);
            LW=ncread(ncf,'strd');
            era.LW=LW./(60.^2);
            S_TOA = ncread(ncf, 'tisr');
            era.S_TOA = S_TOA./(60.^2);
            P=ncread(ncf,'tp');
            era.P=P.*1e3; % m/hour to mm/hour

            % Scale output.
            era.wind_sf=1e-2;
            era.q_sf=1e-6;
            era.ps_sf=1e2;
            era.rad_sf=1e-1;
            era.T_sf=1e-2;
            era.P_sf=1e-2;
            
            era.u10=int16(era.u10./era.wind_sf);
            era.v10=int16(era.v10./era.wind_sf);
            era.ps=uint16(era.ps./era.ps_sf);
            era.SW=uint16(era.SW./era.rad_sf);
            era.LW=uint16(era.LW./era.rad_sf); % Do this in the time loop.
            era.S_TOA=uint16(era.S_TOA./era.rad_sf);
            era.T2=int16((era.T2-273.15)./era.T_sf);
            era.Td2=int16((era.Td2-273.15)./era.T_sf);
            era.P=uint16(era.P./era.P_sf);
%             era.Z=int16(era.Z);
            
            forcing.TEMP.era = era;

            
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

