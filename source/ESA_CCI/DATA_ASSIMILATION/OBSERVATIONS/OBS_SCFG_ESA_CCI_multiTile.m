%made for snow extent DA in ESA CCI, runs for multiple grid cells simultaneously 

%good to the point that PBS produces weights , best-fitting ensemble member
%can be identified and written in outoput file

%must read the correct 300 grid cells from TILE!!!

classdef OBS_SCFG_ESA_CCI_multiTile < matlab.mixin.Copyable

    
    properties
        PARA
        CONST
        STATVAR
        TEMP
    end
    
    methods
        function obs  = provide_PARA(obs )
            
            obs.PARA.deg_tile_list_file = [];
            obs.PARA.deg_tile_list_folder = [];

            obs.PARA.obs_filename = [];
            obs.PARA.obs_folder = [];
            obs.PARA.std_observations = [];
            obs.PARA.start_year = [];
            obs.PARA.cutoff_date = [];
        end
        
        function obs = provide_CONST(obs)
            
        end
        
        function obs = provide_STATVAR(obs)
            
        end
        
        function obs = finalize_init(obs, tile)
            while ~exist('MODIS_deg_list')
                load([obs.PARA.deg_tile_list_folder obs.PARA.deg_tile_list_file], 'MODIS_deg_list');
            end
            obs.TEMP.MODIS_deg_list = MODIS_deg_list;

            obs.TEMP.current_year = obs.PARA.start_year;
        end


        function obs = read_observations(obs, tile)
            obs.STATVAR.observations = [];
            obs.STATVAR.time = [];
            obs.STATVAR.obs_variance = [];

            for i=1:size(obs.TEMP.MODIS_deg_list,1)

                range_in_tile = find(tile.PARA.latitude(:,1) > obs.TEMP.MODIS_deg_list(i,1) & ...
                    tile.PARA.latitude(:,1) < obs.TEMP.MODIS_deg_list(i,2) & tile.PARA.longitude(:,1) > obs.TEMP.MODIS_deg_list(i,3) & ...
                    tile.PARA.longitude(:,1) < obs.TEMP.MODIS_deg_list(i,4));
                dec = 0;
                SCFG_file = [obs.PARA.obs_folder obs.PARA.obs_filename '_' num2str(obs.TEMP.MODIS_deg_list(i,1)) '_' num2str(obs.TEMP.MODIS_deg_list(i,2)) '_' ...
                    num2str(obs.TEMP.MODIS_deg_list(i,3)) '_' num2str(obs.TEMP.MODIS_deg_list(i,4))];
                fname = [SCFG_file '_' num2str(obs.TEMP.current_year) '.nc'];

                if ~isempty(range_in_tile) && exist(fname)==2

                    MODIS_lat = ncread(fname, 'latitude', 1, Inf);
                    MODIS_lon = ncread(fname, 'longitude', 1, Inf);
                    MODIS_delta_lat = ncread(fname, 'delta_lat', 1, Inf);
                    MODIS_delta_lon = ncread(fname, 'delta_lon', 1, Inf);

                    for j=1:size(range_in_tile,1)
                        
                        %find closest pixel in MODIS file
                        [dist, closest_pixel_id] = min((tile.PARA.longitude(range_in_tile(j),1) - MODIS_lon).^2 + (tile.PARA.latitude(range_in_tile(j),1) - MODIS_lat).^2);
                        if dist < 2.*((MODIS_delta_lat).^2 + (MODIS_delta_lon).^2) %pixel actually available in MODIS file
                            
                            observations = double(ncread(fname, 'scfg',  [closest_pixel_id 1], [1 Inf]));

                            if ~dec
                                %make variable
                                obs.STATVAR.observations = zeros(size(tile.PARA.range,1), size(observations,2)) .*NaN;
                                dec = 1;
                            end
                            obs.STATVAR.observations(range_in_tile(j),:) = observations;
                        end
                    end
                end
            end
            if ~isempty(obs.STATVAR.observations)
                obs.STATVAR.observations = obs.STATVAR.observations';
                obs.STATVAR.observations(obs.STATVAR.observations>100) = NaN;
                obs.STATVAR.observations = obs.STATVAR.observations ./ 100;
                obs.STATVAR.time = [datenum(year(tile.t),9,1)+0.5:datenum(year(tile.t),9,1)+size(obs.STATVAR.observations,1)-0.5]';
                if ~isempty(obs.PARA.cutoff_date)
                    cutoff_date = datenum([obs.PARA.cutoff_date num2str(year(tile.t))], 'dd.mm.yyyy');
                    if cutoff_date < tile.t
                        cutoff_date = datenum([obs.PARA.cutoff_date num2str(year(tile.t)+1)], 'dd.mm.yyyy');
                    end
                    delete_obs = find(obs.STATVAR.time(:,1)<cutoff_date);
                    obs.STATVAR.time(delete_obs,:) = [];
                    obs.STATVAR.observations(delete_obs,:) = [];
                end
                obs.STATVAR.obs_variance = zeros(size(obs.STATVAR.observations,1), size(obs.STATVAR.observations,2)) + obs.PARA.std_observations.^2;
            end

            obs.TEMP.current_year = obs.TEMP.current_year + 1;
        end
            
    end
end

