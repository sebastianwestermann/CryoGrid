%made for snow extent DA in ESA CCI, runs for multiple grid cells simultaneously 

%good to the point that PBS produces weights , best-fitting ensemble member
%can be identified and written in outoput file

%must read the correct 300 grid cells from TILE!!!

classdef OBS_SCFG_ESA_CCI < matlab.mixin.Copyable

    
    properties
        PARA
        CONST
        STATVAR
        TEMP
    end
    
    methods
        function obs  = provide_PARA(obs )
            
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
            obs.TEMP.current_year = obs.PARA.start_year;
        end

        function obs = read_observations(obs, tile)
            if exist([obs.PARA.obs_folder obs.PARA.obs_filename '_' num2str(obs.TEMP.current_year) '.nc'])==2
                % obs.STATVAR.observations = double(ncread([obs.PARA.obs_folder obs.PARA.obs_filename '_' num2str(obs.TEMP.current_year) '.nc'], 'scfg', [tile.PARA.range(1) 1], [tile.PARA.range(end)-tile.PARA.range(1)+1 Inf], [1 1]));
                obs.STATVAR.observations = double(ncread([obs.PARA.obs_folder obs.PARA.obs_filename '_' num2str(obs.TEMP.current_year) '.nc'], 'scfg'));
                obs.STATVAR.observations = obs.STATVAR.observations(tile.PARA.key_internal,:); %select range of tile, must be key_internal in case of clipping.
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
            else
                obs.STATVAR.observations = [];
                obs.STATVAR.time = [];
                obs.STATVAR.obs_variance = [];
            end
            obs.TEMP.current_year = obs.TEMP.current_year + 1;
        end

    end
end

