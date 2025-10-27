

classdef DA_OBS_read_mat < matlab.mixin.Copyable

    
    properties
        PARA
        CONST
        STATVAR
        TEMP
    end
    
    methods

        function obs  = provide_PARA(obs)

            obs.PARA.obs_filename = [];
            obs.PARA.obs_folder = [];

        end
        
        function obs = provide_CONST(obs)
            
        end
        
        function obs = provide_STATVAR(obs)
            
        end
        
        function obs = finalize_init(obs, tile)
            obs.TEMP.obs_loaded = 0;
        end

        function obs = read_observations(obs, tile)

            if ~obs.TEMP.obs_loaded
                temp=load([obs.PARA.obs_folder obs.PARA.obs_filename], 'OBS');
                obs.STATVAR.time = temp.OBS.time;
                obs.STATVAR.observations = temp.OBS.observations;
                obs.STATVAR.obs_variance = temp.OBS.obs_variance;

                obs.TEMP.obs_loaded = 1;
                
            end
        end

    end
end

