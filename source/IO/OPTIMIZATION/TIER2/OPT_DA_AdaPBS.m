classdef OPT_DA_AdaPBS < OPT_DA_FUNCTIONS
    
    properties

    end
    
    methods
        function da = provide_PARA(da)

            da.PARA.run_mode = 'OPT_DA_IO_TILE_sequential'; %or DA_IO_MULTITILE  -> not necessary, this is TILE/MULTITLE in run_info
            da.PARA.scratchfolder = [];
            da.PARA.observation_classes = [];
            da.PARA.observation_classes_index = [];
            da.PARA.observable_classes = [];
            da.PARA.observable_classes_index = []; %must all have the same length, i.e. each observational data set requires one observable class 
            da.PARA.assimilation_frequency = []; %year, month or day
            da.PARA.assimilation_interval = []; %number of years, months or days
            da.PARA.assimilation_date = []; %specific date in case of years or months
            da.PARA.start_assimilation_period = []; %Hlist, date from when the assimilation is started, i.e. the initial state
            da.PARA.ensemble_variables = [];
            da.PARA.ensemble_variable_id = [];
            da.PARA.learning_coefficient = [];
            da.PARA.min_ensemble_diversity = [];
            da.PARA.max_iterations = [];
            da.PARA.store_format = [];
            da.PARA.store_file_tag = [];
            da.PARA.new_init_tile = 0;
            da.PARA.recalculate_stratigraphy = 0;

        end
        
        function da = provide_CONST(da)
            
        end
        
        function da = provide_STATVAR(da)
            
        end
        
        function da = finalize_init(da, run_info)
            da.RUN_INFO = run_info;
            da.TEMP.num_iterations = 0; %was 1 before!!
            da.TEMP.recalculate_stratigraphy_now = 0;
            da.TEMP.recalculate_forcing_now = 0;

            da.ENSEMBLE.modeled_obs = []; %Yp in Kris code, size N_obs x N_ens x N_iterations
            da.ENSEMBLE.value_gaussian = [];  %Xp in Kris code, size N_param x N_ens x N_iterations

            run_mode_class = str2func(da.PARA.run_mode);
            da.IO = run_mode_class();
            %da = make_scratchfolder(da.IO, da, tile);

            %is this necessary here, when everything is stored in ENSEMBLE?
            %->do this later in DA_STEP !

            da = get_ensemble_info(da, run_info); 
            % 
             da.TEMP.old_mean_gaussian = run_info.ENSEMBLE.TEMP.mean_gaussian(da.TEMP.pos_in_ensemble,1);
             da.TEMP.old_std_gaussian = run_info.ENSEMBLE.TEMP.std_gaussian(da.TEMP.pos_in_ensemble,1);
            % %fill the necessary inforation for first AMIS
            da.TEMP.cov_gaussian_resampled = diag(run_info.ENSEMBLE.TEMP.std_gaussian(da.TEMP.pos_in_ensemble,1).^2);
            da.TEMP.mean_gaussian_resampled = run_info.ENSEMBLE.TEMP.mean_gaussian(da.TEMP.pos_in_ensemble,1);

            da = load_observations(da, run_info);

            da.TEMP.ACTIVE = 1;
            da.TEMP.new_init_tile = 1; %1 initialze from scratch; 0 use stored state
            da.TEMP.new_iteration = 0; %1 new iteration needed, back to the start; 0 very first start or go to next DA window
        end

        function [da, new_tile] = get_tile_class(da, run_info)
            if da.TEMP.new_init_tile == 1 
                new_tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_class{1,1}){run_info.PARA.tile_class_index(1,1),1});
                for jj=1:size(da.PARA.observable_classes,1)
                    obs_class = da.PARA.observable_classes{jj,1}; %this is actually OUT classes
                    obs_class_index = da.PARA.observable_classes_index(jj,1);
                    run_info.PPROVIDER.CLASSES.(obs_class){obs_class_index,1}.PARA.timestamps = da.STATVAR.obs_time{jj,1};
                    new_tile = add_OUT_class(new_tile, obs_class, obs_class_index);
                end
                new_tile.RUN_INFO = run_info;
                new_tile = finalize_init(new_tile);
            else
                new_tile = copy(run_info.PPROVIDER.STORAGE.saved_tile{da.TEMP.realization_number,1});
                new_tile.RUN_INFO = run_info; %2x needed!
                new_tile = finalize_init(new_tile);
                new_tile.RUN_INFO = run_info;
            end
        end

        function [da, tile] = add_save_state_classes(da, tile)
            if da.TEMP.new_init_tile == 1
                save_state_class = generate_prior_save_state_class(da, da.TEMP.realization_number, tile.PARA.start_time);
                save_state_class = finalize_init(save_state_class, tile);
                tile = add_OUT_class2(tile, save_state_class);

                save_state_class = generate_final_save_state_class(da, da.TEMP.realization_number);
                save_state_class = finalize_init(save_state_class, tile);
                tile = add_OUT_class2(tile, save_state_class);
            else
                tile = change_output_time_prior(da, tile, da.TEMP.realization_number); %go through existing OUT classes and look for identifier classes
                tile = change_output_time_final(da, tile, da.TEMP.realization_number);
            end
        end

        function [da, tile] = reset_observable_classes(da, tile)
            for jj=1:size(da.PARA.observable_classes,1) %reset modeled observations
                for i=1:size(tile.OUT,1)
                    if strcmp(class(tile.OUT{i,1}), da.PARA.observable_classes{jj,1}) && tile.OUT{i,1}.PARA.class_index == da.PARA.observable_classes_index(jj,1)
                        tile.OUT{i,1} = finalize_init(tile.OUT{i,1}, tile);
                    end
                end
            end
        end

        function da = move_obs2opt(da, tile)
            if tile.t >= da.DA_STEP_TIME
                da = obs2opt(da.IO, da, tile);
            end
        end

        %need something to call before TILE is initialized, this adds the
        %OUT classes to tile (must be known to OPT class) 
                
        function da = DA_step(da, run_info)
            if run_info.TILE.t >= da.DA_STEP_TIME

                da.TEMP.ACTIVE = 1;

                if  da.PARA.new_init_tile == 1
                    da.TEMP.new_init_tile = 1; %1 initialze from scratch; 0 use stored state
                else
                    da.TEMP.new_init_tile = 0;
                end

                da.TEMP.ensemble_size_DA = max(run_info.ENSEMBLE.STATVAR.(da.PARA.ensemble_variable_id));
                %can be different from ensemble_size in ENSEMBLE, since there
                %can be further ENSEMBLE classes, over which the modelled
                %observations integrate

                da.TEMP.num_iterations = da.TEMP.num_iterations + 1; %num_iterations is now the value of the current iteration

                % da = collect_modeled_observations(da.IO, da, run_info);
                % %necessary when parallel

                da = collect_observations(da, run_info.TILE);
                da = collect_modeled_observations(da.IO, da, run_info);

                %---------------
                da = AMIS(da);

                %store DA results, in case user wishes
                da = save_da_results_all(da, run_info);

                %start the iterative loop which either terminates or starts a new iteration, dependent on diversity of ensemmble
                if da.ENSEMBLE.effective_ensemble_size./da.TEMP.ensemble_size_DA >= da.PARA.min_ensemble_diversity || da.TEMP.num_iterations>=da.PARA.max_iterations
                    %terminate and move on in time, i.e. do a normal resampling of model state and resample parameters according to the learning coefficient

                    if da.ENSEMBLE.effective_ensemble_size./da.TEMP.ensemble_size_DA >= da.PARA.min_ensemble_diversity
                        disp('successful, ensemble is sufficiently diverse')
                    else
                        disp('maximum number of iterations reached')
                    end

                    da = save_da_results_final(da, run_info);

                    da.TEMP.new_iteration = 0;
                    da.TEMP.new_init_tile = 0;
                    da.PARA.new_init_tile = 0;

                    da = resample_state_and_parameters(da, run_info);

                    da.TEMP.num_iterations = 0; %reset number of iterations for the next DA period

                    %reset modelled observations and vectors conatining parameters
                    da.ENSEMBLE.modeled_obs = []; %Yp in Kris code, size N_obs x N_ens x N_iterations
                    da.ENSEMBLE.value_gaussian = [];  %Xp in Kris code, size N_param x N_ens x N_iterations
                    da = get_ensemble_info(da, run_info);

                    da.TEMP.old_mean_gaussian = run_info.ENSEMBLE.TEMP.mean_gaussian(da.TEMP.pos_in_ensemble,1);
                    da.TEMP.old_std_gaussian = run_info.ENSEMBLE.TEMP.std_gaussian(da.TEMP.pos_in_ensemble,1);
                    % %fill the necessary inforation for first AMIS
                    da.TEMP.cov_gaussian_resampled = diag(run_info.ENSEMBLE.TEMP.std_gaussian(da.TEMP.pos_in_ensemble,1).^2);
                    da.TEMP.mean_gaussian_resampled = run_info.ENSEMBLE.TEMP.mean_gaussian(da.TEMP.pos_in_ensemble,1);
                else
                    %iterate and go back to start of DA period, resample and inflate again, start over with "old" model states at the beginning of the DA period
                    disp('ensemble degenerate, one more time')

                    da.TEMP.new_iteration = 1;

                    %load "old" state, i.e. from 1st iteration
                    da = assign_prior_state(da, run_info);

                    da = resample_AMIS(da, run_info);

                end
                % if da.PARA.recalculate_stratigraphy == 1
                %     % da.TEMP.recalculate_stratigraphy_now = 1;
                %     da.TEMP.initialize_new = 1;
                % else
                %     da.TEMP.initialize_new = 0;
                % end
            end

        end
    end
end

