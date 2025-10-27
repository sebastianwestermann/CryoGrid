classdef DA_iterative_particle_batch_smoother_AMIS_ML2 < DA_FUNCTIONS
    
    properties
        ML
    end
    
    methods
        function da = provide_PARA(da)

            da.PARA.run_mode = 'DA_IO_TILE'; %or DA_IO_MULTITILE
            da.PARA.scratchfolder = [];
            da.PARA.observation_classes = [];
            da.PARA.observation_classes_index = [];
            da.PARA.observable_classes = [];
            da.PARA.observable_classes_index = []; %must all have the same length, i.e. each observational data set requires one observable class 
            
            da.PARA.ml_class = [];
            da.PARA.ml_class_index = [];
            da.PARA.ml_training_class = [];
            da.PARA.ml_training_class_index = [];
            da.PARA.assimilation_frequency = []; %year, month or day
            da.PARA.assimilation_interval = []; %number of years, months or days
            da.PARA.assimilation_date = []; %specific date in case of years or months
            da.PARA.start_assimilation_period = []; %Hlist, date from when the assimilation is started, i.e. the initial state
            da.PARA.ensemble_variables = [];
            da.PARA.learning_coefficient = [];
            da.PARA.min_ensemble_diversity = [];
            da.PARA.max_iterations = [];
            da.PARA.store_format = [];
            da.PARA.store_file_tag = [];
            da.PARA.recalculate_stratigraphy = 0;

        end
        
        function da = provide_CONST(da)
            
        end
        
        function da = provide_STATVAR(da)
            
        end
        
        function da = finalize_init(da, tile)
            da.TILE = tile;
            da.TEMP.run_name = da.TILE.PARA.run_name; %original results folder w.o worker number
            da.TEMP.num_iterations = 0; %was 1 before!!
            da.TEMP.recalculate_stratigraphy_now = 0;
            da.TEMP.recalculate_forcing_now = 0;

            da.ENSEMBLE.modeled_obs = []; %Yp in Kris code, size N_obs x N_ens x N_iterations
            da.ENSEMBLE.value_gaussian = [];  %Xp in Kris code, size N_param x N_ens x N_iterations

            run_mode_class = str2func(da.PARA.run_mode);
            da.IO = run_mode_class();
            da = make_scratchfolder(da.IO, da, tile);

            da = load_observable_operators(da, tile);

            da = get_ensemble_info(da, tile);
            da.TEMP.old_value_gaussian = tile.ENSEMBLE.TEMP.value_gaussian(da.TEMP.pos_in_ensemble,:);
            da.TEMP.old_mean_gaussian = tile.ENSEMBLE.TEMP.mean_gaussian(da.TEMP.pos_in_ensemble,:);
            da.TEMP.old_std_gaussian = tile.ENSEMBLE.TEMP.std_gaussian(da.TEMP.pos_in_ensemble,:);
            %fill the necessary inforation for first AMIS
            da.TEMP.cov_gaussian_resampled = diag(tile.ENSEMBLE.TEMP.std_gaussian(da.TEMP.pos_in_ensemble,1).^2);
            da.TEMP.mean_gaussian_resampled = tile.ENSEMBLE.TEMP.mean_gaussian(da.TEMP.pos_in_ensemble,1);

            da = init_DA(da,tile, 0); %load observations, determine DA_time and da_STEP_TIME
            if da.TEMP.assimilation_started
                da = save_prior_state(da.IO,da,tile);
            end
        end
        
        
        function da = DA_step(da, tile)
                        
            %save the state and the ensemble variables when the DA begins,
            %this step is redone every time the DA is happy and moves on in time
            if ~da.TEMP.assimilation_started && tile.t>=da.TEMP.last_assimilation_date
                da = init_DA(da,tile, 1);
                da = save_state_prior(da.IO,da,tile);
            end

            if da.TEMP.assimilation_started && tile.t>= da.DA_TIME
                da = get_modeled_observations(da, tile);
            
                %da = get_ancillary_features(da, tile); %TO DO
            end
            
            if tile.t>=da.DA_STEP_TIME

                da.TEMP.num_iterations = da.TEMP.num_iterations + 1; %num_iterations is now the value of the current iteration

                da = collect_modeled_observations(da.IO, da, tile);
                
                %da = collect_ancillary_features(da.IO, tile); %TO DO

                da = collect_observations(da, tile);

%---------------
                da = AMIS(da);
                    
                %store DA results, in case user wishes
                da = save_da_results_all(da.IO, da, tile);
                
                %start the iterative loop which either terminates or starts a new iteration, dependent on diversity of ensemmble
                if da.ENSEMBLE.effective_ensemble_size./tile.PARA.ensemble_size >= da.PARA.min_ensemble_diversity || da.TEMP.num_iterations>=da.PARA.max_iterations
                    %terminate and move on in time, i.e. do a normal resampling of model state and resample parameters according to the learning coefficient
                    
                    if da.ENSEMBLE.effective_ensemble_size./tile.PARA.ensemble_size >= da.PARA.min_ensemble_diversity
                        disp('successful, ensemble is sufficiently diverse')
                    else
                        disp('maximum number of iterations reached')
                    end
                    
                    da = save_da_results_final(da.IO, da, tile);
                    
                    da = save_state(da.IO, da, tile); %save the current run for resampling, index is num_iterations
                   
                    da = resample_state_and_parameters(da, tile);
                    
                    %call the transform function of ensemble
                    da.TILE.ENSEMBLE = recalculate_ensemble_parameters_after_DA(da.TILE.ENSEMBLE, tile, da.PARA.ensemble_variables);
                    
                    %assign next DA_STEP_TIME
                    da = get_DA_step_time(da, tile);

                    da.TEMP.last_assimilation_date = tile.t;
                    da.TEMP.num_iterations = 0; %reset number of iterations for the next DA period

                    da = reset_observation_time_indeces(da, tile);
                    
                    %reset modelled observations and vectors conatining parameters 
                    da.ENSEMBLE.modeled_obs = []; %Yp in Kris code, size N_obs x N_ens x N_iterations
                    da.ENSEMBLE.value_gaussian = [];  %Xp in Kris code, size N_param x N_ens x N_iterations                    
                    
                    %this is only valid for learning coefficient 0! 
                    % CHECK THIS! tile.EMSEMBLE is changed depending on
                    % learning coefficient?
                    da.TEMP.old_value_gaussian = tile.ENSEMBLE.TEMP.value_gaussian(da.TEMP.pos_in_ensemble,1);
                    da.TEMP.cov_gaussian_resampled = diag(tile.ENSEMBLE.TEMP.std_gaussian(da.TEMP.pos_in_ensemble,1).^2);
                    da.TEMP.mean_gaussian_resampled = tile.ENSEMBLE.TEMP.mean_gaussian(da.TEMP.pos_in_ensemble,1);
                
                    da = save_state(da.IO, da, tile); %saves state 0 for next asssimilation period 

                else
                    %iterate and go back to start of DA period, resample and inflate again, start over with "old" model states at the beginning of the DA period
                    disp('ensemble degenerate, one more time')

                    da = save_state(da.IO, da, tile); %save the state with the current iteration
                                       
                    %load "old" state, i.e. from 1st iteration
                    da = load_prior_state(da.IO, da, tile);
                    
                    %1. load feature and training data from OUT class which
                    %have been save in save_state above - this should
                    %contain the OUT class with all the info - feature =
                    %perturbed parameters; target = pseudo_obs
                    target_ML = reshape(da.ENSEMBLE.modeled_obs, size(da.ENSEMBLE.modeled_obs,1), size(da.ENSEMBLE.modeled_obs,2).*size(da.ENSEMBLE.modeled_obs,3))';
                    %this should be the conditional likelyood
            Rinv=(obs_variance').^(-1);
            Inn=repmat(observations,1,size(modeled_obs,2))-modeled_obs;   % Innovation.
            EObj=Rinv*(Inn.^2);                     % [1 x Ne] ensemble objective function.
                            Ypell=da.ENSEMBLE.modeled_obs(:,:,ell);
                residuell=da.ENSEMBLE.observations - Ypell;
                phidell=0.5*Rinv*(residuell.^2);

                    
                    feature_data_in_ML = reshape(da.ENSEMBLE.value_gaussian, size(da.ENSEMBLE.value_gaussian,1), size(da.ENSEMBLE.value_gaussian,2).*size(da.ENSEMBLE.value_gaussian,3))';
                    %here's where existing data from a previous DA could come in
                    da.STATVAR.out_mean = mean(target_ML(:),1);
                    da.STATVAR.out_std = std(target_ML(:),[],1);
                    da.STATVAR.in_mean = mean(feature_data_in_ML,1);
                    da.STATVAR.in_std = std(feature_data_in_ML,[],1);

                    da.STATVAR.out = (target_ML(:) - da.STATVAR.out_mean) ./ da.STATVAR.out_std;
                    da.STATVAR.in = (feature_data_in_ML - da.STATVAR.in_mean) ./ da.STATVAR.in_std;

                    %2. (re)train NN
                    da.ML = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(da.PARA.ml_class){da.PARA.ml_class_index,1});
                    da.ML = finalize_init(da.ML, da);
                    da.TEMP.var_ID = 1; %only one variable
                    training_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(da.PARA.training_class){da.PARA.training_class_index,1});
                    training_class = finalize_init(training_class, da);
                    tile.ML = train_ML(training_class, da);

                    %3. create large ensemble with original ensemble width
                    number_of_members = 2000;
                    large_ensemble4ML = generate_ensemble_gaussian_space(da.TILE.ENSEMBLE, tile, da.PARA.ensemble_variables, number_of_members); %must be coded

                    %4. forward run with NN for these parameters (large pseudo_obs data set)
                    [large_ensemble4ML_modeled_obs, ~] = progapagate_ML(da.ML, large_ensemble4ML_param);

                    %5. simple PBS to find optimal ensemble
                    [weights, effective_ensemble_size] = PBS2(da, observations, large_ensemble4ML_modeled_obs, obs_variance);

                    %6. randsample from this data set and set 
                    rng(da.DA_STEP_TIME.*da.TEMP.num_iterations+3);
                    resample_ID = randsample([1:number_of_members]', da.TILE.PARA.ensemble_size, true, weights);
                    value_gaussian_resampled = large_ensemble4ML(:,resample_ID);
                    da.TILE.ENSEMBLE.TEMP.value_gaussian(da.TEMP.pos_in_ensemble,1) = value_gaussian_resampled(:, tile.PARA.worker_number);

                    % da = resample_AMIS(da, tile);
                    da.TILE.ENSEMBLE = recalculate_ensemble_parameters_after_DA(da.TILE.ENSEMBLE, tile, da.PARA.ensemble_variables);

                    %reset timestamps, no need to reset timestamps in the CG stratigraphy since the old states are read in
                    tile.t = da.TEMP.last_assimilation_date;

                    da = reset_observation_time_indeces(da, tile);

                    % da = save_da_results_all(da.IO, da, tile);

                end
                if da.PARA.recalculate_stratigraphy ==1
                    da.TEMP.recalculate_stratigraphy_now = 1;
                end
            end
            
        end
    
    end
end

