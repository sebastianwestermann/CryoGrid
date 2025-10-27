classdef DA_MULTITILE_iterative_PBS < matlab.mixin.Copyable

    
    properties
        TILE
        OBS_OP
        PARA
        CONST
        STATVAR
        TEMP
        DA_TIME     
        DA_STEP_TIME
        ENSEMBLE
    end
    
    methods
        function da = provide_PARA(da)
            
            da.PARA.observation_files = [];
            da.PARA.observation_paths = [];
            da.PARA.observable_classes = [];
            da.PARA.observable_classes_index = []; %must all have the same length, i.e. each observational data set requires one observable class 
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
            da.PARA.recalculate_forcing = [];
        end
        
        function da = provide_CONST(da)
            
        end
        
        function da = provide_STATVAR(da)
            
        end
        
        function da = finalize_init(da, tile)
            da.TILE = tile;
            da.TEMP.run_name = da.TILE.PARA.run_name; %original results folder w.o worker number
            da.TEMP.first_obs_index =[];
            da.TEMP.index_next_obs = [];
            da.TEMP.time_next_obs = [];
            for i=1:size(da.PARA.observation_files,1)
                temp=load([da.PARA.observation_paths{i,1} da.PARA.observation_files{i,1}], 'OBS');
                da.STATVAR.obs_time{i,1} = temp.OBS.time;
                da.STATVAR.observations{i,1} = temp.OBS.observations;
                da.STATVAR.obs_variance{i,1} = temp.OBS.obs_variance;
                da.STATVAR.modeled_obs{i,1} = repmat(da.STATVAR.observations{i,1}.*NaN,1,tile.PARA.ensemble_size);
                da.ENSEMBLE.weights = repmat(1./da.TILE.PARA.ensemble_size, 1, da.TILE.PARA.ensemble_size); %set equal weights before 1st DA step
                
                da.OBS_OP{i,1} = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(da.PARA.observable_classes{i,1}){da.PARA.observable_classes_index(i,1)});     
                da.TEMP.first_obs_index = [da.TEMP.first_obs_index; find(da.STATVAR.obs_time{i,1} > tile.t, 1)];
                da.TEMP.index_next_obs = [da.TEMP.index_next_obs; da.TEMP.first_obs_index(i,1)]; %start with 
                da.TEMP.time_next_obs = [da.TEMP.time_next_obs; da.STATVAR.obs_time{i,1}(da.TEMP.first_obs_index(i,1),1)];
            end

            da.DA_TIME = min(da.TEMP.time_next_obs);
            
            if strcmp(da.PARA.assimilation_frequency, 'year')
                current_year = str2num(datestr(tile.t, 'yyyy'));
                da.DA_STEP_TIME = datenum([da.PARA.assimilation_date num2str(current_year + da.PARA.assimilation_interval)], 'dd.mm.yyyy');
            elseif strcmp(da.PARA.assimilation_frequency, 'month')
                current_year = str2num(datestr(tile.t, 'yyyy'));
                current_month = str2num(datestr(tile.t, 'mm'));
                da.DA_STEP_TIME = datenum(current_year, current_month + da.PARA.assimilation_interval, da.PARA.assimilation_date);
            elseif strcmp(da.PARA.assimilation_frequency, 'day')
                da.DA_STEP_TIME = tile.t + da.PARA.assimilation_interval;
            end
            
            if isempty(da.PARA.start_assimilation_period) || sum(isnan(da.PARA.start_assimilation_period))>0
                da = save_state(da, tile);
                da.TEMP.last_assimilation_date = tile.t+1; %start very close to initial state
                da.TEMP.assimilation_started = 0;
            else
                da.TEMP.assimilation_started = 0;
                da.TEMP.last_assimilation_date = datenum(da.PARA.start_assimilation_period(1,1), da.PARA.start_assimilation_period(2,1), da.PARA.start_assimilation_period(3,1));
            end
            
            da.TEMP.num_iterations = 1;
            
            %vector of positions to convert between the list of variables
            %that are changed by the DA to the full list of perturbed
            %variables in ENSEMBLE
            pos_in_ensemble = [];
            for i=1:size(da.PARA.ensemble_variables,1)
                 pos_in_ensemble = [pos_in_ensemble; find(strcmp(da.PARA.ensemble_variables{i,1}, tile.ENSEMBLE.TEMP.variable_name))];
            end
            da.TEMP.pos_in_ensemble = pos_in_ensemble;

            da.TEMP.recalculate_stratigraphy_now = 0;
            da.TEMP.recalculate_forcing_now = 0;
        end
        
        
        function da = DA_step(da, tile)
            %save the state and the ensemble variables when the DA begins,
            %this step is redone every time the DA is happy and moves on in time
            if ~da.TEMP.assimilation_started && tile.t>=da.TEMP.last_assimilation_date
                da.TEMP.last_assimilation_date = tile.t;
                da.TEMP.assimilation_started = 1;
                da.ENSEMBLE.weights_old = da.ENSEMBLE.weights;
                da = save_state(da, tile); %save states at the start of the assimilation period 
                da.TEMP.old_value_gaussian = tile.ENSEMBLE.TEMP.value_gaussian(da.TEMP.pos_in_ensemble,:);
                da.TEMP.old_mean_gaussian = tile.ENSEMBLE.TEMP.mean_gaussian(da.TEMP.pos_in_ensemble,:);
                da.TEMP.old_std_gaussian = tile.ENSEMBLE.TEMP.std_gaussian(da.TEMP.pos_in_ensemble,:);
            end
            
            
            if da.TEMP.assimilation_started && tile.t>= da.DA_TIME
                %loop over all observation data sets
                for i=1:size(da.STATVAR.obs_time,1)
                    if tile.t>= da.TEMP.time_next_obs(i,1)
                        da.STATVAR.modeled_obs{i,1}(da.TEMP.index_next_obs(i,1),:) = observable_operator(da.OBS_OP{i,1}, tile);
                        %disp(da.STATVAR.modeled_obs{i,1}(da.TEMP.index_next_obs(i,1),:) )
                        disp('collecting synthetic observations')
                        if da.TEMP.index_next_obs(i,1) < size(da.STATVAR.observations{i,1}, 1) %end of observations reached
                            da.TEMP.index_next_obs(i,1) = da.TEMP.index_next_obs(i,1) + 1;
                            da.TEMP.time_next_obs(i,1) = da.STATVAR.obs_time{i,1}(da.TEMP.index_next_obs(i,1),1);
                        else
                            da.TEMP.time_next_obs(i,1) = tile.FORCING.PARA.end_time +1;
                        end
                    end
                    
                end
                
                da.DA_TIME = min(da.TEMP.time_next_obs);
                
            end
            
            if tile.t>=da.DA_STEP_TIME
                %synchronize
                data_package = [];
                
                modeled_obs = [];%gather modeled observations in one vector
                for i=1:size(da.STATVAR.modeled_obs,1) 
                    da.ENSEMBLE.modeled_obs = [modeled_obs; da.STATVAR.modeled_obs{i,1}(da.TEMP.first_obs_index(i,1):da.TEMP.index_next_obs(i,1)-1,:)];  %ONLY USE THE PART IN THE OBS INTERVAL, ALSO MAKE A SIMILAR VECTOR FOR OBSERVATIONS AND VARIANCES
                end
                
                da.ENSEMBLE.value_gaussian = tile.ENSEMBLE.TEMP.value_gaussian(da.TEMP.pos_in_ensemble,:);
                               
%                 ensemble_param = []; %gather ensemble parameters in one vector
%                 variables = fieldnames(da.TILE.ENSEMBLE.PARA);
%                 for j=1:size(variables,1)
%                     da.ENSEMBLE.ensemble_param = [ensemble_param; da.TILE.ENSEMBLE.PARA.(variables{j,1})];
%                 end

                if da.TEMP.num_iterations == 1
                    da = PBS(da);
                else
                    da = adaptive_PBS(da);
                end
                
                if strcmp(da.PARA.store_format, 'all')
                    da_store = copy(da);
                    da_store.TILE = [];
                    if isempty(da.PARA.store_file_tag) || isnan(da.PARA.store_file_tag)
                        save([tile.PARA.result_path tile.PARA.run_name '/' 'da_store_'  datestr(tile.t, 'yyyymmdd') '_' num2str(da.TEMP.num_iterations)  '.mat'], 'da_store')
                    else
                        save([tile.PARA.result_path tile.PARA.run_name '/' 'da_store_'  datestr(tile.t, 'yyyymmdd') '_' num2str(da.TEMP.num_iterations) '_' da.PARA.store_file_tag '.mat'], 'da_store')
                    end
                end
                
                if da.ENSEMBLE.effective_ensemble_size./tile.PARA.ensemble_size >= da.PARA.min_ensemble_diversity || da.TEMP.num_iterations>=da.PARA.max_iterations
                    %do not iterate, but move on in time, do a normal
                    %resampling of model state and resample
                    %parameters according to the learning coefficient
                    
                    
                    if strcmp(da.PARA.store_format, 'final') 
                        da_store = copy(da);
                        da_store.TILE = [];
                        if isempty(da.PARA.store_file_tag) || isnan(da.PARA.store_file_tag)
                            save([tile.PARA.result_path tile.PARA.run_name '/' 'da_store_'  datestr(tile.t, 'yyyymmdd') '.mat'], 'da_store')
                        else
                            save([tile.PARA.result_path tile.PARA.run_name '/' 'da_store_' datestr(tile.t, 'yyyymmdd') '_' da.PARA.store_file_tag '.mat'], 'da_store')
                        end
                    end
                    
                    da.TEMP.num_iterations = 1;
                    resample_ID = randsample(da.TILE.PARA.ensemble_size, da.TILE.PARA.ensemble_size, true, da.ENSEMBLE.weights); %replaces "get_from_new_worker"
                    da = save_state(da, tile);
                                        
                    %read the new stratigraphy and info from file
%                     temp=load([tile.PARA.result_path  da.TEMP.run_name '/tile_' num2str(resample_ID(tile.PARA.worker_number,1)) '.mat']);
                    temp=load([tile.PARA.result_path  da.TEMP.run_name '/saved_state.mat']);
                    variables = fieldnames(temp.state);
                    for i=1:size(variables,1)
                        if size(temp.state.(variables{i,1}),2) == tile.PARA.ensemble_size
                            for j=1:tile.PARA.ensemble_size
                                tile.SUBSURFACE_CLASS.STATVAR.(variables{i,1})(:,j) = temp.state.(variables{i,1})(:,resample_ID(j,1));
                            end
                        end
                    end
                    
                    rand_sequence = rand(1,tile.PARA.ensemble_size);
                    %learning, use the inflation
                    value_gaussian_resampled = da.ENSEMBLE.value_gaussian(:,resample_ID);  %
                    mean_gaussian_resampled = mean(value_gaussian_resampled,2);   %propm=mean(thetap,2);

                    deviation = value_gaussian_resampled - mean_gaussian_resampled; % A=thetap-propm;
                    proposal_cov =  (1./tile.PARA.ensemble_size) .* (deviation*deviation');  %   propc=(1/Ne).*(A*A'); % Proposal covariance
                    % Inflate proposal covariance to avoid overconfidence

                    %proposal_cov = proposal_cov + 0.1 .* (1 - da.ENSEMBLE.effective_ensemble_size./ tile.PARA.ensemble_size).* diag(da.TEMP.old_std_gaussian); % propc=propc+0.1.*(1-diversity).*pric;
                    proposal_cov = proposal_cov + 0.1 .* (1 - min(1, da.ENSEMBLE.effective_ensemble_size./ tile.PARA.ensemble_size./da.PARA.min_ensemble_diversity)).* diag(da.TEMP.old_std_gaussian); % propc=propc+0.1.*(1-diversity).*pric;
                    % Gaussian resampling using the Cholesky decomposition
                    L=chol(proposal_cov,'lower');
                    Z = randn(size(mean_gaussian_resampled,1), tile.PARA.ensemble_size); %Z=randn(Np,Ne);
                    value_gaussian_resampled = mean_gaussian_resampled + L*Z;

                    for j=1:tile.PARA.ensemble_size
                        if da.PARA.learning_coefficient >= rand_sequence(1, j)
                            tile.ENSEMBLE.TEMP.value_gaussian(da.TEMP.pos_in_ensemble,j) = value_gaussian_resampled(:,j);
                        else
                            tile.ENSEMBLE.TEMP.value_gaussian(da.TEMP.pos_in_ensemble,j) = da.TEMP.old_value_gaussian(:,j);
                        end
                    end
                    
                    %call the transform function of ensemble
                    
                    da.TILE.ENSEMBLE = recalculate_ensemble_parameters_after_DA(da.TILE.ENSEMBLE, tile, da.PARA.ensemble_variables);
                    
                    %assign next DA_STEP_TIME
                    if strcmp(da.PARA.assimilation_frequency, 'year')
                        current_year = str2num(datestr(da.DA_STEP_TIME, 'yyyy'));
                        da.DA_STEP_TIME = datenum([da.PARA.assimilation_date num2str(current_year + da.PARA.assimilation_interval)], 'dd.mm.yyyy');
                    elseif strcmp(da.PARA.assimilation_frequency, 'month')
                        current_year = str2num(datestr(da.DA_STEP_TIME, 'yyyy'));
                        current_month = str2num(datestr(da.DA_STEP_TIME, 'mm'));
                        da.DA_STEP_TIME = datenum(current_year, current_month + da.PARA.assimilation_interval, da.PARA.assimilation_date);
                    elseif strcmp(da.PARA.assimilation_frequency, 'day')
                        da.DA_STEP_TIME = da.DA_STEP_TIME + da.PARA.assimilation_interval;
                    end
                    
                    da.TEMP.first_obs_index =[];
                    da.TEMP.index_next_obs = [];
                    da.TEMP.time_next_obs = [];
                    for i=1:size(da.PARA.observation_files,1)
                        if ~isempty(find(da.STATVAR.obs_time{i,1} > tile.t, 1))
                            da.TEMP.first_obs_index = [da.TEMP.first_obs_index; find(da.STATVAR.obs_time{i,1} > tile.t, 1)];
                            da.TEMP.index_next_obs = [da.TEMP.index_next_obs; da.TEMP.first_obs_index(i,1)]; %start with
                            da.TEMP.time_next_obs = [da.TEMP.time_next_obs; da.STATVAR.obs_time{i,1}(da.TEMP.first_obs_index(i,1),1)];
                        else
                            da.TEMP.first_obs_index = [da.TEMP.first_obs_index; size(da.STATVAR.observations{i,1}, 1)];
                            da.TEMP.index_next_obs = [da.TEMP.index_next_obs; size(da.STATVAR.observations{i,1}, 1)];
                            da.TEMP.time_next_obs = [da.TEMP.time_next_obs; tile.FORCING.PARA.end_time + 1];
                        end
                    end
                    
                    da.DA_TIME = min(da.TEMP.time_next_obs);
                    da.TEMP.last_assimilation_date = tile.t;
                    da.ENSEMBLE.weights_old = da.ENSEMBLE.weights;
                    da = save_state(da, tile); %save new states at the start of the new assimilation period, so that it can be read again if ensemble is degenerate
%                     for i=1:size(da.PARA.ensemble_variables,1)
%                         da.TEMP.old_ensemble.(da.PARA.ensemble_variables{i,1}) = get_variable_info(tile.ENSEMBLE, da.PARA.ensemble_variables{i,1});
%                     end
                    
                    da.TEMP.old_value_gaussian = tile.ENSEMBLE.TEMP.value_gaussian(da.TEMP.pos_in_ensemble,:);
%                     da.TEMP.old_mean_gaussian = tile.ENSEMBLE.TEMP.mean_gaussian(logical(da.TEMP.change_yes_no),1);
%                     da.TEMP.old_std_gaussian = tile.ENSEMBLE.TEMP.std_gaussian(logical(da.TEMP.change_yes_no),1);
                    %assimilation successful, ensemble is not degenerate,
                    %move on in time
                    
                else
                    %iterate and go back to start of DA period, resample
                    %and inflate again, start over with "old" model states
                    %at the beginning of the DA period
                    disp('ensemble degenerate, one more time')
                    
                    da.TEMP.num_iterations = da.TEMP.num_iterations + 1;
                    
                    %load "old" state            
                    temp=load([tile.PARA.result_path  da.TEMP.run_name '/saved_state.mat']);
                    variables = fieldnames(temp.state);
                    for i=1:size(variables,1)
                        if size(temp.state.(variables{i,1}),2) == tile.PARA.ensemble_size
                            for j=1:tile.PARA.ensemble_size
                                tile.SUBSURFACE_CLASS.STATVAR.(variables{i,1})= temp.state.(variables{i,1});
                            end
                        end
                    end
                    
                    resample_ID = randsample(da.TILE.PARA.ensemble_size, da.TILE.PARA.ensemble_size, true, da.ENSEMBLE.weights); %replaces "get_from_new_worker"
                    
                    value_gaussian_resampled = da.ENSEMBLE.value_gaussian(:,resample_ID);  %    thetap=thetap(:,resample);
                    mean_gaussian_resampled = mean(value_gaussian_resampled,2);   %propm=mean(thetap,2);
                    
                    deviation = value_gaussian_resampled - mean_gaussian_resampled; % A=thetap-propm;
                    proposal_cov =  (1./tile.PARA.ensemble_size) .* (deviation*deviation');  %   propc=(1/Ne).*(A*A'); % Proposal covariance
                    % Inflate proposal covariance to avoid overconfidence
                    
                    %proposal_cov = proposal_cov + 0.1 .* (1 - da.ENSEMBLE.effective_ensemble_size./ tile.PARA.ensemble_size).* diag(da.TEMP.old_std_gaussian); % propc=propc+0.1.*(1-diversity).*pric;
                    proposal_cov = proposal_cov + 0.1 .* (1 - min(1, da.ENSEMBLE.effective_ensemble_size./ tile.PARA.ensemble_size./da.PARA.min_ensemble_diversity)).* diag(da.TEMP.old_std_gaussian); % propc=propc+0.1.*(1-diversity).*pric;
                    % Gaussian resampling using the Cholesky decomposition
                    L=chol(proposal_cov,'lower');
                    Z = randn(size(mean_gaussian_resampled,1), tile.PARA.ensemble_size); %Z=randn(Np,Ne);
                    value_gaussian_resampled = mean_gaussian_resampled + L*Z;

                    tile.ENSEMBLE.TEMP.value_gaussian(da.TEMP.pos_in_ensemble,:) = value_gaussian_resampled;
                              
                    da.TILE.ENSEMBLE = recalculate_ensemble_parameters_after_DA(da.TILE.ENSEMBLE, tile, da.PARA.ensemble_variables);

                    %reset timestamps, no need to reset timestamps in the
                    %CG stratigraphy since the old states are read in
                    tile.t = da.TEMP.last_assimilation_date;

                    da.TEMP.first_obs_index =[];
                    da.TEMP.index_next_obs = [];
                    da.TEMP.time_next_obs = [];
                    for i=1:size(da.PARA.observation_files,1)
                        if ~isempty(find(da.STATVAR.obs_time{i,1} > tile.t, 1))
                            da.TEMP.first_obs_index = [da.TEMP.first_obs_index; find(da.STATVAR.obs_time{i,1} > tile.t, 1)];
                            da.TEMP.index_next_obs = [da.TEMP.index_next_obs; da.TEMP.first_obs_index(i,1)]; %start with
                            da.TEMP.time_next_obs = [da.TEMP.time_next_obs; da.STATVAR.obs_time{i,1}(da.TEMP.first_obs_index(i,1),1)];
                        else
                            da.TEMP.first_obs_index = [da.TEMP.first_obs_index; size(da.STATVAR.observations{i,1}, 1)];
                            da.TEMP.index_next_obs = [da.TEMP.index_next_obs; size(da.STATVAR.observations{i,1}, 1)];
                            da.TEMP.time_next_obs = [da.TEMP.time_next_obs; tile.FORCING.PARA.end_time + 1];
                        end
                    end
                    
                    da.DA_TIME = min(da.TEMP.time_next_obs);
                    tile.OUT = reset_timestamp_out(tile.OUT,tile);
                end
                if da.PARA.recalculate_stratigraphy ==1
                    da.TEMP.recalculate_stratigraphy_now = 1;
                end
                if da.PARA.recalculate_forcing ==1
                    da.TEMP.recalculate_forcing_now = 1;
                end
            end
%                 
%                 %resampling from weights
%                 rng(tile.t+25) %use current time as seed for randum number generator, i.e. same sequence of random numbers will be generated by each worker
%                 [weights_sorted, posi] = sort(da.ENSEMBLE.weights);
%                 weights_sorted_cum = cumsum(weights_sorted);
%                 weights_sorted_cum = [0 weights_sorted_cum(1, 1:end-1)];
%                 resample_posi=sum(double(repmat(rand(da.TILE.PARA.ensemble_size,1),1,da.TILE.PARA.ensemble_size) > repmat(weights_sorted_cum, da.TILE.PARA.ensemble_size,1)),2);
%                 get_from_new_worker = [];
%                 for i=1:da.TILE.PARA.ensemble_size
%                     get_from_new_worker=[get_from_new_worker; posi(resample_posi(i))];
%                 end
%                 
%                 variables = fieldnames(tile.SUBSURFACE_CLASS.STATVAR);
%                 temp = tile.SUBSURFACE_CLASS.STATVAR;
%                 for i=1:da.TILE.PARA.ensemble_size
%                     for j=1:size(variables,1)
%                         if size(tile.SUBSURFACE_CLASS.STATVAR.(variables{j,1}),2) == tile.PARA.ensemble_size
%                             tile.SUBSURFACE_CLASS.STATVAR.(variables{j,1})(:,i) = temp.(variables{j,1})(:,get_from_new_worker(i,1));
%                         end
%                     end
%                 end
%                 


%                 %assign next DA_STEP_TIME
%                 if strcmp(da.PARA.assimilation_frequency, 'year')
%                     current_year = str2num(datestr(da.DA_STEP_TIME, 'yyyy'));
%                     da.DA_STEP_TIME = datenum([da.PARA.assimilation_date num2str(current_year + da.PARA.assimilation_interval)], 'dd.mm.yyyy');
%                 elseif strcmp(da.PARA.assimilation_frequency, 'month')
%                     current_year = str2num(datestr(da.DA_STEP_TIME, 'yyyy'));
%                     current_month = str2num(datestr(da.DA_STEP_TIME, 'mm'));
%                     da.DA_STEP_TIME = datenum(current_year, current_month + da.PARA.assimilation_interval, da.PARA.assimilation_date);
%                 elseif strcmp(da.PARA.assimilation_frequency, 'day')
%                     da.DA_STEP_TIME = da.DA_STEP_TIME + da.PARA.assimilation_interval;
%                 end
%                 
%                 
%                 da.TEMP.first_obs_index =[];
%                 da.TEMP.index_next_obs = [];
%                 da.TEMP.time_next_obs = [];
%                 for i=1:size(da.PARA.observation_files,1)
%                     if ~isempty(find(da.STATVAR.obs_time{i,1} > tile.t, 1))
%                         da.TEMP.first_obs_index = [da.TEMP.first_obs_index; find(da.STATVAR.obs_time{i,1} > tile.t, 1)];
%                         da.TEMP.index_next_obs = [da.TEMP.index_next_obs; da.TEMP.first_obs_index(i,1)]; %start with
%                         da.TEMP.time_next_obs = [da.TEMP.time_next_obs; da.STATVAR.obs_time{i,1}(da.TEMP.first_obs_index(i,1),1)];
%                     else
%                         da.TEMP.first_obs_index = [da.TEMP.first_obs_index; size(da.STATVAR.observations{i,1}, 1)];
%                         da.TEMP.index_next_obs = [da.TEMP.index_next_obs; size(da.STATVAR.observations{i,1}, 1)];
%                         da.TEMP.time_next_obs = [da.TEMP.time_next_obs; tile.FORCING.PARA.end_time + 1];
%                     end
%                 end
%                 
%                 da.DA_TIME = min(da.TEMP.time_next_obs);


            
        end
        
       function da = PBS(da)
           % w = PBS( HX,Y,R )
            % Efficient implementation of the Particle Batch Smoother
            % presented in Margulis et al. (2015; JHM).
            % N.B. The observation errors are assumed to be uncorrelated (diagonal R)
            % and Gaussian.
            %
            % Dimensions: No = Number of observations in the batch to assimilate.
            %             Np = Number of parameters to update.
            %             Ne = Number of ensemble members (particles).
            %
            % -----------------------------------------------------------------------
            % Inputs:
            %
            %
            % HX   => No x Ne matrix containing an ensemble of Ne predicted
            %         observation column vectors each with No entries.
            %
            % Y     => No x 1 vector containing the batch of (unperturbed) observations.
            %
            % R     => No x No observation error variance matrix; this may also be
            %         specified as a scalar corresponding to the constant variance of
            %         all the observations in the case that these are all from the same
            %         instrument.
            %
            % -----------------------------------------------------------------------
            % Outputs:
            %
            % w     => 1 x Ne vector containing the ensemble of posterior weights,
            %         the prior weights are implicitly 1/N_e.
            %
            % -----------------------------------------------------------------------
            % See e.g. https://jblevins.org/log/log-sum-exp for underflow issue.
            %
            % Code by Kristoffer Aalstad (Feb. 2019)
            
            % Calculate the diagonal of the inverse obs. error covariance.
            observations = [];
            obs_variance = [];
            for i=1:size(da.STATVAR.observations,1)
                observations = [observations; da.STATVAR.observations{i,1}(da.TEMP.first_obs_index(i,1):da.TEMP.index_next_obs(i,1)-1,1)];
                obs_variance = [obs_variance; da.STATVAR.obs_variance{i,1}(da.TEMP.first_obs_index(i,1):da.TEMP.index_next_obs(i,1)-1,1)];
            end
            da.ENSEMBLE.observations = observations;
            da.ENSEMBLE.obs_variance = obs_variance;
            
            No=size(observations,1);
            Rinv=(obs_variance').^(-1);
            
            % Calculate the likelihood.
            Inn=repmat(observations,1,size(da.ENSEMBLE.modeled_obs,2))-da.ENSEMBLE.modeled_obs;   % Innovation.
            EObj=Rinv*(Inn.^2);                     % [1 x Ne] ensemble objective function.
            LLH=-0.5.*EObj; % log-likelihoods.
            normc=logsumexp(da, LLH,2);
            
            
            % NB! The likelihood coefficient (1/sqrt(2*pi...)) is
            % omitted because it drops out in the normalization
            % of the likelihood. Including it (very small term) would lead
            % to problems with FP division.
            
            
            % Calculate the posterior weights as the normalized likelihood.
            logw = LLH-normc;
            da.ENSEMBLE.weights = exp(logw); % Posterior weights.
            
            % Need "log-sum-exp" trick to overcome numerical issues for small R/large
            % number of obs.
            da.ENSEMBLE.effective_ensemble_size = 1./sum(da.ENSEMBLE.weights.^2);
        end
            

            
        function s = logsumexp(da, a, dim)
            % Returns log(sum(exp(a),dim)) while avoiding numerical underflow.
            % Default is dim = 1 (columns).
            % logsumexp(a, 2) will sum across rows instead of columns.
            % Unlike matlab's "sum", it will not switch the summing direction
            % if you provide a row vector.
            
            % Written by Tom Minka
            % (c) Microsoft Corporation. All rights reserved.
            
            if nargin < 2
                dim = 1;
            end
            
            % subtract the largest in each column
            y = max(a,[],dim);
            dims = ones(1,ndims(a));
            dims(dim) = size(a,dim);
            a = a - repmat(y, dims);
            s = y + log(sum(exp(a),dim));
            i = find(~isfinite(y));
            if ~isempty(i)
                s(i) = y(i);
            end
        end
        
        function da = adaptive_PBS(da)
            %[w,Neff]= AdaPBS( Ypred, yobs, R, prim, pricov, proposal )
            % An Adaptive Particle Batch Smoother using a Gaussian proposal
            % N.B. The observation errors are assumed to be uncorrelated (diagonal R)
            % and Gaussian. This can easily be changed.
            %
            % Dimensions: No = Number of observations in the batch to assimilate.
            %             Np = Number of parameters to update.
            %             Ne = Number of ensemble members (particles).
            %
            % -----------------------------------------------------------------------
            % Inputs:
            %
            %
            % Ypred   => No x Ne matrix containing an ensemble of Ne predicted
            %         observation column vectors each with No entries.
            %
            % y     => No x 1 vector containing the batch of (unperturbed) observations.
            %
            % R     => No x 1 observation error variance vector; this may also be
            %         specified as a scalar corresponding to the constant variance of
            %         all the observations in the case that these are all from the same
            %         instrument.
            % prim  => Np x 1 Prior mean vector.
            %
            % pricov => Np x Np Prior covariance matrix of the paramters - typically
            % diagonal matrix w. variances of each parameter
            %
            % proposal => Np x Ne Samples from the proposal (assumed to be Gaussian) -
            % comes out of the loop
            %
            % -----------------------------------------------------------------------
            % Outputs:
            %
            % w     => 1 x Ne vector containing the ensemble of posterior weights,
            %         the prior weights are implicitly 1/N_e.
            %
            % -----------------------------------------------------------------------
            % See e.g. https://jblevins.org/log/log-sum-exp for underflow issue.
            %
            % Code by Kristoffer Aalstad (June 2023)

            
            observations = [];
            obs_variance = [];
            for i=1:size(da.STATVAR.observations,1)
                observations = [observations; da.STATVAR.observations{i,1}(da.TEMP.first_obs_index(i,1):da.TEMP.index_next_obs(i,1)-1,1)];
                obs_variance = [obs_variance; da.STATVAR.obs_variance{i,1}(da.TEMP.first_obs_index(i,1):da.TEMP.index_next_obs(i,1)-1,1)];
            end
            da.ENSEMBLE.observations = observations;
            da.ENSEMBLE.obs_variance = obs_variance;
            
            No=size(observations,1);
            Rinv=(obs_variance').^(-1);
            
            proposal = da.ENSEMBLE.value_gaussian;
            prim = da.TEMP.old_mean_gaussian;
            pricov = diag(da.TEMP.old_std_gaussian);
            
           % None of these three exist yet, compute here!
            % Prior term
            A0=proposal-prim;
            Phi0=-0.5.*(A0')*(pricov\A0);
            Phi0=diag(Phi0);
            Phi0=Phi0';

            % Proposal term
            A=proposal-mean(proposal,2);
            propcov=(1./da.TILE.PARA.ensemble_size).*(A*A');
            Phip=-0.5.*(A')*(propcov\A);
            Phip=diag(Phip);
            Phip=Phip';

            % Likelihood term
            residual = repmat(observations,1,size(da.ENSEMBLE.modeled_obs,2))-da.ENSEMBLE.modeled_obs; 
            Phid=-0.5.*Rinv*(residual.^2);

            Phi=Phid+Phi0-Phip;
            Phimax=max(Phi);
            Phis=Phi-Phimax; % Scaled to avoid numerical overflow (see Chopin book).

            w=exp(Phis);
            w=w./sum(w);

            Neff=1./(sum(w.^2));

            da.ENSEMBLE.weights = w;
            da.ENSEMBLE.effective_ensemble_size = Neff;
        end

        
        
        function data_package = pack(da, data_package, var_name, var) %transform into column vector ready to send
                %variables{i,1}
                data_package=[data_package; size(var_name,2); double(var_name)']; % # of characters followed by characters as doubles
                data_package=[data_package; size(var,1); var]; % # of entries followed by values
        end
        
        function da = unpack(da, data_package, received_from_worker) %read received column vector and transform into STATVAR
            i=1;
            while i<=size(data_package,1)
               variable_name = char(data_package(i+1:i+data_package(i,1),1)');
               i = i + data_package(i,1) + 1;
               da.ENSEMBLE.(variable_name)(:,received_from_worker) = data_package(i+1:i+data_package(i,1),1);
               i = i + data_package(i,1) + 1;
            end
        end

        
        function da = save_state(da, tile)
            state = tile.SUBSURFACE_CLASS.STATVAR;
            
            save([tile.PARA.result_path da.TEMP.run_name '/saved_state.mat'], 'state');
            
        end
    end
end

