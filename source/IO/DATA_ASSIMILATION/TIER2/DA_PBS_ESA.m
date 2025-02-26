%made for snow extent DA in ESA CCI, runs for multiple grid cells simultaneously 

%good to the point that PBS produces weights , bst-fitting ensemble member
%can be identified and written in outoput file

classdef DA_PBS_ESA < DA_FUNCTIONS


    properties

    end
    
    methods
        function da = provide_PARA(da)
            da.PARA.run_mode = 'DA_IO_TILE'; %or DA_IO_MULTITILE
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
            da.PARA.learning_coefficient = [];
            da.PARA.min_ensemble_diversity = [];
            da.PARA.max_iterations = [];
            da.PARA.save_options = []; %PROVIDER or FILE
            da.PARA.save_best_parameters_class = [];
            da.PARA.save_best_parameters_class_index = [];
            da.PARA.save_best_parameters_folder = [];
            da.PARA.recalculate_stratigraphy = 0;
            da.PARA.recalculate_forcing = 0;
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
            % da.TEMP.cov_gaussian_resampled = []; %size N_param x N_param x N_iterations
            % da.TEMP.mean_gaussian_resampled = []; %size N_param x N_iterations

            da_IO = str2func(da.PARA.run_mode);
            da.IO = da_IO();
            da = make_scratchfolder(da.IO, da, tile);

            da = load_observable_operators(da, tile);

            da = get_ensemble_info(da, tile);
            da.TEMP.old_value_gaussian = tile.ENSEMBLE.TEMP.value_gaussian(da.TEMP.pos_in_ensemble,:);
            da.TEMP.old_mean_gaussian = tile.ENSEMBLE.TEMP.mean_gaussian(da.TEMP.pos_in_ensemble,:);
            da.TEMP.old_std_gaussian = tile.ENSEMBLE.TEMP.std_gaussian(da.TEMP.pos_in_ensemble,:);

            da = init_DA(da,tile, 0); %load observations, determine DA_time and da_STEP_TIME
            %prepare store options
            if ~isempty(da.PARA.save_options)
                for i=1:size(da.PARA.save_options,1)
                    if strcmp(da.PARA.save_options{i,1}, 'PROVIDER')
                        for j=1:size(da.PARA.ensemble_variables,1)
                            tile.RUN_INFO.PPROVIDER.CLASSES.(da.PARA.save_best_parameters_class){da.PARA.save_best_parameters_class_index,1}.PARA.scale_parameters = {};
                            % tile.STORE.(['DA_' da.PARA.ensemble_variables{j,1} '_years']) = [];% out in PROVIDER instead, assign to the target class a PARA and make this target class scale_FORCING_DATA, which accepts the fields and the scaling variables as PARA
                            % tile.STORE.(['DA_' da.PARA.ensemble_variables{j,1}]) = []; % add PARA here, so that the target class can be selected
                        end
                    end
                end
            end
        end
        
        
        function da = DA_step(da, tile)
            if ~da.TEMP.assimilation_started && tile.t>=da.TEMP.last_assimilation_date
                da = init_DA(da,tile, 1);
            end
            
            if da.TEMP.assimilation_started && tile.t>= da.DA_TIME
                da = get_modeled_observations(da, tile);
            end
            
            if tile.t>=da.DA_STEP_TIME

                da.TEMP.num_iterations = da.TEMP.num_iterations + 1; %num_iterations is now the value of the current iteration

                da = collect_modeled_observations(da.IO, da, tile);
                da = collect_observations(da, tile);
                
%                 da.ENSEMBLE.value_gaussian = tile.ENSEMBLE.TEMP.value_gaussian(da.TEMP.pos_in_ensemble,:);
                               

                if da.TEMP.num_iterations == 1
                    da = PBS(da, tile);
                else
                    da = adaptive_PBS(da, tile);
                end
                
                best_parameters = get_best_fitting_parameters(tile.ENSEMBLE, tile, da.PARA.ensemble_variables, da.ENSEMBLE.weights);
                if ~isempty(da.PARA.save_options)
                    for i=1:size(da.PARA.save_options,1)
                        if strcmp(da.PARA.save_options{i,1}, 'PROVIDER')
                            best_parameters.year = year(tile.t)-1;
                            tile.RUN_INFO.PPROVIDER.CLASSES.(da.PARA.save_best_parameters_class){da.PARA.save_best_parameters_class_index,1}.PARA.scale_parameters = ...
                                [tile.RUN_INFO.PPROVIDER.CLASSES.(da.PARA.save_best_parameters_class){da.PARA.save_best_parameters_class_index,1}.PARA.scale_parameters;  best_parameters];
                        elseif strcmp(da.PARA.save_options{i,1}, 'FILE')
                            best_parameters.effective_ensemble_size = da.ENSEMBLE.effective_ensemble_size;
                            best_parameters.timestamp = tile.t;
                            if exist([da.PARA.save_best_parameters_folder 'best_parameters_' num2str(tile.PARA.range(1)) '_' num2str(tile.PARA.range(end)) '.mat'])==2
                                load([da.PARA.save_best_parameters_folder 'best_parameters_' num2str(tile.PARA.range(1)) '_' num2str(tile.PARA.range(end)) '.mat'])
                            else
                                 da_param= {};
                            end
                            da_param = [da_param; best_parameters];
                            save([da.PARA.save_best_parameters_folder 'best_parameters_' num2str(tile.PARA.range(1)) '_' num2str(tile.PARA.range(end)) '.mat'], 'da_param')
                        elseif strcmp(da.PARA.save_options{i,1}, 'FORCE_OUT')
                            tile.OUT = save_out_ensemble(tile.OUT, tile, da.ENSEMBLE.weights, best_parameters);
                        end
                    end
                end

                %reset SWE to 0, or read from the successful members, an
                %open new observation data set to continue
                
                
%                 if strcmp(da.PARA.store_format, 'all')
%                     da_store = copy(da);
%                     da_store.TILE = [];
%                     if isempty(da.PARA.store_file_tag) || isnan(da.PARA.store_file_tag)
%                         save([tile.PARA.result_path tile.PARA.run_name '/' 'da_store_'  datestr(tile.t, 'yyyymmdd') '_' num2str(da.TEMP.num_iterations)  '.mat'], 'da_store')
%                     else
%                         save([tile.PARA.result_path tile.PARA.run_name '/' 'da_store_'  datestr(tile.t, 'yyyymmdd') '_' num2str(da.TEMP.num_iterations) '_' da.PARA.store_file_tag '.mat'], 'da_store')
%                     end
%                 end
%                 
%                 if da.TEMP.num_iterations>=da.PARA.max_iterations %da.ENSEMBLE.effective_ensemble_size./tile.ENSEMBLE.PARA.grid_ensemble_size  >= da.PARA.min_ensemble_diversity || 
%                     %do not iterate, but move on in time, do a normal
%                     %resampling of model state and resample
%                     %parameters according to the learning coefficient
%                     
%                     
%                     if strcmp(da.PARA.store_format, 'final') 
%                         da_store = copy(da);
%                         da_store.TILE = [];
%                         if isempty(da.PARA.store_file_tag) || isnan(da.PARA.store_file_tag)
%                             save([tile.PARA.result_path tile.PARA.run_name '/' 'da_store_'  datestr(tile.t, 'yyyymmdd') '.mat'], 'da_store')
%                         else
%                             save([tile.PARA.result_path tile.PARA.run_name '/' 'da_store_' datestr(tile.t, 'yyyymmdd') '_' da.PARA.store_file_tag '.mat'], 'da_store')
%                         end
%                     end
%                     
%                     da.TEMP.num_iterations = 1;
%                     
%                     resample_ID = randsample(tile.ENSEMBLE.PARA.grid_ensemble_size , tile.ENSEMBLE.PARA.grid_ensemble_size , true, da.ENSEMBLE.weights); %replaces "get_from_new_worker"
%                     da = save_state(da, tile);
%                                         
%                     %read the new stratigraphy and info from file
%                     temp=load([tile.PARA.result_path  da.TEMP.run_name '/saved_state.mat']);
%                     variables = fieldnames(temp.state);
%                     for i=1:size(variables,1)
%                         if size(temp.state.(variables{i,1}),2) == tile.ENSEMBLE.PARA.grid_ensemble_size 
%                             for j=1:tile.ENSEMBLE.PARA.grid_ensemble_size 
%                                 tile.SUBSURFACE_CLASS.STATVAR.(variables{i,1})(:,j) = temp.state.(variables{i,1})(:,resample_ID(j,1));
%                             end
%                         end
%                     end
%                     
%                     rand_sequence = rand(1,tile.ENSEMBLE.PARA.grid_ensemble_size );
%                     %learning, use the inflation
%                     value_gaussian_resampled = da.ENSEMBLE.value_gaussian(:,resample_ID);  %
%                     mean_gaussian_resampled = mean(value_gaussian_resampled,2);   %propm=mean(thetap,2);
% 
%                     deviation = value_gaussian_resampled - mean_gaussian_resampled; % A=thetap-propm;
%                     proposal_cov =  (1./tile.ENSEMBLE.PARA.grid_ensemble_size ) .* (deviation*deviation');  %   propc=(1/Ne).*(A*A'); % Proposal covariance
%                     % Inflate proposal covariance to avoid overconfidence
% 
%                     %proposal_cov = proposal_cov + 0.1 .* (1 - da.ENSEMBLE.effective_ensemble_size./ tile.ENSEMBLE.PARA.grid_ensemble_size ).* diag(da.TEMP.old_std_gaussian); % propc=propc+0.1.*(1-diversity).*pric;
%                     proposal_cov = proposal_cov + 0.1 .* (1 - min(1, da.ENSEMBLE.effective_ensemble_size./ tile.ENSEMBLE.PARA.grid_ensemble_size ./da.PARA.min_ensemble_diversity)).* diag(da.TEMP.old_std_gaussian); % propc=propc+0.1.*(1-diversity).*pric;
%                     % Gaussian resampling using the Cholesky decomposition
%                     L=chol(proposal_cov,'lower');
%                     Z = randn(size(mean_gaussian_resampled,1), tile.ENSEMBLE.PARA.grid_ensemble_size ); %Z=randn(Np,Ne);
%                     value_gaussian_resampled = mean_gaussian_resampled + L*Z;
% 
%                     for j=1:tile.ENSEMBLE.PARA.grid_ensemble_size 
%                         if da.PARA.learning_coefficient >= rand_sequence(1, j)
%                             tile.ENSEMBLE.TEMP.value_gaussian(da.TEMP.pos_in_ensemble,j) = value_gaussian_resampled(:,j);
%                         else
%                             tile.ENSEMBLE.TEMP.value_gaussian(da.TEMP.pos_in_ensemble,j) = da.TEMP.old_value_gaussian(:,j);
%                         end
%                     end
%                     
%                     %call the transform function of ensemble
%                     
%                     da.TILE.ENSEMBLE = recalculate_ensemble_parameters_after_DA(da.TILE.ENSEMBLE, tile, da.PARA.ensemble_variables);
%                     
                    %reassign successful states?
                    da.TEMP.last_assimilation_date = tile.t;

                    da.TEMP.num_iterations = 0;
                    da = get_DA_step_time(da, tile);
                    
                    da = reset_observation_time_indeces(da, tile);

%                    da.ENSEMBLE.weights_old = da.ENSEMBLE.weights;
                    da = save_state(da, tile); %save new states at the start of the new assimilation period, so that it can be read again if ensemble is degenerate
                    
                   % da.TEMP.old_value_gaussian = tile.ENSEMBLE.TEMP.value_gaussian(da.TEMP.pos_in_ensemble,:);
% 
%                 else
%                     %not implemented yet, only one iteration allowed!
%                     
%                     %iterate and go back to start of DA period, resample
%                     %and inflate again, start over with "old" model states
%                     %at the beginning of the DA period
%                     disp('ensemble degenerate, one more time')
%                     
%                     da.TEMP.num_iterations = da.TEMP.num_iterations + 1;
%                     
%                     %load "old" state            
%                     temp=load([tile.PARA.result_path  da.TEMP.run_name '/saved_state.mat']);
%                     variables = fieldnames(temp.state);
%                     for i=1:size(variables,1)
%                         if size(temp.state.(variables{i,1}),2) == tile.ENSEMBLE.PARA.grid_ensemble_size 
%                             for j=1:tile.ENSEMBLE.PARA.grid_ensemble_size 
%                                 tile.SUBSURFACE_CLASS.STATVAR.(variables{i,1})= temp.state.(variables{i,1});
%                             end
%                         end
%                     end
%                     
%                     resample_ID = randsample(tile.ENSEMBLE.PARA.grid_ensemble_size , tile.ENSEMBLE.PARA.grid_ensemble_size , true, da.ENSEMBLE.weights); %replaces "get_from_new_worker"
%                     
%                     value_gaussian_resampled = da.ENSEMBLE.value_gaussian(:,resample_ID);  %    thetap=thetap(:,resample);
%                     mean_gaussian_resampled = mean(value_gaussian_resampled,2);   %propm=mean(thetap,2);
%                     
%                     deviation = value_gaussian_resampled - mean_gaussian_resampled; % A=thetap-propm;
%                     proposal_cov =  (1./tile.ENSEMBLE.PARA.grid_ensemble_size ) .* (deviation*deviation');  %   propc=(1/Ne).*(A*A'); % Proposal covariance
%                     % Inflate proposal covariance to avoid overconfidence
%                     
%                     %proposal_cov = proposal_cov + 0.1 .* (1 - da.ENSEMBLE.effective_ensemble_size./ tile.ENSEMBLE.PARA.grid_ensemble_size ).* diag(da.TEMP.old_std_gaussian); % propc=propc+0.1.*(1-diversity).*pric;
%                     proposal_cov = proposal_cov + 0.1 .* (1 - min(1, da.ENSEMBLE.effective_ensemble_size./ tile.ENSEMBLE.PARA.grid_ensemble_size ./da.PARA.min_ensemble_diversity)).* diag(da.TEMP.old_std_gaussian); % propc=propc+0.1.*(1-diversity).*pric;
%                     % Gaussian resampling using the Cholesky decomposition
%                     L=chol(proposal_cov,'lower');
%                     Z = randn(size(mean_gaussian_resampled,1), tile.ENSEMBLE.PARA.grid_ensemble_size ); %Z=randn(Np,Ne);
%                     value_gaussian_resampled = mean_gaussian_resampled + L*Z;
% 
%                     tile.ENSEMBLE.TEMP.value_gaussian(da.TEMP.pos_in_ensemble,:) = value_gaussian_resampled;
%                               
%                     da.TILE.ENSEMBLE = recalculate_ensemble_parameters_after_DA(da.TILE.ENSEMBLE, tile, da.PARA.ensemble_variables);
% 
%                     %reset timestamps, no need to reset timestamps in the
%                     %CG stratigraphy since the old states are read in
%                     tile.t = da.TEMP.last_assimilation_date;
% 
%                     da.TEMP.first_obs_index =[];
%                     da.TEMP.index_next_obs = [];
%                     da.TEMP.time_next_obs = [];
%                     for i=1:size(da.PARA.observation_files,1)
%                         if ~isempty(find(da.STATVAR.obs_time{i,1} > tile.t, 1))
%                             da.TEMP.first_obs_index = [da.TEMP.first_obs_index; find(da.STATVAR.obs_time{i,1} > tile.t, 1)];
%                             da.TEMP.index_next_obs = [da.TEMP.index_next_obs; da.TEMP.first_obs_index(i,1)]; %start with
%                             da.TEMP.time_next_obs = [da.TEMP.time_next_obs; da.STATVAR.obs_time{i,1}(da.TEMP.first_obs_index(i,1),1)];
%                         else
%                             da.TEMP.first_obs_index = [da.TEMP.first_obs_index; size(da.STATVAR.observations{i,1}, 1)];
%                             da.TEMP.index_next_obs = [da.TEMP.index_next_obs; size(da.STATVAR.observations{i,1}, 1)];
%                             da.TEMP.time_next_obs = [da.TEMP.time_next_obs; tile.FORCING.PARA.end_time + 1];
%                         end
%                     end
%                     
%                     da.DA_TIME = min(da.TEMP.time_next_obs);
%                     tile.OUT = reset_timestamp_out(tile.OUT,tile);
%                 end
%                 if da.PARA.recalculate_stratigraphy ==1
%                     da.TEMP.recalculate_stratigraphy_now = 1;
%                 end
%                 if da.PARA.recalculate_forcing ==1
%                     da.TEMP.recalculate_forcing_now = 1;
%                 end
            end

        end
        
        function da = save_state(da, tile)
            state = tile.SUBSURFACE_CLASS.STATVAR;
            
            save([tile.PARA.result_path da.TEMP.run_name '/saved_state.mat'], 'state');
            
        end
    end
end

