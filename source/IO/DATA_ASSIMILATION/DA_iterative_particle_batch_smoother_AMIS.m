classdef DA_iterative_particle_batch_smoother_AMIS < matlab.mixin.Copyable

    
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
        end
        
        function da = provide_CONST(da)
            
        end
        
        function da = provide_STATVAR(da)
            
        end
        
        function da = finalize_init(da, tile)
            da.TILE = tile;
%             offset = length(num2str(da.TILE.PARA.worker_number))+1;
            da.TEMP.run_name = da.TILE.PARA.run_name;%(1, 1:length(da.TILE.PARA.run_name)-offset); %original results folder w.o worker number
            da.TEMP.first_obs_index =[];
            da.TEMP.index_next_obs = [];
            da.TEMP.time_next_obs = [];
            for i=1:size(da.PARA.observation_files,1)
                temp=load([da.PARA.observation_paths{i,1} da.PARA.observation_files{i,1}], 'OBS');
                da.STATVAR.obs_time{i,1} = temp.OBS.time;
                da.STATVAR.observations{i,1} = temp.OBS.observations;
                da.STATVAR.obs_variance{i,1} = temp.OBS.obs_variance;
                da.STATVAR.modeled_obs{i,1} = da.STATVAR.observations{i,1}.*NaN;
                da.ENSEMBLE.weights = repmat(1./da.TILE.PARA.ensemble_size, 1, da.TILE.PARA.ensemble_size); %set equal weights before 1st DA step
                
                da.OBS_OP{i,1} = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(da.PARA.observable_classes{i,1}){da.PARA.observable_classes_index(i,1)});
                da.OBS_OP{i,1} = finalize_init(da.OBS_OP{i,1}, tile);
                da.TEMP.first_obs_index = [da.TEMP.first_obs_index; find(da.STATVAR.obs_time{i,1} > tile.t, 1)];
                da.TEMP.index_next_obs = [da.TEMP.index_next_obs; da.TEMP.first_obs_index(i,1)]; %start with 
                da.TEMP.time_next_obs = [da.TEMP.time_next_obs; da.STATVAR.obs_time{i,1}(da.TEMP.first_obs_index(i,1),1)];
            end

            da.TEMP.num_iterations = 0;

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
                        
            %vector of positions to convert between the list of variables
            %that are changed by the DA to the full list of perturbed
            %variables in ENSEMBLE
            pos_in_ensemble = [];
            for i=1:size(da.PARA.ensemble_variables,1)
                 pos_in_ensemble = [pos_in_ensemble; find(strcmp(da.PARA.ensemble_variables{i,1}, tile.ENSEMBLE.TEMP.variable_name))];
            end
            da.TEMP.pos_in_ensemble = pos_in_ensemble;

            da.TEMP.recalculate_stratigraphy_now = 0;
            
            da.ENSEMBLE.modeled_obs = []; %Yp in Kris code, size N_obs x N_ens x N_iterations
            da.ENSEMBLE.value_gaussian = [];  %Xp in Kris code, size N_param x N_ens x N_iterations
            da.TEMP.cov_gaussian_resampled = []; %size N_param x N_param x N_iterations
            da.TEMP.mean_gaussian_resampled = []; %size N_param x N_iterations
            
        end
        
        
        function da = DA_step(da, tile)
            %save the state and the ensemble variables when the DA begins,
            %this step is redone every time the DA is happy and moves on in time
            if ~da.TEMP.assimilation_started && tile.t>=da.TEMP.last_assimilation_date
                da.TEMP.last_assimilation_date = tile.t;
                da.TEMP.assimilation_started = 1;
                da.ENSEMBLE.weights_old = da.ENSEMBLE.weights;
                da = save_state(da, tile); %save states at the start of the assimilation period, index is 0
                da.TEMP.old_value_gaussian = tile.ENSEMBLE.TEMP.value_gaussian(da.TEMP.pos_in_ensemble,1);
                da.TEMP.old_mean_gaussian = tile.ENSEMBLE.TEMP.mean_gaussian(da.TEMP.pos_in_ensemble,1);
                da.TEMP.old_std_gaussian = tile.ENSEMBLE.TEMP.std_gaussian(da.TEMP.pos_in_ensemble,1);
                
                %fill the necessary inforatio for first AMIS
                da.TEMP.cov_gaussian_resampled = diag(tile.ENSEMBLE.TEMP.std_gaussian(da.TEMP.pos_in_ensemble,1).^2);
                da.TEMP.mean_gaussian_resampled = tile.ENSEMBLE.TEMP.mean_gaussian(da.TEMP.pos_in_ensemble,1);
                
%                 labBarrier;
%                 value_gaussian = tile.ENSEMBLE.TEMP.value_gaussian(da.TEMP.pos_in_ensemble,1); 
%                 data_package = [];
%                 data_package = pack(da, data_package, 'value_gaussian', value_gaussian);
%                 
%                 da.ENSEMBLE.value_gaussian = cat(3, da.ENSEMBLE.value_gaussian, repmat(value_gaussian, 1, da.TILE.PARA.ensemble_size) .* NaN);
%                 da.ENSEMBLE.value_gaussian(:, da.TILE.PARA.worker_number, end) = value_gaussian;
%                 
%                 %send
%                 for i = 1:da.TILE.PARA.ensemble_size
%                     if i~=da.TILE.PARA.worker_number
%                         labSend(data_package, i, 1);
%                     end
%                 end
%                 
%                 for i = 1:da.TILE.PARA.ensemble_size
%                     if i~=da.TILE.PARA.worker_number
%                         data_package_in = labReceive(i, 1);
%                         if ~isempty(data_package_in)
%                             da = unpack(da, data_package_in, i); %read received column vector and transform into ENSEMBLE
%                         end
%                     end
%                 end
%                 labBarrier;
            end
            
            if da.TEMP.assimilation_started && tile.t>= da.DA_TIME
                %loop over all observation data sets
                for i=1:size(da.STATVAR.obs_time,1)
                    if tile.t>= da.TEMP.time_next_obs(i,1)
                        da.STATVAR.modeled_obs{i,1}(da.TEMP.index_next_obs(i,1),1) = observable_operator(da.OBS_OP{i,1}, tile);
                        %disp(da.STATVAR.modeled_obs{i,1}(da.TEMP.index_next_obs(i,1),1) )
                        if da.TEMP.index_next_obs(i,1) < size(da.STATVAR.observations{i,1}, 1) 
                            da.TEMP.index_next_obs(i,1) = da.TEMP.index_next_obs(i,1) + 1;
                            da.TEMP.time_next_obs(i,1) = da.STATVAR.obs_time{i,1}(da.TEMP.index_next_obs(i,1),1);
                        else %end of observations reached
                            da.TEMP.time_next_obs(i,1) = tile.FORCING.PARA.end_time +1;
                        end
                    end
                    
                end
                
                da.DA_TIME = min(da.TEMP.time_next_obs);
                
            end
            
            if tile.t>=da.DA_STEP_TIME
                labBarrier;
                %synchronize
                
                da.TEMP.num_iterations = da.TEMP.num_iterations + 1; %num_iertaions is now the value of the current iteration, i.e. the 

                data_package = [];
                
                modeled_obs = [];%gather modeled observations in one vector
                for i=1:size(da.STATVAR.modeled_obs,1) 
                    modeled_obs = [modeled_obs; da.STATVAR.modeled_obs{i,1}(da.TEMP.first_obs_index(i,1):da.TEMP.index_next_obs(i,1)-1,1)];  %ONLY USE THE PART IN THE OBS INTERVAL, ALSO MAKE A SIMILAR VECTOR FOR OBSERVATIONS AND VARIANCES
                end
                data_package = pack(da, data_package, 'modeled_obs', modeled_obs);
                da.ENSEMBLE.modeled_obs = cat(3, da.ENSEMBLE.modeled_obs, repmat(modeled_obs, 1, da.TILE.PARA.ensemble_size) .* NaN);
                da.ENSEMBLE.modeled_obs(:, da.TILE.PARA.worker_number, end) = modeled_obs;
                
                value_gaussian = tile.ENSEMBLE.TEMP.value_gaussian(da.TEMP.pos_in_ensemble,1); 
                data_package = pack(da, data_package, 'value_gaussian', value_gaussian);

                da.ENSEMBLE.value_gaussian = cat(3, da.ENSEMBLE.value_gaussian, repmat(value_gaussian, 1, da.TILE.PARA.ensemble_size) .* NaN);
                da.ENSEMBLE.value_gaussian(:, da.TILE.PARA.worker_number, end) = value_gaussian;
                
                %send
                for i = 1:da.TILE.PARA.ensemble_size
                    if i~=da.TILE.PARA.worker_number
                        labSend(data_package, i, 1);
                    end
                end
                
                for i = 1:da.TILE.PARA.ensemble_size
                    if i~=da.TILE.PARA.worker_number
                        data_package_in = labReceive(i, 1);
                        if ~isempty(data_package_in)
                            da = unpack(da, data_package_in, i); %read received column vector and transform into ENSEMBLE
                        end
                    end
                end
                
                labBarrier;
                %synchronize
                
%                 save(['save_all_1_' num2str(tile.PARA.worker_number) '.mat'])
               
                %actual DA 
%                 if da.TEMP.num_iterations == 1
%                     da = PBS(da);
%                 else
%                     da = AMIS(da);
%                 end


                da = AMIS(da);
                N_weights = numel(da.ENSEMBLE.weights(:)); %same as number of ensemble members x number of iterations
                    
                %store
                if strcmp(da.PARA.store_format, 'all') %&& da.TILE.PARA.worker_number==1
                    da_store = copy(da);
                    da_store.TILE = [];
                    if isempty(da.PARA.store_file_tag) || isnan(da.PARA.store_file_tag)
                        save([tile.PARA.result_path tile.PARA.run_name(1:end-2) '/' 'da_store_'  datestr(tile.t, 'yyyymmdd') '_' num2str(da.TEMP.num_iterations) '_' num2str(da.TILE.PARA.worker_number) '.mat'], 'da_store')
                    else
                        save([tile.PARA.result_path tile.PARA.run_name(1:end-2) '/' 'da_store_'  datestr(tile.t, 'yyyymmdd') '_' num2str(da.TEMP.num_iterations) '_' num2str(da.TILE.PARA.worker_number) '_' da.PARA.store_file_tag '.mat'], 'da_store')
                    end
                end
                %end store
                
                %start the iterative loop which either terminates or starts a new iteration, dependent on diversity of ensemmble
                if da.ENSEMBLE.effective_ensemble_size./tile.PARA.ensemble_size >= da.PARA.min_ensemble_diversity || da.TEMP.num_iterations>=da.PARA.max_iterations
                    %terminate and move on in time, i.e. do a normal resampling of model state and resample parameters according to the learning coefficient
                    
                    if da.ENSEMBLE.effective_ensemble_size./tile.PARA.ensemble_size >= da.PARA.min_ensemble_diversity
                        disp('successful, ensemble is sufficiently diverse')
                    else
                        disp('maximum number of iterations reached')
                    end
                    
                    %store
                    if strcmp(da.PARA.store_format, 'final') && da.TILE.PARA.worker_number==1
                        da_store = copy(da);
                        da_store.TILE = [];
                        if isempty(da.PARA.store_file_tag) || isnan(da.PARA.store_file_tag)
                            save([tile.PARA.result_path tile.PARA.run_name(1:end-2) '/' 'da_store_'  datestr(tile.t, 'yyyymmdd') '.mat'], 'da_store')
                        else
                            save([tile.PARA.result_path tile.PARA.run_name(1:end-2) '/' 'da_store_' datestr(tile.t, 'yyyymmdd') '_' da.PARA.store_file_tag '.mat'], 'da_store')
                        end
                    end
                    %end store
                    da = save_state(da, tile); %save the current run for resampling, index is num_iterations
                   
                    rng(da.DA_TIME.*da.TEMP.num_iterations+1);                 
                    resample_ID = randsample(da.TILE.PARA.ensemble_size.*da.TEMP.num_iterations, da.TILE.PARA.ensemble_size, true, da.ENSEMBLE.weights(:)); 
                    %find the correct ID of the suviving ensemble member
                    [ensemble_number, iteration_number] = meshgrid([1:da.TILE.PARA.ensemble_size], [1:da.TEMP.num_iterations]);
                    ensemble_number = ensemble_number';
                    iteration_number = iteration_number';
                                                            
                    %read the new stratigraphy and info from file
                    %IMPORTANT: this must now load from all iterations!!!  
                    temp=load([tile.PARA.result_path  da.TEMP.run_name '/tile_' ...
                        num2str(ensemble_number(resample_ID(tile.PARA.worker_number,1))) '_' num2str(iteration_number(resample_ID(tile.PARA.worker_number,1)))  '.mat']);
                    variables = fieldnames(temp.state);
                    for i=1:size(variables,1)
                        if ~isempty(temp.state.(variables{i,1}))
                            tile.(variables{i,1}) = temp.state.(variables{i,1});
                        end
                    end
                    
                    rand_sequence = rand(1,tile.PARA.ensemble_size);
                    if da.PARA.learning_coefficient >= rand_sequence(1, tile.PARA.worker_number)
                        %learning
                        value_gaussian_resampled = da.ENSEMBLE.value_gaussian; %Xp(:,:,sell);
                        value_gaussian_resampled=reshape(value_gaussian_resampled,[size(value_gaussian_resampled,1), da.TILE.PARA.ensemble_size.*da.TEMP.num_iterations]); % Concatenate across all weights (past proposals)
                        value_gaussian_resampled = value_gaussian_resampled(:,resample_ID);
                        %only resampling is done as the ensemble is assumed to not be degenerate
                        tile.ENSEMBLE.TEMP.value_gaussian(da.TEMP.pos_in_ensemble,1) = value_gaussian_resampled(:, tile.PARA.worker_number);
                    else
                        tile.ENSEMBLE.TEMP.value_gaussian(da.TEMP.pos_in_ensemble,1) = da.TEMP.old_value_gaussian;
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
                    
                    %reset modelled observations and vectors conatining parameters 
                    da.ENSEMBLE.modeled_obs = []; %Yp in Kris code, size N_obs x N_ens x N_iterations
                    da.ENSEMBLE.value_gaussian = [];  %Xp in Kris code, size N_param x N_ens x N_iterations
                    da.TEMP.cov_gaussian_resampled = []; %size N_param x N_param x N_iterations
                    da.TEMP.mean_gaussian_resampled = []; %size N_param x N_iterations
                    
                    %set new times
                    da.DA_TIME = min(da.TEMP.time_next_obs);
                    da.TEMP.last_assimilation_date = tile.t;
                    da.ENSEMBLE.weights_old = da.ENSEMBLE.weights;
                    da = save_state(da, tile); %save new states at the start of the new assimilation period, so that it can be read again if ensemble is degenerate
                    
                    %this is only valid for learning coefficient 0!
                    da.TEMP.old_value_gaussian = tile.ENSEMBLE.TEMP.value_gaussian(da.TEMP.pos_in_ensemble,1);
                    da.TEMP.cov_gaussian_resampled = diag(tile.ENSEMBLE.TEMP.std_gaussian(da.TEMP.pos_in_ensemble,1).^2);
                    da.TEMP.mean_gaussian_resampled = tile.ENSEMBLE.TEMP.mean_gaussian(da.TEMP.pos_in_ensemble,1);
                
                    
                    %assimilation successful, ensemble is not degenerate, move on in time
                    da.TEMP.num_iterations = 0; %reset number of iterations for the next DA period
                    da = save_state(da, tile); %saves state 0 for next asssimilation period 
                else
                    %iterate and go back to start of DA period, resample and inflate again, start over with "old" model states at the beginning of the DA period
%                     save(['save_all' num2str(tile.PARA.worker_number) '.mat'])
                    disp('ensemble degenerate, one more time')

                    da = save_state(da, tile); %save the state with the current iteration
                                       
                    %load "old" state, i.e. from 1st iteration
                    temp=load([tile.PARA.result_path da.TEMP.run_name '/tile_' num2str(da.TILE.PARA.worker_number) '_0.mat']);
                    variables = fieldnames(temp.state);
                    for i=1:size(variables,1)
                        if ~isempty(temp.state.(variables{i,1}))
                            tile.(variables{i,1}) = temp.state.(variables{i,1});
                        end
                    end
                    
                    rng(da.DA_TIME.*da.TEMP.num_iterations+2);                 
                    resample_ID = randsample(da.TILE.PARA.ensemble_size.*da.TEMP.num_iterations, da.TILE.PARA.ensemble_size, true, da.ENSEMBLE.weights(:)); 
                    value_gaussian_resampled = da.ENSEMBLE.value_gaussian; %Xp(:,:,sell);
                    thetapc=value_gaussian_resampled; % For clipping potentially
                    value_gaussian_resampled=reshape(value_gaussian_resampled,[size(value_gaussian_resampled,1), da.TILE.PARA.ensemble_size.*da.TEMP.num_iterations]); % Concatenate across all weights (past proposals)
                    value_gaussian_resampled = value_gaussian_resampled(:,resample_ID);
                    mean_gaussian_resampled = mean(value_gaussian_resampled,2); % proposal mean for next iteration
                    A=value_gaussian_resampled-mean_gaussian_resampled;
                    cov_gaussian_resampled=(1./da.TILE.PARA.ensemble_size).*(A*A'); % proposal covariance for next iteration
                    d = min(da.ENSEMBLE.effective_ensemble_size./da.TILE.PARA.ensemble_size./da.PARA.min_ensemble_diversity,1); %d=min(diversity/adapt_thresh,1);

                    clip = round(da.PARA.min_ensemble_diversity.*da.TILE.PARA.ensemble_size); %clip=round(adapt_thresh*Ne);
                    if sum(da.ENSEMBLE.weights(:) > 1/(10*da.TILE.PARA.ensemble_size))>clip  %should have a threshold, is 0 in Kris original code
                        disp('clipping')
                        w = da.ENSEMBLE.weights(:);
                        ws=sort(w,'descend');
                        wc=ws(clip);
                        w(w>wc)=wc; % > truncated ("clipped") < anti-truncated
                        w=w./sum(w); %renomralize
                        resamplec_ID = randsample(da.TILE.PARA.ensemble_size.*da.TEMP.num_iterations, da.TILE.PARA.ensemble_size, true, da.ENSEMBLE.weights(:)); 
                        thetapc=reshape(thetapc,[size(thetapc,1), da.TILE.PARA.ensemble_size.*da.TEMP.num_iterations]); %reshape(thetapc,[Np,Nw]);
                        thetapc=thetapc(:,resamplec_ID);
                        pmc=mean(value_gaussian_resampled,2); % proposal mean for next iteration
                        Ac=thetapc-pmc;
                        pcc=(1./da.TILE.PARA.ensemble_size).*(Ac*Ac');
                        mean_gaussian_resampled = pmc;
                        cov_gaussian_resampled=pcc;
                    else
                        disp('no clipping')
                        pric = diag(da.TEMP.old_std_gaussian.^2);
                        cov_gaussian_resampled = d.*cov_gaussian_resampled+(1-d).*pric;
                    end

                    da.TEMP.cov_gaussian_resampled = cat(3, da.TEMP.cov_gaussian_resampled, cov_gaussian_resampled); %propc(:,:,ell)=pc;
                    da.TEMP.mean_gaussian_resampled = cat(3, da.TEMP.mean_gaussian_resampled, mean_gaussian_resampled); %propm(:,ell)=pm; %SEB: here, propm and the others are expaned by one

                    rng(da.DA_TIME.*da.TEMP.num_iterations+3);
                    Z = randn(size(mean_gaussian_resampled,1), tile.PARA.ensemble_size); %Z=randn(Np,Ne);
                    L=chol(cov_gaussian_resampled,'lower');
                    value_gaussian_resampled= mean_gaussian_resampled+L*Z;
    %                da.ENSEMBLE.value_gaussian = cat(3,
    %                da.ENSEMBLE.value_gaussian, value_gaussian_resampled);
    %                %Xp(:,:,ell)=thetap; %SEB: Xp is the last argument in
    %                AMIS %SEB: this step is done at the beginning of the
    %                next iteration!
                    

%                     value_gaussian_resampled = da.ENSEMBLE.value_gaussian(:,resample_ID);  %    thetap=thetap(:,resample);
%                     mean_gaussian_resampled = mean(value_gaussian_resampled,2);   %propm=mean(thetap,2);
%                     
%                     deviation = value_gaussian_resampled - mean_gaussian_resampled; % A=thetap-propm;
%                     proposal_cov =  (1./tile.PARA.ensemble_size) .* (deviation*deviation');  %   propc=(1/Ne).*(A*A'); % Proposal covariance
%                     % Inflate proposal covariance to avoid overconfidence
%                     
%                     proposal_cov = proposal_cov + 0.1 .* (1 - da.ENSEMBLE.effective_ensemble_size./ tile.PARA.ensemble_size).* diag(da.TEMP.old_std_gaussian); % propc=propc+0.1.*(1-diversity).*pric;
%                     % Gaussian resampling using the Cholesky decomposition
%                     L=chol(proposal_cov,'lower');
%                     Z = randn(size(mean_gaussian_resampled,1), tile.PARA.ensemble_size); %Z=randn(Np,Ne);
%                     value_gaussian_resampled = mean_gaussian_resampled + L*Z;

                    tile.ENSEMBLE.TEMP.value_gaussian(da.TEMP.pos_in_ensemble,1) = value_gaussian_resampled(:, tile.PARA.worker_number);
                              
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
                   % tile.OUT = reset_timestamp_out(tile.OUT,tile);
                end
                if da.PARA.recalculate_stratigraphy ==1
                    da.TEMP.recalculate_stratigraphy_now = 1;
                end
            end
            
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
            
            
        
        function da = AMIS(da)
            %[w,Neff]= AMIS( Ypred, yobs, R, prim, pricov, propm, propc, props )
            
            % Adaptive Multiple Importance Sampling (AMIS)
            % N.B. The observation errors are assumed to be uncorrelated (diagonal R)
            % and Gaussian. This can easily be changed.
            %
            % Dimensions: No = Number of observations in the batch to assimilate.
            %             Np = Number of parameters to update.
            %             Ne = Number of ensemble members (particles) per iteration.
            %             Nl = Number of iterations so far in AMIS
            % -----------------------------------------------------------------------
            % Inputs:
            %
            %
            % Ypred   => No x Ne x Nl matrix containing an ensemble of Ne predicted
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
            % pricov => Np x Np Prior covariance matrix
            %
            % propm => Np x Nl Mean of the 'deterministic mixture' (DM) of proposals
            %
            % propc => Np x Np x Nl Covariance of the DM of proposals
            %
            % props => Np x Ne x Nl Samples from DM of proposals
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
            
            %proposal = da.ENSEMBLE.value_gaussian;
            prim = da.TEMP.old_mean_gaussian;
            pricov = diag(da.TEMP.old_std_gaussian.^2);
            %these must be made new!
            propm = da.TEMP.mean_gaussian_resampled; %mean of DM
            propc = da.TEMP.cov_gaussian_resampled; %cov of DM
            props = da.ENSEMBLE.value_gaussian; %da.ENSEMBLE.samples_deterministic_mixture_proposals; %full 3D matrix containing all the sampled parameters
                        
            % Neglog of the target term, the unnormalized posterior, up to constants
            phi=zeros(size(da.ENSEMBLE.modeled_obs,2), size(da.ENSEMBLE.modeled_obs,3)); %Ne,Nl); % Group by iterations
            % Log-sum-exp of the "DM" of proposals including normalizing constants
            lsepsi=zeros(size(da.ENSEMBLE.modeled_obs,2), size(da.ENSEMBLE.modeled_obs,3)); %Ne,Nl);
            
            for ell=1:size(da.ENSEMBLE.modeled_obs,3)
                sampell=props(:,:,ell); % Samples from proposal ell "sampell" (not a typo)
                A0ell=sampell-prim;
                phi0ell=0.5.*(A0ell')*(pricov\A0ell);
                phi0ell=diag(phi0ell);
                phi0ell=phi0ell';
                Ypell=da.ENSEMBLE.modeled_obs(:,:,ell);
                residuell=da.ENSEMBLE.observations - Ypell;
                phidell=0.5*Rinv*(residuell.^2);
                phiell=phi0ell+phidell;
                phi(:,ell)=phiell;
                
                % "DM" of proposals, the denominator term
                psij=zeros(size(da.ENSEMBLE.modeled_obs,2),size(da.ENSEMBLE.modeled_obs,3));
                
                for j=1:size(da.ENSEMBLE.modeled_obs,3)
                    mj=propm(:,j);
                    Cj=propc(:,:,j);
                    cj=det(2*pi*Cj).^(-1/2); % Normalizing constant of proposal j
                    lcj=log(cj);
                    Aj=sampell-mj;
                    psi=0.5*(Aj')*(Cj\Aj);
                    psi=diag(psi);
                    psi=psi';
                    psi=psi-lcj; % Must include the constant in this case.
                    psij(:,j)=psi;
                end
                lsepsiell=logsumexp(da, -psij,2);
                lsepsi(:,ell)=lsepsiell;
            end
        
            logw=-phi-lsepsi;
            logw=logw - max(logw(:));
            da.ENSEMBLE.weights = exp(logw)./sum(exp(logw(:))); %dimension N_ens x N_iterations
            da.ENSEMBLE.effective_ensemble_size = 1./sum(da.ENSEMBLE.weights(:).^2);
            
        end
        
        function lse=logsumexp(da, a,dim)
                % A simple implementation of the "log-sum-exp" function that is
                % valid for matrices.
                amax=max(a,[],dim);
                ad=a-amax;
                lse=amax+log(sum(exp(ad),dim));
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
               da.ENSEMBLE.(variable_name)(:,received_from_worker, da.TEMP.num_iterations) = data_package(i+1:i+data_package(i,1),1);
               i = i + data_package(i,1) + 1;
            end
        end
        
        function da = save_state(da, tile)
            state = copy(tile);
            variables = fieldnames(state);
            for i=1:size(variables,1)
                if ~strcmp(variables{i,1}, 'LATERAL') && ~strcmp(variables{i,1}, 'TOP') && ~strcmp(variables{i,1}, 'BOTTOM') && ~strcmp(variables{i,1}, 'TOP_CLASS') && ~strcmp(variables{i,1}, 'BOTTOM_CLASS') && ~strcmp(variables{i,1}, 'OUT')
                    state.(variables{i,1}) = [];
                end
            end
            %each ensemble member also stores its OUT field, so in this
            %case OUT is with each state and then is passed on, so that it gets stored directly after DA step -> for this to work, the OUT save_interval must be the same or larrger as the DA interval 
            %consider doing the same for the other DA class!
            save([tile.PARA.result_path da.TEMP.run_name '/tile_' num2str(da.TILE.PARA.worker_number) '_' num2str(da.TEMP.num_iterations) '.mat'], 'state');
            
            labBarrier;

        end

    end
end

