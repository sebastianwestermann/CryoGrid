classdef DA_FUNCTIONS < matlab.mixin.Copyable

    
    properties
        TILE
        OBS_OP
        OBSERVATIONS
        PARA
        CONST
        STATVAR
        TEMP
        DA_TIME     
        DA_STEP_TIME
        ENSEMBLE
        IO
    end
    
    methods

        function da = load_observable_operators(da, tile)
            for i=1:size(da.PARA.observable_classes,1)
                da.OBS_OP{i,1} = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(da.PARA.observable_classes{i,1}){da.PARA.observable_classes_index(i,1)});
                da.OBS_OP{i,1} = finalize_init(da.OBS_OP{i,1}, tile);
            end
        end

        function da = load_observations(da, tile)

            da.TEMP.first_obs_index =[];
            da.TEMP.index_next_obs = [];
            da.TEMP.time_next_obs = [];
            for i=1:size(da.PARA.observation_classes,1)

                da.OBSERVATIONS{i,1} = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(da.PARA.observation_classes{i,1}){da.PARA.observation_classes_index(i,1)});
                da.OBSERVATIONS{i,1}  = finalize_init(da.OBSERVATIONS{i,1}, tile);
                da.OBSERVATIONS{i,1}  = read_observations(da.OBSERVATIONS{i,1}, tile);

                da.STATVAR.obs_time{i,1} = da.OBSERVATIONS{i,1}.STATVAR.time;
                da.STATVAR.observations{i,1} = da.OBSERVATIONS{i,1}.STATVAR.observations;
                da.STATVAR.obs_variance{i,1} = da.OBSERVATIONS{i,1}.STATVAR.obs_variance;
                da.STATVAR.modeled_obs{i,1} = repmat(da.STATVAR.observations{i,1}.*NaN,1, get_size_of_modeled_obs(da.IO, da, tile));
                da.OBSERVATIONS{i,1}.STATVAR = [];

                %delete observations that are no longer relevant for DA
                delete_obs = find(da.STATVAR.obs_time{i,1}(:,1) < da.TEMP.last_assimilation_date);
                da.STATVAR.obs_time{i,1}(delete_obs,:) = [];
                da.STATVAR.observations{i,1}(delete_obs,:) = [];
                da.STATVAR.obs_variance{i,1}(delete_obs,:) = [];
                da.STATVAR.modeled_obs{i,1}(delete_obs,:) = [];

                da.TEMP.first_obs_index = [da.TEMP.first_obs_index; find(da.STATVAR.obs_time{i,1} > tile.t, 1)];
                da.TEMP.index_next_obs = [da.TEMP.index_next_obs; da.TEMP.first_obs_index(i,1)]; %start with
                da.TEMP.time_next_obs = [da.TEMP.time_next_obs; da.STATVAR.obs_time{i,1}(da.TEMP.first_obs_index(i,1),1)];

            end
            da.DA_TIME = min(da.TEMP.time_next_obs);
           
        end

        function da = reset_observation_time_indeces(da, tile)


            %Check this against old code
            % da.ENSEMBLE.modeled_obs = [];
            % da.ENSEMBLE.value_gaussian = [];
            %end check this

            da.TEMP.first_obs_index =[];
            da.TEMP.index_next_obs = [];
            da.TEMP.time_next_obs = [];
            for i=1:size(da.PARA.observation_classes,1)

                %delete observations that are no longer relevant for DA
                delete_obs = find(da.STATVAR.obs_time{i,1}(:,1) < da.TEMP.last_assimilation_date);
                da.STATVAR.obs_time{i,1}(delete_obs,:) = [];
                da.STATVAR.observations{i,1}(delete_obs,:) = [];
                da.STATVAR.obs_variance{i,1}(delete_obs,:) = [];
                da.STATVAR.modeled_obs{i,1}(delete_obs,:) = [];

                if ~isempty(find(da.STATVAR.obs_time{i,1} > tile.t, 1))
                    da.TEMP.first_obs_index = [da.TEMP.first_obs_index; find(da.STATVAR.obs_time{i,1} > tile.t, 1)];
                    da.TEMP.index_next_obs = [da.TEMP.index_next_obs; da.TEMP.first_obs_index(i,1)]; %start with
                    da.TEMP.time_next_obs = [da.TEMP.time_next_obs; da.STATVAR.obs_time{i,1}(da.TEMP.first_obs_index(i,1),1)];
                else
                    da.TEMP.first_obs_index = [da.TEMP.first_obs_index; size(da.STATVAR.observations{i,1}, 1)];
                    da.TEMP.index_next_obs = [da.TEMP.index_next_obs; size(da.STATVAR.observations{i,1}, 1)];
                    da.TEMP.time_next_obs = [da.TEMP.time_next_obs; Inf];
                end

            end
            da.DA_TIME = min(da.TEMP.time_next_obs);
        end
        


        function da = get_DA_step_time(da, tile)
            if strcmp(da.PARA.assimilation_frequency, 'year')
                current_year = str2num(datestr(tile.t, 'yyyy'));
                da.DA_STEP_TIME = datenum([da.PARA.assimilation_date num2str(current_year + da.PARA.assimilation_interval)], 'dd.mm.yyyy');
            elseif strcmp(da.PARA.assimilation_frequency, 'month')
                current_year = str2num(datestr(tile.t, 'yyyy'));
                current_month = str2num(datestr(tile.t, 'mm'));
                da.DA_STEP_TIME = datenum(current_year, current_month + da.PARA.assimilation_interval, da.PARA.assimilation_date);
            elseif strcmp(da.PARA.assimilation_frequency, 'day')
                da.DA_STEP_TIME = tile.t + da.PARA.assimilation_interval;
            elseif strcmp(da.PARA.assimilation_frequency, 'next_obs')
                da.DA_STEP_TIME = da.DA_TIME; 
            end
        end

        function da = init_DA(da,tile, dec)
            if dec || isempty(da.PARA.start_assimilation_period) || sum(isnan(da.PARA.start_assimilation_period))>0
                da.TEMP.last_assimilation_date = tile.t; %start with initial state
                da.TEMP.assimilation_started = 1;
                da = load_observations(da, tile);
                da = get_DA_step_time(da, tile);
                da = save_state(da.IO,da,tile);
            else
                da.TEMP.assimilation_started = 0;
                da.TEMP.last_assimilation_date = datenum(da.PARA.start_assimilation_period(1,1), da.PARA.start_assimilation_period(2,1), da.PARA.start_assimilation_period(3,1));
            end
        end

        function da = get_ensemble_info(da, tile)
            %vector of positions to convert between the list of variables
            %that are changed by the DA to the full list of perturbed 
            %variables in ENSEMBLE
            pos_in_ensemble = [];
            for i=1:size(da.PARA.ensemble_variables,1)
                 pos_in_ensemble = [pos_in_ensemble; find(strcmp(da.PARA.ensemble_variables{i,1}, tile.ENSEMBLE.TEMP.variable_name))];
            end
            da.TEMP.pos_in_ensemble = pos_in_ensemble;

        end


        function da = get_modeled_observations(da, tile)
            for i=1:size(da.STATVAR.obs_time,1)
                if tile.t>= da.TEMP.time_next_obs(i,1)
                    da.STATVAR.modeled_obs{i,1}(da.TEMP.index_next_obs(i,1),:) = observable_operator(da.OBS_OP{i,1}, tile);
                    % disp('collecting synthetic observations')
                    if da.TEMP.index_next_obs(i,1) < size(da.STATVAR.observations{i,1}, 1) %end of observations not  reached
                        da.TEMP.index_next_obs(i,1) = da.TEMP.index_next_obs(i,1) + 1;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        1;
                        da.TEMP.time_next_obs(i,1) = da.STATVAR.obs_time{i,1}(da.TEMP.index_next_obs(i,1),1);
                    else
                        %load new observations
                        da.OBSERVATIONS{i,1}  = read_observations(da.OBSERVATIONS{i,1}, tile);
                        if ~isempty(da.OBSERVATIONS{i,1}.STATVAR)
                            da.STATVAR.obs_time{i,1} = [da.STATVAR.obs_time{i,1}; da.OBSERVATIONS{i,1}.STATVAR.time];
                            da.STATVAR.observations{i,1} = [da.STATVAR.observations{i,1}; da.OBSERVATIONS{i,1}.STATVAR.observations];
                            da.STATVAR.obs_variance{i,1} = [da.STATVAR.obs_variance{i,1}; da.OBSERVATIONS{i,1}.STATVAR.obs_variance];
                            da.STATVAR.modeled_obs{i,1} = [da.STATVAR.modeled_obs{i,1}; repmat(da.OBSERVATIONS{i,1}.STATVAR.observations.*NaN,1, get_size_of_modeled_obs(da.IO,da, tile))];

                            da.TEMP.index_next_obs(i,1) = min(da.TEMP.index_next_obs(i,1) + 1, size(da.STATVAR.obs_time{i,1},1));
                            da.TEMP.time_next_obs(i,1) = da.STATVAR.obs_time{i,1}(da.TEMP.index_next_obs(i,1),1);
                            if tile.t > da.TEMP.time_next_obs(i,1)
                                da.TEMP.time_next_obs(i,1) = Inf;
                            end
                        else
                            da.TEMP.time_next_obs(i,1) = Inf;
                        end
                        da.OBSERVATIONS{i,1}.STATVAR = [];
                    end
                end
            end
            da.DA_TIME = min(da.TEMP.time_next_obs);
        end


        function da = collect_observations(da, tile)
            observations = [];
            obs_variance = [];
            for i=1:size(da.STATVAR.observations,1)
                if isinf(da.TEMP.time_next_obs)
                    observations = [observations; da.STATVAR.observations{i,1}(da.TEMP.first_obs_index(i,1):da.TEMP.index_next_obs(i,1),:)];
                    obs_variance = [obs_variance; da.STATVAR.obs_variance{i,1}(da.TEMP.first_obs_index(i,1):da.TEMP.index_next_obs(i,1),:)];
                else
                    observations = [observations; da.STATVAR.observations{i,1}(da.TEMP.first_obs_index(i,1):da.TEMP.index_next_obs(i,1)-1,:)];
                    obs_variance = [obs_variance; da.STATVAR.obs_variance{i,1}(da.TEMP.first_obs_index(i,1):da.TEMP.index_next_obs(i,1)-1,:)];
                end
            end
            
            da.ENSEMBLE.observations = observations;
            da.ENSEMBLE.obs_variance = obs_variance;
        end

        function da = resample_state_and_parameters(da, tile)
            rng(da.DA_STEP_TIME.*da.TEMP.num_iterations+1);
            % randsample N_ensemblesize values from all members (of
            % all iterations) weighted by their weight
            resample_ID = randsample(da.TILE.PARA.ensemble_size.*da.TEMP.num_iterations, da.TILE.PARA.ensemble_size, true, da.ENSEMBLE.weights(:));
           
            da = resample_state(da.IO, da, tile, resample_ID);
            %find the correct ID of the suviving ensemble member
            % [ensemble_number, iteration_number] = meshgrid([1:da.TILE.PARA.ensemble_size], [1:da.TEMP.num_iterations]);
            % ensemble_number = ensemble_number';
            % iteration_number = iteration_number';
            % 
            % %read the new stratigraphy and info from file
            % %IMPORTANT: this must now load from all iterations!!!
            % % each worker draws one of the resampled members and
            % % loads the state
            % temp=load([da.PARA.scratchfolder  da.TEMP.run_name '/tile_' ...
            %     num2str(ensemble_number(resample_ID(tile.PARA.worker_number,1))) '_' num2str(iteration_number(resample_ID(tile.PARA.worker_number,1)))  '.mat']);
            % variables = fieldnames(temp.state);
            % for i=1:size(variables,1)
            %     if ~isempty(temp.state.(variables{i,1}))
            %         da.TILE.(variables{i,1}) = temp.state.(variables{i,1});
            %     end
            % end

            % rand creates pseudorandom values drawn from standard
            % uniform distribution on the open interval (0,1), i.e.
            % with learning_coefficient = 0 it never fulfills the
            % if-criterion to learn
            rand_sequence = rand(1,tile.PARA.ensemble_size);
            if da.PARA.learning_coefficient >= rand_sequence(1, tile.PARA.worker_number)
                %learning
                value_gaussian_resampled = da.ENSEMBLE.value_gaussian; %Xp(:,:,sell);
                value_gaussian_resampled=reshape(value_gaussian_resampled,[size(value_gaussian_resampled,1), da.TILE.PARA.ensemble_size.*da.TEMP.num_iterations]); % Concatenate across all weights (past proposals)
                value_gaussian_resampled = value_gaussian_resampled(:,resample_ID);
                %only resampling is done as the ensemble is assumed to not be degenerate
                da.TILE.ENSEMBLE.TEMP.value_gaussian(da.TEMP.pos_in_ensemble,1) = value_gaussian_resampled(:, tile.PARA.worker_number);
            else
                da.TILE.ENSEMBLE.TEMP.value_gaussian(da.TEMP.pos_in_ensemble,1) = da.TEMP.old_value_gaussian;
            end
        end

        function da = resample_AMIS_old(da, tile)
            rng(da.DA_STEP_TIME.*da.TEMP.num_iterations+2); % seeds the random number generator new for each iteration
            resample_ID = randsample(da.TILE.PARA.ensemble_size.*da.TEMP.num_iterations, da.TILE.PARA.ensemble_size, true, da.ENSEMBLE.weights(:));
            value_gaussian_resampled = da.ENSEMBLE.value_gaussian; %Xp(:,:,sell);
            thetapc=value_gaussian_resampled; % For clipping potentially
            value_gaussian_resampled=reshape(value_gaussian_resampled,[size(value_gaussian_resampled,1), da.TILE.PARA.ensemble_size.*da.TEMP.num_iterations]); % Concatenate across all weights (past proposals)
            value_gaussian_resampled = value_gaussian_resampled(:,resample_ID);
            mean_gaussian_resampled = mean(value_gaussian_resampled,2); % proposal mean for next iteration
            A=value_gaussian_resampled-mean_gaussian_resampled;
            cov_gaussian_resampled=(1./da.TILE.PARA.ensemble_size).*(A*A'); % proposal covariance for next iteration
            d = min(da.ENSEMBLE.effective_ensemble_size./da.TILE.PARA.ensemble_size./da.PARA.min_ensemble_diversity,1); %d=min(diversity/adapt_thresh,1);

            a = rand(1);
            clip = round(da.PARA.min_ensemble_diversity.*da.TILE.PARA.ensemble_size); %clip=round(adapt_thresh*Ne);
            if sum(da.ENSEMBLE.weights(:) > 1/(10*da.TILE.PARA.ensemble_size))>clip  %should have a threshold, is 0 in Kris original code
                disp('clipping')
                w = da.ENSEMBLE.weights(:);
                ws=sort(w,'descend');
                wc=ws(clip);
                w(w>wc)=wc; % > truncated ("clipped") < anti-truncated
                w=w./sum(w); %renomralize
                resamplec_ID = randsample(da.TILE.PARA.ensemble_size.*da.TEMP.num_iterations, da.TILE.PARA.ensemble_size, true, w);
                thetapc=reshape(thetapc,[size(thetapc,1), da.TILE.PARA.ensemble_size.*da.TEMP.num_iterations]); %reshape(thetapc,[Np,Nw]);
                thetapc=thetapc(:,resamplec_ID);
                pmc=mean(value_gaussian_resampled,2); % proposal mean for next iteration
                Ac=thetapc-pmc;
                pcc=(1./da.TILE.PARA.ensemble_size).*(Ac*Ac');
                mean_gaussian_resampled = pmc;
                if all(eig(pcc)>0) %%sum(pcc(:))>0
                    cov_gaussian_resampled=pcc;
                else
                    disp('clipping unsuccessful')
                    pric = diag(da.TEMP.old_std_gaussian.^2);
                    % % 07.11.
                    % cov_gaussian_resampled = d.*cov_gaussian_resampled+(1-d).*((a.*pric)+(1-a).*(d)^(1/2*(da.TEMP.num_iterations-1)).*pric);
                    % % 07.11.
                    cov_gaussian_resampled = d.*cov_gaussian_resampled + (1-d).*(a + (1-a)./da.TEMP.num_iterations) .*pric;
                end
            else
                disp('no clipping')
                pric = diag(da.TEMP.old_std_gaussian.^2);
                % % 07.11.
                % cov_gaussian_resampled = d.*cov_gaussian_resampled+(1-d).*((a.*pric)+(1-a).*(d)^(1/2*(da.TEMP.num_iterations-1)).*pric);
                % % 07.11.
                cov_gaussian_resampled = d.*cov_gaussian_resampled + (1-d).*(a + (1-a)./da.TEMP.num_iterations) .*pric;
            end

            da.TEMP.cov_gaussian_resampled = cat(3, da.TEMP.cov_gaussian_resampled, cov_gaussian_resampled); %propc(:,:,ell)=pc;
            da.TEMP.mean_gaussian_resampled = cat(3, da.TEMP.mean_gaussian_resampled, mean_gaussian_resampled); %propm(:,ell)=pm; %SEB: here, propm and the others are expaned by one

            rng(da.DA_STEP_TIME.*da.TEMP.num_iterations+3);
            Z = randn(size(mean_gaussian_resampled,1), tile.PARA.ensemble_size); %Z=randn(Np,Ne);
            L=chol(cov_gaussian_resampled,'lower');
            value_gaussian_resampled= mean_gaussian_resampled+L*Z;
            da.TEMP.value_gaussian_resampled = value_gaussian_resampled; %this needs to be the same as da.ENSEMBLE.value_gaussian
            da.TILE.ENSEMBLE.TEMP.value_gaussian(da.TEMP.pos_in_ensemble,1) = value_gaussian_resampled(:, tile.PARA.worker_number);
        end
            

        function da = resample_AMIS(da, tile)
            rng(da.DA_STEP_TIME.*da.TEMP.num_iterations+2); % seeds the random number generator new for each iteration
            resample_ID = randsample(da.TILE.PARA.ensemble_size.*da.TEMP.num_iterations, da.TILE.PARA.ensemble_size, true, da.ENSEMBLE.weights(:));
            value_gaussian_resampled = da.ENSEMBLE.value_gaussian; %Xp(:,:,sell);
            thetapc=value_gaussian_resampled; % For clipping potentially
            value_gaussian_resampled=reshape(value_gaussian_resampled,[size(value_gaussian_resampled,1), da.TILE.PARA.ensemble_size.*da.TEMP.num_iterations]); % Concatenate across all weights (past proposals)
            value_gaussian_resampled = value_gaussian_resampled(:,resample_ID);
            mean_gaussian_resampled = mean(value_gaussian_resampled,2); % proposal mean for next iteration
            A=value_gaussian_resampled-mean_gaussian_resampled;
            cov_gaussian_resampled=(1./da.TILE.PARA.ensemble_size).*(A*A'); % proposal covariance for next iteration

            d = min(2.*(da.ENSEMBLE.effective_ensemble_size-1)./da.TILE.PARA.ensemble_size./da.PARA.min_ensemble_diversity,1); %d=min(diversity/adapt_thresh,1);

            a = rand(1);
            %clip = round(da.PARA.min_ensemble_diversity.*da.TILE.PARA.ensemble_size); %clip=round(adapt_thresh*Ne);
            clip = sum(double(da.ENSEMBLE.weights(:) > 1e-3));
            if clip>=2 %da.ENSEMBLE.effective_ensemble_size>2 %sum(da.ENSEMBLE.weights(:) > 1/(10*da.TILE.PARA.ensemble_size))>clip  %should have a threshold, is 0 in Kris original code
                disp('clipping')
                w = da.ENSEMBLE.weights(:);
                ws=sort(w,'descend');
                wc=ws(max(clip, 4));
                w(w>wc)=wc; % > truncated ("clipped") < anti-truncated
                w=w./sum(w); %renomralize
                resamplec_ID = randsample(da.TILE.PARA.ensemble_size.*da.TEMP.num_iterations, da.TILE.PARA.ensemble_size, true, w);
                thetapc=reshape(thetapc,[size(thetapc,1), da.TILE.PARA.ensemble_size.*da.TEMP.num_iterations]); %reshape(thetapc,[Np,Nw]);
                thetapc=thetapc(:,resamplec_ID);
                pmc=mean(value_gaussian_resampled,2); % proposal mean for next iteration
                Ac=thetapc-pmc;
                pcc=(1./da.TILE.PARA.ensemble_size).*(Ac*Ac');
                mean_gaussian_resampled = pmc;
                if all(eig(pcc)>0) && clip >= da.ENSEMBLE.effective_ensemble_size 
                    cov_gaussian_resampled=pcc;
                else
                    disp('clipping unsuccessful')
                    pric = diag(da.TEMP.old_std_gaussian.^2);
                    cov_gaussian_resampled = d.*pcc + (1-d).*(a + (1-a)./da.TEMP.num_iterations) .*pric;
                end
            else
                disp('no clipping')
                pric = diag(da.TEMP.old_std_gaussian.^2);
                cov_gaussian_resampled = d.*cov_gaussian_resampled + (1-d).*(a + (1-a)./da.TEMP.num_iterations) .*pric;
            end

            da.TEMP.cov_gaussian_resampled = cat(3, da.TEMP.cov_gaussian_resampled, cov_gaussian_resampled); %propc(:,:,ell)=pc;
            da.TEMP.mean_gaussian_resampled = cat(3, da.TEMP.mean_gaussian_resampled, mean_gaussian_resampled); %propm(:,ell)=pm; %SEB: here, propm and the others are expaned by one

            Z = randn(size(mean_gaussian_resampled,1), tile.PARA.ensemble_size); %Z=randn(Np,Ne);
            L=chol(cov_gaussian_resampled,'lower');
            value_gaussian_resampled= mean_gaussian_resampled+L*Z;
            da.TEMP.value_gaussian_resampled = value_gaussian_resampled; %this needs to be the same as da.ENSEMBLE.value_gaussian
            da.TILE.ENSEMBLE.TEMP.value_gaussian(da.TEMP.pos_in_ensemble,1) = value_gaussian_resampled(:, tile.PARA.worker_number);
        end

%-------------------------------------------------------------------------
% actual DA functions

        function da = PBS(da, tile)
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
            da.ENSEMBLE.weights = [];
            da.ENSEMBLE.effective_ensemble_size = [];
            
            for gridcell=1:tile.PARA.number_of_realizations
%                 modeled_obs = [];%gather modeled observations in one vector
%                 for i=1:size(da.STATVAR.modeled_obs,1)
%                     modeled_obs = [modeled_obs; da.STATVAR.modeled_obs{i,1}(:,gridcell:tile.PARA.number_of_realizations:size(da.STATVAR.modeled_obs{i,1},2))];  %ONLY USE THE PART IN THE OBS INTERVAL, ALSO MAKE A SIMILAR VECTOR FOR OBSERVATIONS AND VARIANCES
%                 end
                modeled_obs = da.ENSEMBLE.modeled_obs(:,gridcell:tile.PARA.number_of_realizations:size(da.ENSEMBLE.modeled_obs,2));
                observations = da.ENSEMBLE.observations(:,gridcell);
                obs_variance = da.ENSEMBLE.obs_variance(:,gridcell);
                
                valid = ~isnan(observations);
                modeled_obs = modeled_obs(valid,:);
                observations = observations(valid,:);             
                obs_variance = obs_variance(valid,:);
                if ~isempty(observations)
                    
                    No=size(observations,1);
                    Rinv=(obs_variance').^(-1);
                    
                    % Calculate the likelihood.
                    Inn=repmat(observations,1,size(modeled_obs,2))-modeled_obs;   % Innovation.
                    EObj=Rinv*(Inn.^2);                     % [1 x Ne] ensemble objective function.
                    LLH=-0.5.*EObj; % log-likelihoods.
                    normc=logsumexp(da, LLH,2);
                    
                    
                    % NB! The likelihood coefficient (1/sqrt(2*pi...)) is
                    % omitted because it drops out in the normalization
                    % of the likelihood. Including it (very small term) would lead
                    % to problems with FP division.
                    
                    % Calculate the posterior weights as the normalized likelihood.
                    logw = LLH-normc;
                    weights = exp(logw); % Posterior weights.
                    da.ENSEMBLE.weights = [da.ENSEMBLE.weights; weights];
                    
                    % Need "log-sum-exp" trick to overcome numerical issues for small R/large
                    % number of obs.
                    da.ENSEMBLE.effective_ensemble_size = [da.ENSEMBLE.effective_ensemble_size; 1./sum(weights.^2)];
                else
                    da.ENSEMBLE.weights = [da.ENSEMBLE.weights; repmat(NaN, 1, tile.ENSEMBLE.PARA.grid_ensemble_size)];
                    da.ENSEMBLE.effective_ensemble_size = [da.ENSEMBLE.effective_ensemble_size; NaN];
                end
            end
        end
            
        function [weights, effective_ensemble_size] = PBS2(da, observations, modeled_obs, obs_variance)
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

            No=size(observations,1);
            Rinv=(obs_variance').^(-1);

            % Calculate the likelihood.
            Inn=repmat(observations,1,size(modeled_obs,2))-modeled_obs;   % Innovation.
            EObj=Rinv*(Inn.^2);                     % [1 x Ne] ensemble objective function.
            LLH=-0.5.*EObj; % log-likelihoods.
            normc=logsumexp(da, LLH,2);


            % NB! The likelihood coefficient (1/sqrt(2*pi...)) is
            % omitted because it drops out in the normalization
            % of the likelihood. Including it (very small term) would lead
            % to problems with FP division.

            % Calculate the posterior weights as the normalized likelihood.
            logw = LLH-normc;
            weights = exp(logw); % Posterior weights.

            % Need "log-sum-exp" trick to overcome numerical issues for small R/large
            % number of obs.
            effective_ensemble_size = 1./sum(weights.^2);
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
        
        function da = adaptive_PBS(da, tile)
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
            propcov=(1./tile.ENSEMBLE.PARA.grid_ensemble_size ).*(A*A');
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
            
            
            No=size(da.ENSEMBLE.observations,1);
            Rinv=(da.ENSEMBLE.obs_variance').^(-1);
            
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

    end

end

