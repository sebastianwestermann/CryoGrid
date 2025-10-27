%IO functions for DA for normal TILEs, including communication between
%workers


classdef DA_IO_TILE < matlab.mixin.Copyable

    
    properties

    end
    
    methods

        function da = collect_modeled_observations(da_IO, da, tile)

            spmdBarrier;
            %synchronize

            data_package = [];

            modeled_obs = [];%gather modeled observations in one vector
            for i=1:size(da.STATVAR.modeled_obs,1)
                if isinf(da.TEMP.time_next_obs)
                    modeled_obs = [modeled_obs; da.STATVAR.modeled_obs{i,1}(da.TEMP.first_obs_index(i,1):da.TEMP.index_next_obs(i,1),1)];  %ONLY USE THE PART IN THE OBS INTERVAL, ALSO MAKE A SIMILAR VECTOR FOR OBSERVATIONS AND VARIANCES
                else
                    modeled_obs = [modeled_obs; da.STATVAR.modeled_obs{i,1}(da.TEMP.first_obs_index(i,1):da.TEMP.index_next_obs(i,1)-1,1)];  %ONLY USE THE PART IN THE OBS INTERVAL, ALSO MAKE A SIMILAR VECTOR FOR OBSERVATIONS AND VARIANCE
                end
            end
            data_package = pack(da_IO, data_package, 'modeled_obs', modeled_obs);
            da.ENSEMBLE.modeled_obs = cat(3, da.ENSEMBLE.modeled_obs, repmat(modeled_obs, 1, da.TILE.PARA.ensemble_size) .* NaN);
            da.ENSEMBLE.modeled_obs(:, da.TILE.PARA.worker_number, end) = modeled_obs;

            value_gaussian = tile.ENSEMBLE.TEMP.value_gaussian(da.TEMP.pos_in_ensemble,1);
            data_package = pack(da_IO, data_package, 'value_gaussian', value_gaussian);

            da.ENSEMBLE.value_gaussian = cat(3, da.ENSEMBLE.value_gaussian, repmat(value_gaussian, 1, da.TILE.PARA.ensemble_size) .* NaN);
            da.ENSEMBLE.value_gaussian(:, da.TILE.PARA.worker_number, end) = value_gaussian;

            %send
            for i = 1:da.TILE.PARA.ensemble_size
                if i~=da.TILE.PARA.worker_number
                    spmdSend(data_package, i, 1);
                end
            end

            for i = 1:da.TILE.PARA.ensemble_size
                if i~=da.TILE.PARA.worker_number
                    data_package_in = spmdReceive(i, 1);
                    % if tile.PARA.worker_number == 3 && i==5
                    %     save('test.mat', 'data_package')
                    % end
                    if ~isempty(data_package_in)
                        da = unpack(da_IO, da, data_package_in, i); %read received column vector and transform into ENSEMBLE
                    end
                end
            end
            % if tile.PARA.worker_number == 3
            %     save('test3.mat', 'da')
            % end
            spmdBarrier;
        end


        function da = resample_state(da_IO, da, tile, resample_ID)

            %find the correct ID of the suviving ensemble member
            [ensemble_number, iteration_number] = meshgrid([1:da.TILE.PARA.ensemble_size], [1:da.TEMP.num_iterations]);
            ensemble_number = ensemble_number';
            iteration_number = iteration_number';

            %read the new stratigraphy and info from file
            %IMPORTANT: this must now load from all iterations!!!
            % each worker draws one of the resampled members and
            % loads the state
            for ii=1:tile.PARA.ensemble_size
                if ii==tile.PARA.worker_number
                    temp=load([da.PARA.scratchfolder  da.TEMP.run_name '/tile_' ...
                        num2str(ensemble_number(resample_ID(tile.PARA.worker_number,1))) '_' num2str(iteration_number(resample_ID(tile.PARA.worker_number,1)))  '.mat']);
                end
                spmdBarrier;
            end

            variables = fieldnames(temp.state);
            for i=1:size(variables,1)
                if ~isempty(temp.state.(variables{i,1}))
                    da.TILE.(variables{i,1}) = temp.state.(variables{i,1});
                end
            end
        end

        
        function da = load_prior_state(da_IO, da, tile)
            temp=load([da.PARA.scratchfolder da.TEMP.run_name '/tile_' num2str(da.TILE.PARA.worker_number) '_0.mat']);
            variables = fieldnames(temp.state);
            for i=1:size(variables,1)
                if ~isempty(temp.state.(variables{i,1}))
                    da.TILE.(variables{i,1}) = temp.state.(variables{i,1});
                end
            end

        end


        function da = make_scratchfolder(da_IO, da, tile)
            if isempty(da.PARA.scratchfolder) || sum(isnan(da.PARA.scratchfolder))>0
                da.PARA.scratchfolder = tile.PARA.result_path;
            else
                if ~exist([da.PARA.scratchfolder da.TEMP.run_name], 'dir')
                    mkdir([da.PARA.scratchfolder da.TEMP.run_name]);
                end
            end
        end

        function res = get_size_of_modeled_obs(da_IO, da, tile)
            res = 1;
        end


        function data_package = pack(da_IO, data_package, var_name, var) %transform into column vector ready to send
                %variables{i,1}
                data_package=[data_package; size(var_name,2); double(var_name)']; % # of characters followed by characters as doubles
                data_package=[data_package; size(var,1); var]; % # of entries followed by values
        end
        

        function da = unpack(da_IO, da, data_package, received_from_worker) %read received column vector and transform into STATVAR
            i=1;
            while i<=size(data_package,1)
               variable_name = char(data_package(i+1:i+data_package(i,1),1)');
               i = i + data_package(i,1) + 1;
               da.ENSEMBLE.(variable_name)(:,received_from_worker, da.TEMP.num_iterations) = data_package(i+1:i+data_package(i,1),1);
               i = i + data_package(i,1) + 1;
            end
        end


         function da = save_state(da_IO, da, tile)
            state = copy(tile);
            variables = fieldnames(state);
            for i=1:size(variables,1)
                if ~strcmp(variables{i,1}, 'LATERAL') && ~strcmp(variables{i,1}, 'TOP') && ~strcmp(variables{i,1}, 'BOTTOM') && ~strcmp(variables{i,1}, 'TOP_CLASS') && ~strcmp(variables{i,1}, 'BOTTOM_CLASS') && ~strcmp(variables{i,1}, 'OUT')
                    state.(variables{i,1}) = [];
                end
            end

           save([da.PARA.scratchfolder da.TEMP.run_name '/tile_' num2str(da.TILE.PARA.worker_number) '_' num2str(da.TEMP.num_iterations) '.mat'], 'state');
           spmdBarrier

         end

         function da = save_prior_state(da_IO, da, tile)
          % save iteration 1 aditionally as prior of each year's DA, so
          % this is not overwritten
          if da.TEMP.num_iterations == 1

              state = copy(tile);
              variables = fieldnames(state);
              for i=1:size(variables,1)
                  if ~strcmp(variables{i,1}, 'LATERAL') && ~strcmp(variables{i,1}, 'TOP') && ~strcmp(variables{i,1}, 'BOTTOM') && ~strcmp(variables{i,1}, 'TOP_CLASS') && ~strcmp(variables{i,1}, 'BOTTOM_CLASS') && ~strcmp(variables{i,1}, 'OUT')
                      state.(variables{i,1}) = [];
                  end
              end
              save([tile.PARA.result_path da.TEMP.run_name '/tile_' num2str(da.TILE.PARA.worker_number) '_prior_' datestr(tile.t, 'yyyy') '.mat'], 'state');
              spmdBarrier
          end

         end

         function da = save_da_results_all(da_IO, da, tile)
             if da.TILE.PARA.worker_number==3 && strcmp(da.PARA.store_format, 'all')
                 da_store = copy(da);
                 da_store.TILE = [];
                 if isempty(da.PARA.store_file_tag) || isnan(da.PARA.store_file_tag)
                     %  save([tile.PARA.result_path tile.PARA.run_name(1:end-2) '/' 'da_store_'  datestr(tile.t, 'yyyymmdd') '_' num2str(da.TEMP.num_iterations) '_' num2str(da.TILE.PARA.worker_number) '.mat'], 'da_store')
                     save([tile.PARA.result_path tile.PARA.run_name(1:end-2) '/' 'da_store_'  datestr(tile.t, 'yyyymmdd') '_' num2str(da.TEMP.num_iterations) '.mat'], 'da_store')
                 else
                     %  save([tile.PARA.result_path tile.PARA.run_name(1:end-2) '/' 'da_store_'  datestr(tile.t, 'yyyymmdd') '_' num2str(da.TEMP.num_iterations) '_' num2str(da.TILE.PARA.worker_number) '_' da.PARA.store_file_tag '.mat'], 'da_store')
                     save([tile.PARA.result_path tile.PARA.run_name(1:end-2) '/' 'da_store_'  datestr(tile.t, 'yyyymmdd') '_' num2str(da.TEMP.num_iterations) '_' da.PARA.store_file_tag '.mat'], 'da_store')
                 end
             end
         end

         function da = save_da_results_final(da_IO, da, tile)
             if da.TILE.PARA.worker_number==1 && strcmp(da.PARA.store_format, 'final')
                 da_store = copy(da);
                 da_store.TILE = [];
                 if isempty(da.PARA.store_file_tag) || isnan(da.PARA.store_file_tag)
                     save([tile.PARA.result_path tile.PARA.run_name(1:end-2) '/' 'da_store_'  datestr(tile.t, 'yyyymmdd') '.mat'], 'da_store')
                 else
                     save([tile.PARA.result_path tile.PARA.run_name(1:end-2) '/' 'da_store_' datestr(tile.t, 'yyyymmdd') '_' da.PARA.store_file_tag '.mat'], 'da_store')
                 end
             end
         end


    end

end

