%IO functions for DA for MULTITILEs


classdef DA_IO_MULTITILE < matlab.mixin.Copyable

    
    properties

    end
    
    methods

        function da = make_scratchfolder(da_IO, da, tile)
            da.PARA.scratchfolder = tile.PARA.result_path;
        end

        function res = get_size_of_modeled_obs(da_IO, da, tile)
            res = tile.PARA.ensemble_size;
        end


        function da = save_state(da_IO, da, tile)
           %  state = copy(tile);
           %  variables = fieldnames(state);
           %  for i=1:size(variables,1)
           %      if ~strcmp(variables{i,1}, 'SUBSURFACE_CLASS') && ~strcmp(variables{i,1}, 'LATERAL') && ~strcmp(variables{i,1}, 'TOP') && ~strcmp(variables{i,1}, 'BOTTOM') && ~strcmp(variables{i,1}, 'TOP_CLASS') && ~strcmp(variables{i,1}, 'BOTTOM_CLASS') && ~strcmp(variables{i,1}, 'OUT')
           %          state.(variables{i,1}) = [];
           %      end
           %  end
           % 
           % save([da.PARA.scratchfolder da.TEMP.run_name '/tile_' num2str(da.TEMP.num_iterations) '.mat'], 'state');
           % 
        end

         function da = save_prior_state(da_IO, da, tile)
             % save iteration 1 aditionally as prior of each year's DA, so
             % this is not overwritten
             % if da.TEMP.num_iterations == 1
             % 
             %     state = copy(tile);
             %     variables = fieldnames(state);
             %     for i=1:size(variables,1)
             %         if ~strcmp(variables{i,1}, 'SUBSURFACE_CLASS') && ~strcmp(variables{i,1}, 'LATERAL') && ~strcmp(variables{i,1}, 'TOP') && ~strcmp(variables{i,1}, 'BOTTOM') && ~strcmp(variables{i,1}, 'TOP_CLASS') && ~strcmp(variables{i,1}, 'BOTTOM_CLASS') && ~strcmp(variables{i,1}, 'OUT')
             %             state.(variables{i,1}) = [];
             %         end
             %     end
             %     save([tile.PARA.result_path da.TEMP.run_name '/tile_' num2str(da.TILE.PARA.worker_number) '_prior_' datestr(tile.t, 'yyyy') '.mat'], 'state');
             % end
         end

         function da = collect_modeled_observations(da_IO, da, tile)

             modeled_obs = [];
             for i=1:size(da.STATVAR.modeled_obs,1)

                 modeled_obs = [modeled_obs; da.STATVAR.modeled_obs{i,1}(da.TEMP.first_obs_index(i,1):da.TEMP.index_next_obs(i,1)-1,:)];  %ONLY USE THE PART IN THE OBS INTERVAL, ALSO MAKE A SIMILAR VECTOR FOR OBSERVATIONS AND VARIANCES

                 % da.STATVAR.modeled_obs{i,1}(1:da.TEMP.index_next_obs(i,1)-1,:) = []; %delete everything until then
                 % da.STATVAR.obs_time{i,1}(1:da.TEMP.index_next_obs(i,1)-1,:) = [];
                 % da.STATVAR.observations{i,1}(1:da.TEMP.index_next_obs(i,1)-1,:) = [];
                 % da.STATVAR.obs_variance{i,1}(1:da.TEMP.index_next_obs(i,1)-1,:) = [];
                 %
                 % da.TEMP.first_obs_index(i,1) = da.TEMP.index_next_obs(i,1);

             end

             da.ENSEMBLE.modeled_obs = cat(3,  da.ENSEMBLE.modeled_obs, modeled_obs);
             da.ENSEMBLE.value_gaussian = cat(3, da.ENSEMBLE.value_gaussian, tile.ENSEMBLE.TEMP.value_gaussian);
         
         end


         function da = save_da_results_all(da_IO, da, tile)
             if strcmp(da.PARA.store_format, 'all') 
                 da_store = copy(da);
                 da_store.TILE = [];
                 if isempty(da.PARA.store_file_tag) || isnan(da.PARA.store_file_tag)
                     save([tile.PARA.result_path tile.PARA.run_name(1:end-2) '/' 'da_store_'  datestr(tile.t, 'yyyymmdd') '_' num2str(da.TEMP.num_iterations) '.mat'], 'da_store')
                 else
                     save([tile.PARA.result_path tile.PARA.run_name(1:end-2) '/' 'da_store_'  datestr(tile.t, 'yyyymmdd') '_' num2str(da.TEMP.num_iterations) '_' da.PARA.store_file_tag '.mat'], 'da_store')
                 end
             end
         end

         function da = save_da_results_final(da_IO, da, tile)
             if strcmp(da.PARA.store_format, 'final') 
                 da_store = copy(da);
                 da_store.TILE = [];
                 if isempty(da.PARA.store_file_tag) || isnan(da.PARA.store_file_tag)
                     save([tile.PARA.result_path tile.PARA.run_name(1:end-2) '/' 'da_store_'  datestr(tile.t, 'yyyymmdd') '_' num2str(da.TEMP.num_iterations) '.mat'], 'da_store')
                 else
                     save([tile.PARA.result_path tile.PARA.run_name(1:end-2) '/' 'da_store_'  datestr(tile.t, 'yyyymmdd') '_' num2str(da.TEMP.num_iterations) '_' da.PARA.store_file_tag '.mat'], 'da_store')
                 end
             end
         end
         
        % function da = save_state(da_IO, da, tile)
        %     state = tile.SUBSURFACE_CLASS.STATVAR;
        % 
        %     save([tile.PARA.result_path da.TEMP.run_name '/saved_state.mat'], 'state');
        % 
        % end
    end

end

