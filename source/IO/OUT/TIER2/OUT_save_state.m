%========================================================================
% CryoGrid OUT class OUT_last_timestep
% can be used to periodically store the entire CryoGrid stratigraphy, with 
% the ouput file getting overwritten each time. This is useful for being 
% able to recover and restart simulations (e.g. when the cluster is shut 
% down due to maintenances) in case of a long runtime. It is also posisble 
% to only store the final state of the CryoGrid stratigraphy, e.g. to run
% ensembles starting from the final state of an existing TILE class. Use
% TILE_BUILDER class "restart_OUT_last_timestep" to restart simulations.
% S. Westermann, October 2025
%========================================================================


classdef OUT_save_state < OUT_BASE

    properties

        STRATIGRAPHY
    
    end
    
    
    methods
        
        
        function out = provide_PARA(out)         

            out.PARA.save_timestamp = []; %assigned by optimization action classes, corresponds to start and end times of the TILE class
            out.PARA.out_folder = [];
            out.PARA.identifier = [];
            out.PARA.tag = [];
            out.PARA.tag2 = [];

        end

		
		function out = finalize_init(out, tile)
		
            if ~isempty(out.PARA.tag) && sum(isnan(out.PARA.tag))>0
                out.PARA.tag = [];
            end

            out.SAVE_TIME = out.PARA.save_timestamp(find(out.PARA.save_timestamp(:,1)-tile.PARA.start_time>=0, 1, 'first'), 1);
            out.OUTPUT_TIME = out.SAVE_TIME;
        end
        
        %-------time integration----------------
        
        function out = store_OUT(out, tile)

            %out_tag = out.PARA.tag;
            out_tag = [out.PARA.tag '_' out.PARA.tag2];
            if strcmp(out_tag(end), '_')
                out_tag = out_tag(1:end-1);
            end
            if ~isempty(out_tag) && strcmp(out_tag(1), '_')
                out_tag = out_tag(2:end);
            end
            
            if tile.t>=out.SAVE_TIME 
                                
                run_name = tile.PARA.run_name; 
                if isempty(out.PARA.out_folder) || sum(isnan(out.PARA.out_folder))>0
                    result_path = tile.PARA.result_path;
                else
                    result_path = out.PARA.out_folder;
                end
                out.STRATIGRAPHY = copy(tile);
                out.STRATIGRAPHY.RUN_INFO = [];
                out.STRATIGRAPHY.BUILDER = [];

                if ~(exist([result_path run_name])==7)
                    mkdir([result_path run_name])
                end
                if isempty(out_tag) || all(isnan(out_tag))
                    save([result_path run_name '/' run_name '_state.mat'], 'out')
                else
                    save([result_path run_name '/' run_name '_state_' out_tag '.mat'], 'out')
                end
                
                out.SAVE_TIME = out.PARA.save_timestamp(find(out.PARA.save_timestamp(:,1)-tile.t>0, 1, 'first'), 1);
                if isempty(out.SAVE_TIME)
                    out.SAVE_TIME = Inf;
                end
                out.OUTPUT_TIME = out.SAVE_TIME;
            end
            out.STRATIGRAPHY = [];
        end
               
                % %-------------param file generation-----
                % function out = param_file_info(out)
                %     out = provide_PARA(out);
                % 
                %     out.PARA.STATVAR = [];
                %     out.PARA.options = [];
                %     out.PARA.class_category = 'OUT';
                % 
                %     out.PARA.default_value.save_timestep = {''};
                %     out.PARA.comment.save_timestep = {'in days, if empty save final state at the end of the run, so that it can serve as initial condition for new runs'};
                % 
                %     out.PARA.default_value.tag = {''};
                %     out.PARA.comment.tag = {'additional tag added to file name'};
                % end

        
    end
end