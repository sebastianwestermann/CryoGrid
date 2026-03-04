%========================================================================

% S. Westermann, Feb 2026
%========================================================================


classdef OUT_BASE < matlab.mixin.Copyable
 

    properties
        TEMP
        PARA
        CONST
        OUTPUT_TIME
        SAVE_TIME
    end
    
    
    methods
        
        function out = provide_PARA(out)         

        end

        
        function out = provide_CONST(out)

        end

        
        function out = provide_STATVAR(out)

        end

        
        function out = finalize_init(out, tile)

            if ~isempty(out.PARA.tag) && sum(isnan(out.PARA.tag))>0
                out.PARA.tag = [];
            end

            % Set the next (first) output time. This is the next (first) time output
            % is collected (in memory) for later storage to disk.
            out.OUTPUT_TIME = tile.FORCING.PARA.start_time + out.PARA.output_timestep;

            % Set the next (first) save time. This is the next (first) time all the
            % collected output is saved to disk.
            if isempty(out.PARA.save_interval) || sum(isnan(out.PARA.save_interval))>0 
                out.SAVE_TIME = tile.FORCING.PARA.end_time;
            else
                out.SAVE_TIME = min(tile.FORCING.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(tile.FORCING.PARA.start_time,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
            end
                        
        end
        

        function out = state2out_all(out, tile)

        end

        function out = out2file_all(out, tile)

            out.TEMP.tag = ['_' out.PARA.tag '_' out.PARA.tag2 '_'];
            out.TEMP.tag = strrep(out.TEMP.tag, '___', '_');
            out.TEMP.tag = strrep(out.TEMP.tag, '__', '_');
            out.TEMP.identifier = tile.RUN_INFO.PPROVIDER.PARA.identifier;

            if ~(exist([tile.PARA.result_path tile.PARA.run_name])==7)
                mkdir([tile.PARA.result_path tile.PARA.run_name])
            end

            out.TEMP.identifier = tile.RUN_INFO.PPROVIDER.PARA.identifier;
            save([tile.PARA.result_path tile.PARA.run_name '/' tile.PARA.run_name out.TEMP.tag datestr(tile.t,'yyyymmdd') '.mat'], 'out')

            out.STRATIGRAPHY=[];
            if isfield(out, 'LATERAL')
                out.LATERAL=[];
            end
            out.TIMESTAMP=[];
            
            % if ~isnan(out.PARA.save_interval)
            %     % If save_interval is defined, uptate SAVE_TIME for next save opertion
            %     out.SAVE_TIME = min(forcing.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(out.SAVE_TIME,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
            %     % If save_interval is not defined, we will save at the very end of the model run
            %     % and thus do not need to update SAVE_TIME (update would fail because save_interval is nan)
            % end

            if ~isempty(out.PARA.save_interval) && sum(isnan(out.PARA.save_interval))==0 && out.PARA.save_interval > 0
                out.SAVE_TIME = min(tile.FORCING.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(out.SAVE_TIME,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
            end
        end

        function out = out2file_CG(out, tile)

            out.TEMP.tag = ['_' out.PARA.tag '_' out.PARA.tag2 '_'];
            out.TEMP.tag = strrep(out.TEMP.tag, '___', '_');
            out.TEMP.tag = strrep(out.TEMP.tag, '__', '_');

            if ~(exist([tile.PARA.result_path tile.PARA.run_name])==7)
                mkdir([tile.PARA.result_path tile.PARA.run_name])
            end

            fn = fieldnames(out.STATVAR);
            for i=1:size(fn,1)
                a.(out.TEMP.output_var).(fn{i,1}) = out.STATVAR.(fn{i,1});
            end
            % CG_out.timestamp = out.TIMESTAMP';
            a.(out.TEMP.output_var).identifier = tile.RUN_INFO.PPROVIDER.PARA.identifier;
            save([tile.PARA.result_path tile.PARA.run_name '/' tile.PARA.run_name '_' out.TEMP.keyword out.TEMP.tag datestr(tile.t,'yyyymmdd') '.mat'], '-struct', 'a', 'CG_*')

            % Clear the out structure
            for i=1:size(fn,1)
                if iscell(out.STATVAR.(fn{i,1}))
                    out.STATVAR.(fn{i,1}) = {};
                else
                    out.STATVAR.(fn{i,1}) = [];
                end
            end

            if ~isempty(out.PARA.save_interval) && sum(isnan(out.PARA.save_interval))==0 && out.PARA.save_interval > 0
                out.SAVE_TIME = min(tile.FORCING.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(out.SAVE_TIME,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
            end
        end

        function out = write_out_init(out, tile)

        end

        function out = write_out_final(out, tile)

        end
        
        % 
        % %-------------param file generation-----
        % function out = param_file_info(out)
        %     out = provide_PARA(out);
        % 
        %     out.PARA.STATVAR = [];
        %     out.PARA.options = [];
        %     out.PARA.class_category = 'OUT';
        % 
        %     out.PARA.comment.variables = {'select output variables: T, water, ice, waterIce, XwaterIce, Xwater, Xice, saltConc are supported'};
        %     out.PARA.options.variables.name = 'H_LIST';
        %     out.PARA.options.variables.entries_x = {'T' 'water' 'ice' 'waterIce'};
        % 
        %     out.PARA.comment.upper_elevation = {'upper elevation of output domain, in m a.s.l.; must match altitude in TILE!'};
        % 
        %     out.PARA.comment.lower_elevation = {'lower elevation of output domain, in m a.s.l.; must match altitude in TILE!'};
        % 
        %     out.PARA.default_value.target_grid_size = {0.02};
        %     out.PARA.comment.target_grid_size = {'cell size of regular grid that values are interpolated to'};
        % 
        %     out.PARA.default_value.output_timestep = {0.25};
        %     out.PARA.comment.output_timestep = {'timestep of output [days]'};
        % 
        %     out.PARA.default_value.save_date = {'01.09.'};
        %     out.PARA.comment.save_date = {'date (dd.mm.) when output file is written'};
        % 
        %     out.PARA.default_value.save_interval = {1};
        %     out.PARA.comment.save_interval = {'interval of output files [years]'};
        % 
        %     out.PARA.default_value.tag = {''};
        %     out.PARA.comment.tag = {'additional tag added to file name'};
        % end
        

    end
end