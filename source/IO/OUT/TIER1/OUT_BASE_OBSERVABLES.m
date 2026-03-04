%========================================================================

% S. Westermann, Feb 2026
%========================================================================


classdef OUT_BASE_OBSERVABLES < OUT_BASE 
 

    properties

    end
    
    
    methods


        function out = finalize_init(out, tile)

            if ~isempty(out.PARA.tag) && sum(isnan(out.PARA.tag))>0
                out.PARA.tag = [];
            end
            out = set_out_save_time_init(out, tile);
                        
        end
        
        function out = set_out_save_time_init(out, tile)

            if ~isempty(out.PARA.timestamps) %DA mode, defined timestamps written by DA
                timestamps = [out.PARA.timestamps; Inf];
                out.OUTPUT_TIME = timestamps(find(timestamps(:,1)-tile.PARA.start_time > 0, 1, 'first'), 1);
                %out.OUTPUT_TIME = out.PARA.timestamps(find(out.PARA.timestamps(:,1)-tile.PARA.start_time > 0, 1, 'first'), 1);
                out.SAVE_TIME = Inf;
                out.TEMP.write_out_init = 0;
                out.TEMP.write_out_final = 0;
            else
                out.OUTPUT_TIME = tile.PARA.start_time + out.PARA.output_timestep;
                if isempty(out.PARA.save_interval) || sum(isnan(out.PARA.save_interval)) > 0 %write at the very end of the simulations
                    out.SAVE_TIME = tile.FORCING.PARA.end_time;
                    out.TEMP.write_out_init = 0;
                    out.TEMP.write_out_final = 0;
                else
                    if out.PARA.save_interval == 0 %write for each DA
                        out.SAVE_TIME = Inf;
                        out.TEMP.write_out_init = double(out.PARA.write_file_mode==1);
                        out.TEMP.write_out_final = double(out.PARA.write_file_mode==2);
                    else %write for DA after a certain interval
                        out.SAVE_TIME = min(tile.FORCING.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(tile.PARA.start_time,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
                        out.TEMP.write_out_init = 0;
                        out.TEMP.write_out_final = 0;
                    end
                end
            end
            if isempty(out.OUTPUT_TIME)
                out.OUTPUT_TIME = Inf;
            end

        end

        function out = set_output_time(out, tile)

            if ~isempty(out.PARA.timestamps)
                timestamps = [out.PARA.timestamps; Inf];
                out.OUTPUT_TIME = timestamps(find(timestamps(:,1)-tile.t > 0, 1, 'first'), 1);
                %out.OUTPUT_TIME = out.PARA.timestamps(find(out.PARA.timestamps(:,1)-tile.t>0, 1, 'first'), 1);
            else
                out.OUTPUT_TIME = min(out.SAVE_TIME, out.OUTPUT_TIME + out.PARA.output_timestep);
            end

        end

        function out = set_output_time_after_write(out, tile)

            if ~isempty(out.PARA.timestamps)
                timestamps = [out.PARA.timestamps; Inf];
                out.OUTPUT_TIME = timestamps(find(timestamps(:,1)-tile.t > 0, 1, 'first'), 1);
                %out.OUTPUT_TIME = out.PARA.timestamps(find(out.PARA.timestamps(:,1)-tile.t>0, 1, 'first'), 1);
                out.TEMP.write_out_init = 0;
                out.TEMP.write_out_final = 0;
            else
                out.OUTPUT_TIME = min(out.SAVE_TIME, out.OUTPUT_TIME + out.PARA.output_timestep);
                if ~isempty(out.PARA.save_interval) && sum(isnan(out.PARA.save_interval))==0 && out.PARA.save_interval>0
                    %reset write_out variables to be ready for next period
                    out.TEMP.write_out_init = 0;
                    out.TEMP.write_out_final = 0;
                end
            end
            if isempty(out.OUTPUT_TIME)
                out.OUTPUT_TIME = Inf;
            end
        end

        function out = out2file_CG(out, tile)

            out = out2file_CG@OUT_BASE(out, tile);
            out = set_output_time_after_write(out, tile);

        end

        function out = write_out_final(out, tile) %make these functions in the OPT version of all "normal" OUT classes
            if out.TEMP.write_out_final == 1
                out = out2file_CG(out, tile);
            end
        end

        function out = write_out_init(out, tile)
            if out.TEMP.write_out_init == 1
                out = out2file_CG(out, tile);
            end
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