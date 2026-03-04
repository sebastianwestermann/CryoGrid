%========================================================================
% CryoGrid OUT class OUT_snow delivering simple snow depths
% S. Westermann, Sep 2025
%========================================================================


classdef OUT_snow < matlab.mixin.Copyable
 

    properties
        MISC
        TEMP
        PARA
        OUTPUT_TIME
        SAVE_TIME
        CONST
        result
        timestamp
        
    end
    
    
    methods
        
        %initialization
        
        function out = provide_PARA(out)         

            out.PARA.output_timestep = [];
            out.PARA.save_date = [];
            out.PARA.save_interval = [];
            out.PARA.tag = [];
        end

        
        function out = provide_CONST(out)

        end

        
        function out = provide_STATVAR(out)

        end

        
        function out = finalize_init(out, tile)

            forcing = tile.FORCING;
            
            % Set the next (first) output time. This is the next (first) time output
            % is collected (in memory) for later storage to disk.
            out.OUTPUT_TIME = forcing.PARA.start_time + out.PARA.output_timestep;
            out.timestamp = [];
    
            out.result.snow_depth = [];

            % Set the next (first) save time. This is the next (first) time all the
            % collected output is saved to disk.
            if isempty(out.PARA.save_interval) || isnan(out.PARA.save_interval) 
                out.SAVE_TIME = forcing.PARA.end_time;
            else
                out.SAVE_TIME = min(forcing.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(forcing.PARA.start_time,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
            end
                        
        end
        
        %---------------time integration-------------
                    
        function out = store_OUT(out, tile)           
            
            t = tile.t;
            TOP = tile.TOP; 
            BOTTOM = tile.BOTTOM;
            forcing = tile.FORCING;
            run_name = tile.PARA.run_name; %tile.RUN_NUMBER;
            result_path = tile.PARA.result_path;            
            timestep = tile.timestep;
            out_tag = out.PARA.tag;
            
            if t>=out.OUTPUT_TIME
                % It is time to collect output
                % Store the current state of the model in the out structure.

                out.timestamp = [out.timestamp t];
                
                snow_depth = 0;
                CURRENT = TOP.NEXT;
                if isprop(CURRENT, 'CHILD') && CURRENT.CHILD ~= 0
                    class_name = class(CURRENT.CHILD);
                    if strcmp(class_name(1:4), 'SNOW')
                        snow_depth = CURRENT.CHILD.STATVAR.area .* CURRENT.CHILD.STATVAR.layerThick ./ CURRENT.STATVAR.area(1,1);
                    end
                end
                class_name = class(CURRENT);
                if strcmp(class_name(1:4), 'SNOW')
                    snow_depth = sum(CURRENT.STATVAR.layerThick);
                end                

               out.result.snow_depth = [out.result.snow_depth ; snow_depth];

                
                % Set the next OUTPUT_TIME
                out.OUTPUT_TIME = min(out.SAVE_TIME, out.OUTPUT_TIME + out.PARA.output_timestep);
                
                if t>=out.SAVE_TIME
                    % It is time to save all the collected model output to disk
                     
                    if ~(exist([result_path run_name])==7)
                        mkdir([result_path run_name])
                    end
                    CG_out = out.result;
                    CG_out.timestamp = out.timestamp;
                    if isempty(out_tag) || all(isnan(out_tag))
                        save([result_path run_name '/' run_name '_snow_' datestr(t,'yyyymmdd') '.mat'], 'CG_out')
                    else
                        save([result_path run_name '/' run_name '_snow_' out_tag '_' datestr(t,'yyyymmdd') '.mat'], 'CG_out')
                    end
                    
                    % Clear the out structure
                    out.timestamp = [];
                    out.result.snow_depth = [];
                    
                    if ~isnan(out.PARA.save_interval)
                        % If save_interval is defined, uptate SAVE_TIME for next save opertion 
                        out.SAVE_TIME = min(forcing.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(out.SAVE_TIME,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
                        % If save_interval is not defined, we will save at the very end of the model run
					end
                end
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