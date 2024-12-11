%========================================================================
% CryoGrid OUT class OUT_all
% CryoGrid OUT class defining storage format of the output 
% OUT_all stores identical copies of all GROUND classses (including STATVAR, TEMP, PARA) in the
% stratigraphy for each output timestep, while lateral classes are not stored.
% The user can specify the save date and the save interval (e.g. yearly
% files), as well as the output timestep (e.g. 6 hourly). The output files
% are in Matlab (".mat") format.
% S. Westermann, T. Ingeman-Nielsen, J. Scheer, June 2021
%========================================================================


classdef OUT_scalars < matlab.mixin.Copyable
 

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

            out.PARA.variables = [];
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
            %make struct "result" and initialize all variables defined by
            %the user as empty arrays
            out.timestamp = [];
            for i =  1:size(out.PARA.variables,1)
               out.result.(out.PARA.variables{i,1}) = []; 
            end

            % Set the next (first) save time. This is the next (first) time all the
            % collected output is saved to disk.
            if isempty(out.PARA.save_interval) || isnan(out.PARA.save_interval) 
                out.SAVE_TIME = forcing.PARA.end_time;
            else
                out.SAVE_TIME = min(forcing.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(forcing.PARA.start_time,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
            end
                        
        end
        
        %---------------time integration-------------
        
%         function out = store_OUT(out, t, TOP, BOTTOM, forcing, run_number, timestep, result_path)
            
        function out = store_OUT(out, tile)           
            
            t = tile.t;
            TOP = tile.TOP; 
            BOTTOM = tile.BOTTOM;
            forcing = tile.FORCING;
            run_name = tile.PARA.run_name; %tile.RUN_NUMBER;
            result_path = tile.PARA.result_path;            
            timestep = tile.timestep;
            out_tag = out.PARA.tag;
            variableList = out.PARA.variables;
            
            if t>=out.OUTPUT_TIME

                out.timestamp = [out.timestamp t];

                for k=1:size(variableList,1)
                    if any(strcmp(fieldnames(TOP.NEXT.STATVAR), variableList{k,1}))
                        out.result.(variableList{k,1}) = [out.result.(variableList{k,1}); TOP.NEXT.STATVAR.(variableList{k,1})];
                    else
                        out.result.(variableList{k,1}) = [out.result.(variableList{k,1}); NaN];
                    end
                end

                % Set the next OUTPUT_TIME
                out.OUTPUT_TIME = min(out.SAVE_TIME, out.OUTPUT_TIME + out.PARA.output_timestep);
                
                if t>=out.SAVE_TIME
                    % It is time to save all the collected model output to disk
                     
                    if ~(exist([result_path run_name])==7)
                        mkdir([result_path run_name])
                    end
                    CG_out_scalars = out.result;
                    CG_out_scalars.timestamp = out.timestamp;
                    if isempty(out_tag) || all(isnan(out_tag))
                        save([result_path run_name '/' run_name '_scalars_' datestr(t,'yyyymmdd') '.mat'], 'CG_out_scalars')
                    else
                        save([result_path run_name '/' run_name '_scalars_' out_tag '_' datestr(t,'yyyymmdd') '.mat'], 'CG_out_scalars')
                    end
                    
                    % Clear the out structure
                    %make struct "result" and initialize all variables defined by
                    %the user as empty arrays
                    out.timestamp = [];
                    for i =  1:size(out.PARA.variables,1)
                        out.result.(out.PARA.variables{i,1}) = [];
                    end

                    if ~isnan(out.PARA.save_interval)
                        % If save_interval is defined, uptate SAVE_TIME for next save opertion 
                        out.SAVE_TIME = min(forcing.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(out.SAVE_TIME,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
                        % If save_interval is not defined, we will save at the very end of the model run
                        % and thus do not need to update SAVE_TIME (update would fail because save_interval is nan)
					end
                end
            end
        end
        


        %-------------param file generation-----
        function out = param_file_info(out)
            out = provide_PARA(out);

            out.PARA.STATVAR = [];
            out.PARA.options = [];
            out.PARA.class_category = 'OUT';
            
            out.PARA.comment.variables = {'select output variables: T, water, ice, waterIce, XwaterIce, Xwater, Xice, saltConc are supported'};
            out.PARA.options.variables.name = 'H_LIST';
            out.PARA.options.variables.entries_x = {'T' 'water' 'ice' 'waterIce'};
      
            out.PARA.comment.upper_elevation = {'upper elevation of output domain, in m a.s.l.; must match altitude in TILE!'};
            
            out.PARA.comment.lower_elevation = {'lower elevation of output domain, in m a.s.l.; must match altitude in TILE!'};
            
            out.PARA.default_value.target_grid_size = {0.02};
            out.PARA.comment.target_grid_size = {'cell size of regular grid that values are interpolated to'};
           
            out.PARA.default_value.output_timestep = {0.25};
            out.PARA.comment.output_timestep = {'timestep of output [days]'};

            out.PARA.default_value.save_date = {'01.09.'};
            out.PARA.comment.save_date = {'date (dd.mm.) when output file is written'};
            
            out.PARA.default_value.save_interval = {1};
            out.PARA.comment.save_interval = {'interval of output files [years]'};
            
            out.PARA.default_value.tag = {''};
            out.PARA.comment.tag = {'additional tag added to file name'};
        end
        

    end
end