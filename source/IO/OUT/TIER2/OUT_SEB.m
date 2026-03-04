%========================================================================
% CryoGrid OUT class OUT_SEB for storing surface energy balance variables
%
% S. Westermann, Oct 2022
%========================================================================


classdef OUT_SEB < OUT_BASE
 

    properties

        STATVAR
        
    end
    
    
    methods
                
        function out = provide_PARA(out)         

            out.PARA.variables = [];
            out.PARA.output_timestep = [];
            out.PARA.save_date = [];
            out.PARA.save_interval = [];
            out.PARA.tag = [];
            out.PARA.tag2 = [];
        end

        
        
        function out = finalize_init(out, tile)

            out = finalize_init@OUT_BASE(out, tile);

            for i = 1:size(out.PARA.variables,1)
               out.STATVAR.(out.PARA.variables{i,1}) = []; 
            end
            out.STATVAR.timestamp = [];

            out.TEMP.output_var = 'CG_SEB';
            out.TEMP.keyword = 'SEB';
      
        end
        
        function out = store_OUT(out, tile)

            if tile.t >= out.OUTPUT_TIME

                disp([datestr(tile.t)])

                out = state2out(out, tile);
                out.OUTPUT_TIME = out.OUTPUT_TIME + out.PARA.output_timestep;

                if tile.t >= out.SAVE_TIME
                    out = out2file_CG(out, tile);
                end

            end
        end
            
        function out = state2out(out, tile)

            out.STATVAR.timestamp = [out.STATVAR.timestamp; tile.t];
            CURRENT = tile.TOP.NEXT;
            for i=1:size(out.PARA.variables,1)
                if strcmp(out.PARA.variables{i,1}, 'Sin') || strcmp(out.PARA.variables{i,1}, 'Lin')
                    out.STATVAR.(out.PARA.variables{i,1}) = [out.STATVAR.(out.PARA.variables{i,1}); tile.FORCING.TEMP.(out.PARA.variables{i,1})];
                else
                    out.STATVAR.(out.PARA.variables{i,1}) = [out.STATVAR.(out.PARA.variables{i,1}); CURRENT.STATVAR.(out.PARA.variables{i,1})];
                end
            end

        end

% 
%         function out = store_OUT(out, tile)           
% 
%             t = tile.t;
%             TOP = tile.TOP; 
%             BOTTOM = tile.BOTTOM;
%             forcing = tile.FORCING;
%             run_name = tile.PARA.run_name; %tile.RUN_NUMBER;
%             result_path = tile.PARA.result_path;            
%             timestep = tile.timestep;
%             out_tag = out.PARA.tag;
% 
%             if t>=out.OUTPUT_TIME
%                 % It is time to collect output
%                 % Store the current state of the model in the out structure.
% 
% %                 disp([datestr(t)])                
%                 out.timestamp = [out.timestamp t];
% 
%                 CURRENT =TOP.NEXT;
%                 for i=1:size(out.PARA.variables,1)
%                     if strcmp(out.PARA.variables{i,1}, 'Sin') || strcmp(out.PARA.variables{i,1}, 'Lin')
%                        out.result.(out.PARA.variables{i,1}) = [out.result.(out.PARA.variables{i,1}) tile.FORCING.TEMP.(out.PARA.variables{i,1})];
%                     else
%                         out.result.(out.PARA.variables{i,1}) = [out.result.(out.PARA.variables{i,1}) CURRENT.STATVAR.(out.PARA.variables{i,1})];
%                     end
%                 end
% 
%                 % Set the next OUTPUT_TIME
%                 out.OUTPUT_TIME = min(out.SAVE_TIME, out.OUTPUT_TIME + out.PARA.output_timestep);
% 
%                 if t>=out.SAVE_TIME
%                     % It is time to save all the collected model output to disk
% 
%                     if ~(exist([result_path run_name])==7)
%                         mkdir([result_path run_name])
%                     end
%                     OUT_SEB = out.result;
%                     OUT_SEB.timestamp = out.timestamp;
%                     if isempty(out_tag) || all(isnan(out_tag))
%                         save([result_path run_name '/' run_name '_SEB_' datestr(t,'yyyymmdd') '.mat'], 'OUT_SEB')
%                     else
%                         save([result_path run_name '/' run_name '_SEB_' out_tag '_' datestr(t,'yyyymmdd') '.mat'], 'OUT_SEB')
%                     end
% 
%                     % Clear the out structure
%                     out.MISC=[];
%                     %make struct "result" and initialize all variables defined by
%                     %the user as empty arrays
%                     out.timestamp = [];
%                     for i =  1:size(out.PARA.variables,1)
%                         out.result.(out.PARA.variables{i,1}) = [];
%                     end
%                     out.result.depths = [];
%                     out.result.class_number = [];
% 
%                     if ~isnan(out.PARA.save_interval)
%                         % If save_interval is defined, uptate SAVE_TIME for next save opertion 
%                         out.SAVE_TIME = min(forcing.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(out.SAVE_TIME,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
%                         % If save_interval is not defined, we will save at the very end of the model run
%                         % and thus do not need to update SAVE_TIME (update would fail because save_interval is nan)
% 					end
%                 end
%             end
%         end
        
                

        %-------------param file generation-----
        function out = param_file_info(out)
            out = provide_PARA(out);

            out.PARA.STATVAR = [];
            out.PARA.options = [];
            out.PARA.class_category = 'OUT';
            
            out.PARA.comment.variables = {'select output variables: Sin, Sout, Lin, Lout, Qe, Qh are supported'};
            out.PARA.options.variables.name = 'H_LIST';
            out.PARA.options.variables.entries_x = {'Sin' 'Sout' 'Lin' 'Lout' 'Qh' 'Qe'};
           
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