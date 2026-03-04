%========================================================================
% CryoGrid OUT class OUT_all_lateral
% CryoGrid OUT class defining storage format of the output 
% OUT_all_lateral stores identical copies of all GROUND classses (including STATVAR, TEMP, PARA) in the
% stratigraphy for each output timestep and copies of all LATERAL classes.
% Other than that, it is identical to OUT_all.
% The user can specify the save date and the save interval (e.g. yearly
% files), as well as the output timestep (e.g. 6 hourly). The output files
% are in Matlab (".mat") format.
% S. Westermann, T. Ingeman-Nielsen, J. Scheer, October 2020
%========================================================================


classdef OUT_MB < OUT_BASE

    properties
        STATVAR
	end
    
    
    methods
		
 
        function out = provide_PARA(out)         

            out.PARA.output_timestep = [];
            out.PARA.save_date = [];
            out.PARA.save_interval = [];
            out.PARA.tag = [];
            out.PARA.tag2 = [];
        end
                
        function out = finalize_init(out, tile)

            out = finalize_init@OUT_BASE(out, tile);

            out.STATVAR.timestamp = [];
            out.STATVAR.SMB = [];
            out.STATVAR.runoff = [];

            out.TEMP.keyword = 'MB';
            out.TEMP.output_var = 'CG_MB';
        end

        function out = store_OUT(out, tile)

            if tile.t >= out.OUTPUT_TIME
                
                disp([datestr(tile.t)])

                out = state2out(out, tile);
                out.OUTPUT_TIME = out.OUTPUT_TIME + out.PARA.output_timestep;

                if tile.t >= out.SAVE_TIME
                    out = out2file_CG(out, tile);
                    out.STATVAR.timestamp = [];
                    out.STATVAR.SMB = [];
                    out.STATVAR.runoff = [];
                end
                
            end
        end

        function out = state2out(out, tile)

            disp([datestr(tile.t)])
            out.STATVAR.timestamp=[out.STATVAR.timestamp; tile.t];

            CURRENT = tile.TOP.NEXT;
            [smb,runoff]=deal([]);
            while ~isequal(CURRENT, tile.BOTTOM)
                if isprop(CURRENT, 'CHILD') && CURRENT.CHILD ~= 0
                    smb=[smb; CURRENT.CHILD.STATVAR.smb];
                    runoff=[runoff; CURRENT.CHILD.STATVAR.runoff];
                end
                smb=[smb; CURRENT.STATVAR.smb];
                runoff=[runoff; CURRENT.STATVAR.runoff];

                CURRENT = CURRENT.NEXT;
            end

            out.STATVAR.SMB = [out.STATVAR.SMB; nansum(smb)];
            out.STATVAR.runoff = [out.STATVAR.runoff; nansum(runoff)];
        end

            
%         function out = store_OUT(out, tile)
% 
%              t = tile.t;
%              TOP = tile.TOP; 
%              BOTTOM = tile.BOTTOM;
%              forcing = tile.FORCING;
%              run_name = tile.PARA.run_name;
%              result_path = tile.PARA.result_path;
%              timestep = tile.timestep;
%              out_tag = out.PARA.tag;
% 
%             if t>=out.OUTPUT_TIME
%         				% It is time to collect output
%                 % Store the current state of the model in the out structure.
% 
%                 disp([datestr(t)])
% 
%                 out.TIMESTAMP=[out.TIMESTAMP t];
% 
%                 CURRENT =TOP.NEXT;
% 
%                 [smb,runoff]=deal([]);
%                 while ~isequal(CURRENT, BOTTOM)
%                     if isprop(CURRENT, 'CHILD') && CURRENT.CHILD ~= 0
%                         res=copy(CURRENT.CHILD);
%                         res.NEXT =[]; res.PREVIOUS=[]; res.IA_NEXT=[]; res.IA_PREVIOUS=[];  res.PARENT = []; %cut all dependencies
%                         smb=[smb; res.STATVAR.smb];
%                         runoff=[runoff; res.STATVAR.runoff];
% 
%                     end
%                     res = copy(CURRENT);
% 
%                     smb=[smb; res.STATVAR.smb];
%                     runoff=[runoff; res.STATVAR.runoff];
% 
%                     CURRENT = CURRENT.NEXT;
%                 end
% 
%                 out.SMB = [out.SMB; nansum(smb)];
%                 out.runoff = [out.runoff; nansum(runoff)];
% 
% 				% Set the next OUTPUT_TIME
%                 out.OUTPUT_TIME = out.OUTPUT_TIME + out.PARA.output_timestep;
% 
% 
%                 if t>=out.SAVE_TIME
% 					% It is time to save all the collected model output to disk
% 
% 				    if ~(exist([result_path run_name])==7)
% 				    	mkdir([result_path run_name])
%                     end
%                     if isempty(out_tag) || all(isnan(out_tag))
%                         save([result_path run_name '/' run_name '_' datestr(t,'yyyymmdd') '_MB.mat'], 'out')
%                     else
%                         save([result_path run_name '/' run_name '_' out_tag '_' datestr(t,'yyyymmdd') '_MB.mat'], 'out')
%                     end
% 
% 					% Clear the out structure
% 					out.SMB=[];
% % 				    out.LATERAL=[];
% 				    out.TIMESTAMP=[];
% % 				    out.MISC=[];
% 					if ~isnan(out.PARA.save_interval)
%                         % If save_interval is defined, uptate SAVE_TIME for next save opertion 
%                         out.SAVE_TIME = min(forcing.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(out.SAVE_TIME,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
%                         % If save_interval is not defined, we will save at the very end of the model run
%                         % and thus do not need to update SAVE_TIME (update would fail because save_interval is nan)
% 					end
% 
%                 end
%             end
%         end
        
        %-------------param file generation-----
        function out = param_file_info(out)
%             out = provide_PARA(out);

%             out.PARA.STATVAR = [];
%             out.PARA.options = [];
%             out.PARA.class_category = 'OUT';
%            
%             out.PARA.default_value.output_timestep = {0.25};
%             out.PARA.comment.output_timestep = {'timestep of output [days]'};
% 
%             out.PARA.default_value.save_date = {'01.09.'};
%             out.PARA.comment.save_date = {'date (dd.mm.) when output file is written'};
%             
%             out.PARA.default_value.save_interval = {1};
%             out.PARA.comment.save_interval = {'interval of output files [years]'};
%             
%             out.PARA.default_value.tag = {''};
%             out.PARA.comment.tag = {'additional tag added to file name'};
        end

%         function xls_out = write_excel(out)
% 			% XLS_OUT  Is a cell array corresponding to the class-specific content of the parameter excel file (refer to function write_controlsheet).
% 			
%             xls_out = {'OUT','index',NaN,NaN;'OUT_all',1,NaN,NaN;'output_timestep',0.250000000000000,'[days]',NaN;'save_date','01.09.','provide in format dd.mm.',NaN;'save_interval',1,'[y]','if left empty, the entire output will be written out at the end';'OUT_END',NaN,NaN,NaN};
%         end
        
    end
end