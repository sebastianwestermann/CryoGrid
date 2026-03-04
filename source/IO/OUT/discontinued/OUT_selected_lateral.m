%========================================================================
% CryoGrid OUT class OUT_selected_lateral
% CryoGrid OUT class defining storage format of the output 
% OUT_all_lateral stores identical copies of all GROUND classses (including STATVAR, TEMP, PARA) in the
% stratigraphy for each output timestep and copies of all LATERAL classes.
% Other than that, it is identical to OUT_all.
% The user can specify the save date and the save interval (e.g. yearly
% files), as well as the output timestep (e.g. 6 hourly). The output files
% are in Matlab (".mat") format.
% S. Westermann, T. Ingeman-Nielsen, J. Scheer, October 2020
%========================================================================


classdef OUT_selected_lateral < matlab.mixin.Copyable

    properties
		out_index
        STRATIGRAPHY
        LATERAL
        TIMESTAMP
        MISC
        TEMP
        PARA
        OUTPUT_TIME
        SAVE_TIME
		CONST
	end
    
    
    methods
		


%         function out = initialize_excel(out)
%             
%         end
        
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
			
			out.OUTPUT_TIME = forcing.PARA.start_time + out.PARA.output_timestep;
            if isempty(out.PARA.save_interval) || isnan(out.PARA.save_interval) 
                out.SAVE_TIME = forcing.PARA.end_time;
            else
                out.SAVE_TIME = min(forcing.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(forcing.PARA.start_time,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
            end
            out.TEMP = struct();
        end
        
        %-------time integration----------------
		
		%function out = store_OUT(out, t, TOP, BOTTOM, forcing, run_number, timestep, result_path)
            
        function out = store_OUT(out, tile)
            
             t = tile.t;
             TOP = tile.TOP; 
             BOTTOM = tile.BOTTOM;
             forcing = tile.FORCING;
             run_name = tile.PARA.run_name;
             result_path = tile.PARA.result_path;
             timestep = tile.timestep;
             out_tag = out.PARA.tag;
             

            if t>=out.OUTPUT_TIME
        				% It is time to collect output
                % Store the current state of the model in the out structure.
                
                disp([datestr(t)])

                out.TIMESTAMP=[out.TIMESTAMP t];
                
                CURRENT =TOP.NEXT;
                if isprop(CURRENT, 'CHILD') && CURRENT.CHILD ~= 0
                    out.MISC=[out.MISC [CURRENT.CHILD.STATVAR.T(1,1); CURRENT.CHILD.STATVAR.layerThick(1,1)]]; 
                else
                    out.MISC=[out.MISC [NaN; NaN]];
                end
                result={};
                while ~isequal(CURRENT, BOTTOM)
                    if isprop(CURRENT, 'CHILD') && CURRENT.CHILD ~= 0
                        for ii=1:size(out.PARA.variables,1)
                            res.class_name = class(CURRENT.CHILD);
                            res.is_CHILD = 1;
                            if isfield(CURRENT.CHILD.STATVAR, out.PARA.variables{ii,1})
                                res.(out.PARA.variables{ii,1}) = CURRENT.CHILD.STATVAR.(out.PARA.variables{ii,1});
                            end
                        end
                        result=[result; {res}];
                    end
                    for ii=1:size(out.PARA.variables,1)
                            res.class_name = class(CURRENT);
                            res.is_CHILD = 0;
                            if isfield(CURRENT.STATVAR, out.PARA.variables{ii,1})
                                res.(out.PARA.variables{ii,1}) = CURRENT.STATVAR.(out.PARA.variables{ii,1});
                            end

                        end
                    result=[result; {res}];
                    CURRENT = CURRENT.NEXT;
                end
                out.STRATIGRAPHY{1,size(out.STRATIGRAPHY,2)+1} = result;
                
                %lateral, read only STATVAR and PARA---
                result={};
                ia_classes = tile.LATERAL.IA_CLASSES;
                for i=1:size(ia_classes,1)
                    res = copy(ia_classes{i,1});
                    vars = fieldnames(res);
                    for j=1:size(vars,1)
                        if ~strcmp(vars{j,1}, 'PARA') && ~strcmp(vars{j,1}, 'STATVAR')
                           res.(vars{j,1}) = []; 
                        end
                    end
                    result=[result; {res}];                    
                end
                out.LATERAL{1,size(out.LATERAL, 2)+1} = result;
                %---
                
				% Set the next OUTPUT_TIME
                out.OUTPUT_TIME = out.OUTPUT_TIME + out.PARA.output_timestep;
				
                if t>=out.SAVE_TIME
					% It is time to save all the collected model output to disk
					
				    if ~(exist([result_path run_name])==7)
				    	mkdir([result_path run_name])
                    end
                    if isempty(out_tag) || all(isnan(out_tag))
                        save([result_path run_name '/' run_name '_' datestr(t,'yyyymmdd') '.mat'], 'out')
                    else
                        save([result_path run_name '/' run_name '_' out_tag '_' datestr(t,'yyyymmdd') '.mat'], 'out')
                    end
				    				    
					% Clear the out structure
					out.STRATIGRAPHY=[];
				    out.LATERAL=[];
				    out.TIMESTAMP=[];
				    out.MISC=[];
					if ~isnan(out.PARA.save_interval)
                        % If save_interval is defined, uptate SAVE_TIME for next save opertion 
                        out.SAVE_TIME = min(forcing.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(out.SAVE_TIME,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
                        % If save_interval is not defined, we will save at the very end of the model run
                        % and thus do not need to update SAVE_TIME (update would fail because save_interval is nan)
					end

                end
            end
        end
        
        function out = reset_timestamp_out(out,tile)
            if size(out.TIMESTAMP,2)>0
                delete_timestamps = find(out.TIMESTAMP(1,:)>tile.t);
                if ~isempty(delete_timestamps)
                    
                    out.STRATIGRAPHY(:,delete_timestamps) = [];
                    out.LATERAL(:,delete_timestamps) = [];
                    out.TIMESTAMP(:,delete_timestamps) = [];
                end
            end
            out.OUTPUT_TIME = tile.t + out.PARA.output_timestep;
        end
        
        %-------------param file generation-----
        function out = param_file_info(out)
            out = provide_PARA(out);

            out.PARA.STATVAR = [];
            out.PARA.options = [];
            out.PARA.class_category = 'OUT';
           
            out.PARA.default_value.output_timestep = {0.25};
            out.PARA.comment.output_timestep = {'timestep of output [days]'};

            out.PARA.default_value.save_date = {'01.09.'};
            out.PARA.comment.save_date = {'date (dd.mm.) when output file is written'};
            
            out.PARA.default_value.save_interval = {1};
            out.PARA.comment.save_interval = {'interval of output files [years]'};
            
            out.PARA.default_value.tag = {''};
            out.PARA.comment.tag = {'additional tag added to file name'};
        end

%         function xls_out = write_excel(out)
% 			% XLS_OUT  Is a cell array corresponding to the class-specific content of the parameter excel file (refer to function write_controlsheet).
% 			
%             xls_out = {'OUT','index',NaN,NaN;'OUT_all',1,NaN,NaN;'output_timestep',0.250000000000000,'[days]',NaN;'save_date','01.09.','provide in format dd.mm.',NaN;'save_interval',1,'[y]','if left empty, the entire output will be written out at the end';'OUT_END',NaN,NaN,NaN};
%         end
        
    end
end