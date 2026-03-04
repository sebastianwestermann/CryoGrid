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


classdef OUT_all < OUT_BASE
 

    properties
		STRATIGRAPHY
        TIMESTAMP
    end
    
    
    methods
        
        %initialization
        
        function out = provide_PARA(out)         

            out.PARA.output_timestep = [];
            out.PARA.save_date = [];
            out.PARA.save_interval = [];
            out.PARA.tag = [];
            out.PARA.tag2 = [];

        end
            
        function out = store_OUT(out, tile)

            if tile.t >= out.OUTPUT_TIME
                
                disp([datestr(tile.t)])

                out = state2out(out, tile);
                out.OUTPUT_TIME = out.OUTPUT_TIME + out.PARA.output_timestep;

                if tile.t >= out.SAVE_TIME
                    out = out2file_all(out, tile);
                end
                
            end
        end

        function out = state2out(out, tile)

            out.TIMESTAMP=[out.TIMESTAMP tile.t];
                
                CURRENT = tile.TOP.NEXT;
                result={};
                while ~isequal(CURRENT, tile.BOTTOM)
                    if isprop(CURRENT, 'CHILD') && CURRENT.CHILD ~= 0
                        res=copy(CURRENT.CHILD);
                        res.NEXT =[]; res.PREVIOUS=[]; res.IA_NEXT=[]; res.IA_PREVIOUS=[];  res.PARENT = []; %cut all dependencies
                        result=[result; {res}];
                    end
                    res = copy(CURRENT);
                    if isprop(res, 'LUT')
                        res.LUT =[];  %remove look-up tables, runs out of memory otherwise
                    end
                    if isprop(res, 'READ_OUT')
                        res.READ_OUT =[];  %remove look-up tables, runs out of memory otherwise
                    end
                    if isprop(res, 'STORE')
                        res.STORE = [];
                    end
                    res.NEXT =[]; res.PREVIOUS=[]; res.IA_NEXT=[]; res.IA_PREVIOUS=[];  %cut all dependencies
                    if isprop(res, 'CHILD')
                        res.CHILD = [];
                        res.IA_CHILD =[];
                    end
                    result=[result; {res}];
                    CURRENT = CURRENT.NEXT;
                end
                out.STRATIGRAPHY{1,size(out.STRATIGRAPHY,2)+1} = result;

        end
        
        function out = reset_timestamp_out(out,tile)
            if size(out.TIMESTAMP,2)>0
                delete_timestamps = find(out.TIMESTAMP(1,:)>tile.t);
                if ~isempty(delete_timestamps)
                    
                    out.STRATIGRAPHY(:,delete_timestamps) = [];
                    out.TIMESTAMP(:,delete_timestamps) = [];
                    out.MISC(:,delete_timestamps) = [];
                end
            end
            out.OUTPUT_TIME = tile.t + out.PARA.output_timestep;
            if ~isnan(out.PARA.save_interval)
                potential_save_time =  datenum([out.PARA.save_date num2str(str2num(datestr(out.SAVE_TIME,'yyyy')) - out.PARA.save_interval)], 'dd.mm.yyyy');
                if tile.t <= potential_save_time
                    out.SAVE_TIME = potential_save_time;
                end
            end
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
%             % XLS_OUT  Is a cell array corresponding to the class-specific content of the parameter excel file (refer to function write_controlsheet).
%             
%             xls_out = {'OUT','index',NaN,NaN;'OUT_all',1,NaN,NaN;'output_timestep',0.250000000000000,'[days]',NaN;'save_date','01.09.','provide in format dd.mm.',NaN;'save_interval',1,'[y]','if left empty, the entire output will be written out at the end';'OUT_END',NaN,NaN,NaN};
%         end
        
    end
end