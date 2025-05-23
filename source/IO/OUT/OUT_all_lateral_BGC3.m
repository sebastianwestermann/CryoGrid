%========================================================================
% CryoGrid OUT class defining storage format of the output
% OUT_all_lateral stores identical copies of all GROUND classses (including STATVAR, TEMP, PARA) in the
% stratigraphy for each output timestep and copies of all LATERAL classes.
% Other than that, it is identical to OUT_all.
% The user can specify the save date and the save interval (e.g. yearly
% files), as well as the output timestep (e.g. 6 hourly). The output files
% are in Matlab (".mat") format.
% S. Westermann, T. Ingeman-Nielsen, J. Scheer, October 2020
%========================================================================


classdef OUT_all_lateral_BGC3 < matlab.mixin.Copyable
    
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
            % INITIALIZE_PARA  Initializes PARA structure.
            
            out.PARA.output_timestep = [];
            out.PARA.save_date = [];
            out.PARA.save_interval = [];
            out.PARA.number_of_save_intervals = [];
            out.PARA.number_of_gap_intervals = [];
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
                save_day = str2num(out.PARA.save_date(1:2));
                save_month = str2num(out.PARA.save_date(4:5));
                [year, ~,~] = datevec(forcing.PARA.start_time);
                %out.SAVE_TIME = min(forcing.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(forcing.PARA.start_time,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
                out.SAVE_TIME = min(forcing.PARA.end_time,  datenum(year + out.PARA.save_interval, save_month, save_day));

            end
            out.TEMP = struct();
            
            out.TEMP.save_interval_count = 0;
            out.TEMP.save_active = 1;
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
            
            
            if t>=out.OUTPUT_TIME %|| (tile.timestep <1e-12 && t>datenum(2014,4,1))
                
                [a,b,c] = datevec(t);
                disp(['year: ', num2str(a) ', month: ' num2str(b) ', day: ' num2str(c)])
                
                if out.TEMP.save_active
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
                            res=copy(CURRENT.CHILD);
                            res.NEXT =[]; res.PREVIOUS=[]; res.IA_NEXT=[]; res.IA_PREVIOUS=[];  res.PARENT = []; %cut all dependencies
                            result=[result; {res}];
                        end
                        res = copy(CURRENT);
                        if isprop(res, 'BGC')
                            if t>=out.SAVE_TIME
                                res.STATVAR.BGC = res.BGC.STATVAR;
                            end
                            res.BGC = [];
                            res.IA_BGC = [];
                        end
                        if isprop(res, 'LUT')
                            res.LUT =[];  %remove look-up tables, runs out of memory otherwise
                        end
                        if isprop(res, 'STORE')
                            res.STORE = [];
                        end
                        if isprop(res, 'MODEL_CLASS')
                            res.MODEL_CLASS = [];
                        end
                        if isprop(res, 'READ_OUT')
                            res.READ_OUT =[];  %remove look-up tables, runs out of memory otherwise
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
                    
%                     out.OUTPUT_TIME = out.OUTPUT_TIME + out.PARA.output_timestep;
                    if t>=out.SAVE_TIME  %|| (tile.timestep <1e-12 && t>datenum(2014,4,1))
                        if ~(exist([result_path run_name])==7)
                            mkdir([result_path run_name])
                        end

                        [year, ~,~] = datevec(t);
                        save([result_path run_name '/' run_name '_' num2str(year) '.mat'], 'out')
                        out.STRATIGRAPHY=[];
                        out.LATERAL=[];
                        out.TIMESTAMP=[];
                        out.MISC=[];
                        %out.SAVE_TIME = min(forcing.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(out.SAVE_TIME,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
                    end
                end
                out.OUTPUT_TIME = out.OUTPUT_TIME + out.PARA.output_timestep;
                if t>=out.SAVE_TIME
                    save_day = str2num(out.PARA.save_date(1:2));
                    save_month = str2num(out.PARA.save_date(4:5));
                    [year, ~,~] = datevec(t);
                    out.SAVE_TIME = min(forcing.PARA.end_time, datenum(year + out.PARA.save_interval, save_month, save_day));
                    %out.SAVE_TIME = min(forcing.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(out.SAVE_TIME,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
                    out.TEMP.save_interval_count = out.TEMP.save_interval_count + 1;
                    
                    if out.TEMP.save_interval_count == out.PARA.number_of_save_intervals && out.TEMP.save_active
                        out.TEMP.save_active = 0;
                        out.TEMP.save_interval_count = 0;
                    elseif out.TEMP.save_interval_count == out.PARA.number_of_gap_intervals && ~out.TEMP.save_active
                        out.TEMP.save_active = 1;
                        out.TEMP.save_interval_count = 0;
                    end
                end
                
            end
        end

        function  [uppermost_class, lowermost_class] = restore_stratigraphy_from_OUT(out, restart_time, tile)
            if ~isempty(restart_time)
                pos = find(out.TIMESTAMP(1,:) == restart_time);
            else
                pos = size(out.TIMESTAMP,2);
            end

            uppermost_class_name = class(out.STRATIGRAPHY{1,pos}{1,1});
            if ~strcmp(uppermost_class_name(1:4), 'SNOW')
                uppermost_class = reset_from_OUT(out.STRATIGRAPHY{1,pos}{1,1}, tile);
                next_cell = 2;
            else
                uppermost_class = reset_from_OUT(out.STRATIGRAPHY{1,pos}{2,1}, tile);
                next_cell = 3;
            end
            CURRENT = uppermost_class;
            for i=next_cell:size(out.STRATIGRAPHY{1,pos},1)
                CURRENT.NEXT = reset_from_OUT(out.STRATIGRAPHY{1,pos}{i,1}, tile);
                CURRENT.NEXT.PREVIOUS = CURRENT;
                CURRENT = CURRENT.NEXT;
            end
            lowermost_class = CURRENT;

            %restore BGC 
            CURRENT = uppermost_class;
            for i=next_cell:size(out.STRATIGRAPHY{1,pos},1)
                variables = fieldnames(CURRENT);
                if any(strcmp(variables, 'BGC'))
                    a=str2func(CURRENT.PARA.BGC_CLASS);
                    CURRENT.BGC = a();
                    CURRENT.BGC = provide_PARA(CURRENT.BGC);
                    CURRENT.BGC = provide_CONST(CURRENT.BGC);
                    CURRENT.BGC = finalize_init(CURRENT.BGC, tile);
                    CURRENT.BGC.STATVAR = CURRENT.STATVAR.BGC;

                    CURRENT.STATVAR.BGC = [];
                    CURRENT.IA_BGC = IA_BGC_RichardsEqW_Xice();
                    CURRENT.IA_BGC.BGC = CURRENT.BGC;
                    CURRENT.IA_BGC.GROUND = CURRENT;
                    CURRENT.BGC.IA_BGC = CURRENT.IA_BGC;
                    finalize_init(CURRENT.IA_BGC, tile);
                end
                CURRENT = CURRENT.NEXT;
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
            
%             out.PARA.default_value.tag = {''};
%             out.PARA.comment.tag = {'additional tag added to file name'};
        end
        
    end
end