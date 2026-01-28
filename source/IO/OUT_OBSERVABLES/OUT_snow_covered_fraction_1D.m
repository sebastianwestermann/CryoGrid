classdef OUT_snow_covered_fraction_1D < matlab.mixin.Copyable
    
    properties
        PARA
        CONST
        TEMP
        STATVAR
        TIMESTAMP
        OUTPUT_TIME
        SAVE_TIME
    end
    
    methods

        function out = provide_PARA(out)
            out.PARA.output_timestep = [];
            out.PARA.save_date = [];
            out.PARA.save_interval = [];
            out.PARA.timestamps = [];
            out.PARA.reference_SWE = []; 
            out.PARA.add_Xice = [];
            out.PARA.tag = [];

            out.PARA.timestamps = [];
        end
        
        function out = provide_CONST(out)
            
        end
        
        function out = provide_STATVAR(out)
            
        end
        
        function out = finalize_init(out, tile)
            out.STATVAR.fsca = []; 
            out.TIMESTAMP = [];

            if ~isempty(out.PARA.timestamps) %DA mode, defined timestamps written by DA
                out.OUTPUT_TIME = out.PARA.timestamps(find(out.PARA.timestamps(:,1)-tile.PARA.start_time > 0, 1, 'first'), 1);
                out.SAVE_TIME = Inf;
            else
                out.OUTPUT_TIME = tile.PARA.start_time + out.PARA.output_timestep;
                if isempty(out.PARA.save_interval) || isnan(out.PARA.save_interval)
                    out.SAVE_TIME = tile.FORCING.PARA.end_time;
                else
                    out.SAVE_TIME = min(tile.FORCING.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(tile.PARA.start_time,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
                end
            end
            if isempty(out.OUTPUT_TIME)
                out.OUTPUT_TIME = Inf;
            end
            
        end
        
        function out = store_OUT(out, tile)

            if tile.t >= out.OUTPUT_TIME

                out.TIMESTAMP = [out.TIMESTAMP tile.t];

                if ~isempty(out.PARA.timestamps)
                    out.OUTPUT_TIME = out.PARA.timestamps(find(out.PARA.timestamps(:,1)-tile.t>0, 1, 'first'), 1);
                else
                    out.OUTPUT_TIME = min(out.SAVE_TIME, out.OUTPUT_TIME + out.PARA.output_timestep);
                end
                if isempty(out.OUTPUT_TIME)
                    out.OUTPUT_TIME = Inf;
                end

                result = 0;
                CURRENT = tile.TOP.NEXT;
                while ~(strcmp(class(CURRENT), 'Bottom'))
                    class_name = class(CURRENT);
                    if strcmp(class_name(1:4), 'SNOW')
                        result = 1;
                        break
                    end
                    if out.PARA.add_Xice == 1
                        if any(strcmp(properties(CURRENT), 'CHILD'))
                            if isempty(CURRENT.CHILD) || CURRENT.CHILD == 0
                                ice_height = CURRENT.STATVAR.Xice(1,1) ./ CURRENT.STATVAR.area(1,1);
                                result = min(1,max(0, ice_height./ out.PARA.reference_SWE));
                            else
                                ice_height = (CURRENT.CHILD.STATVAR.ice + CURRENT.STATVAR.Xice(1,1)) ./ CURRENT.STATVAR.area(1,1);
                                result = min(1,max(0, ice_height./ out.PARA.reference_SWE));
                            end
                            break
                        end
                    else
                        if any(strcmp(properties(CURRENT), 'CHILD'))
                            ice_height = CURRENT.CHILD.STATVAR.ice ./ CURRENT.STATVAR.area(1,1);
                            result = min(1,max(0, ice_height./ out.PARA.reference_SWE));
                        end
                        break
                    end
                    CURRENT = CURRENT.NEXT;
                end
                
                out.STATVAR.fsca = [out.STATVAR.fsca; result];
            end
            
            if tile.t>=out.SAVE_TIME

                if ~(exist([tile.PARA.result_path tile.PARA.run_name])==7)
                    mkdir([tile.PARA.result_path tile.PARA.run_name])
                end
                CG_out.fsca = out.STATVAR.fsca;
                CG_out.timestamp = out.TIMESTAMP';
                if isempty(out.PARA.tag) || all(isnan(out.PARA.tag))
                    save([tile.PARA.result_path tile.PARA.run_name '/' tile.PARA.run_name '_FSCA_' datestr(tile.t,'yyyymmdd') '.mat'], 'CG_out')
                else
                    save([tile.PARA.result_path tile.PARA.run_name '/' tile.PARA.run_name '_FSCA_' out.PARA.tag '_' datestr(tile.t,'yyyymmdd') '.mat'], 'CG_out')
                end

                % Clear the out structure
                out.STATVAR.fsca = [];
                out.TIMESTAMP=[];
                if ~isnan(out.PARA.save_interval)
                    out.SAVE_TIME = min(tile.FORCING.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(out.SAVE_TIME,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
                end
            end

        end

        function result = move_out2obs(out)
            result = out.STATVAR.fsca;
        end
            
        
        function out = reset_new_stratigraphy(out, tile)
            
        end
        
    end
end

