

classdef OUT_MULTITILE_snow_simple < matlab.mixin.Copyable

    properties
        
        PARA
        STATVAR
        CONST
        TEMP
        OUTPUT_TIME
        SAVE_TIME
    end
    
    methods

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

            out = reset_STATVAR(out); %initializes STATVAR
            
            if isempty(out.PARA.save_interval) || isnan(out.PARA.save_interval) 
                out.SAVE_TIME = tile.FORCING.PARA.end_time;
            else
                if datenum([out.PARA.save_date datestr(tile.FORCING.PARA.start_time,'yyyy')], 'dd.mm.yyyy') <= tile.FORCING.PARA.start_time + 30
                    out.SAVE_TIME = min(tile.FORCING.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(tile.FORCING.PARA.start_time,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
                else
                    out.SAVE_TIME = min(tile.FORCING.PARA.end_time,  datenum([out.PARA.save_date datestr(tile.FORCING.PARA.start_time,'yyyy')], 'dd.mm.yyyy'));
                end
            end
            
            if isempty(out.PARA.output_timestep) || isnan(out.PARA.output_timestep)
                out.OUTPUT_TIME = out.SAVE_TIME;
            else
                out.OUTPUT_TIME = tile.FORCING.PARA.start_time + out.PARA.output_timestep;
            end
            
        end
        
        function out = store_OUT(out, tile)

            ground = tile.SUBSURFACE_CLASS;
            if tile.t >= out.OUTPUT_TIME

                disp(datestr(tile.t))
                out.TEMP.snow_water = [out.TEMP.snow_water; sum(ground.STATVAR.SWE_water,1)];
                out.TEMP.snow_SWE = [out.TEMP.snow_SWE; sum(ground.STATVAR.SWE,1)];
                out.TEMP.snow_GST = [out.TEMP.snow_GST; ground.STATVAR.T];
                out.TEMP.timestamp = [out.TEMP.timestamp tile.t];

                if tile.t>=out.SAVE_TIME

                    run_name = tile.PARA.run_name;
                    result_path = tile.PARA.result_path;
                    out_tag = out.PARA.tag;
                    t = tile.t;

                    CG_out = out.TEMP;

                    if ~(exist([result_path run_name])==7)
                        mkdir([result_path run_name])
                    end
                    if isempty(out_tag) || all(isnan(out_tag))
                        save([result_path run_name '/' run_name '_' datestr(t,'yyyymmdd') '.mat'], 'CG_out')
                    else
                        save([result_path run_name '/' run_name '_' out_tag '_' datestr(t,'yyyymmdd') '.mat'], 'CG_out')
                    end

                    out = reset_STATVAR(out);
                    out.SAVE_TIME = min(tile.FORCING.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(out.SAVE_TIME,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
                end

                if isempty(out.PARA.output_timestep) || isnan(out.PARA.output_timestep)
                    out.OUTPUT_TIME = out.SAVE_TIME;
                else
                    out.OUTPUT_TIME = min(out.SAVE_TIME, out.OUTPUT_TIME + out.PARA.output_timestep);
                end
            end
        end


            
        function out = reset_STATVAR(out)
            out.TEMP.snow_water = [];
            out.TEMP.snow_SWE = [];
            out.TEMP.snow_GST = [];
            out.TEMP.timestamp = [];
        end

    end
end

