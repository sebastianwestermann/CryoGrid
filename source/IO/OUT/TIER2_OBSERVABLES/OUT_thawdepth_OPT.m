classdef OUT_thawdepth_OPT < OUT_BASE_OBSERVABLES & OUT_thawdepth
    
    properties

    end
    
    methods

        function out = provide_PARA(out)
            out.PARA.output_timestep = [];
            out.PARA.save_date = [];
            out.PARA.save_interval = [];
            out.PARA.tag = [];
            out.PARA.tag2 = [];
            out.PARA.timestamps = [];
            out.PARA.write_file_mode = 0; %empty/default setting: DA mode, do not write any output; % 0: normal behaviour, no DA, get normal output times and run regularly; 1: store after DA; 2: store after each run
      
        end  
        
        function out = finalize_init(out, tile)

            out = finalize_init@OUT_BASE_OBSERVABLES(out, tile);

            out.STATVAR.timestamp = [];
            out.STATVAR.thawdepth = [];

            out.TEMP.keyword = 'thawdepth';
            out.TEMP.output_var = 'CG_thawdepth';

        end


        function out = store_OUT(out, tile)

            if tile.t >= out.OUTPUT_TIME

                disp([datestr(tile.t)])

                out = state2out(out, tile);
                
                out = set_output_time(out, tile);
                
                if tile.t >= out.SAVE_TIME
                    if out.PARA.write_file_mode == 0
                        out = out2file_CG(out, tile);
                    elseif out.PARA.write_file_mode == 1
                        out.SAVE_TIME = min(tile.FORCING.PARA.end_time, datenum([out.PARA.save_date num2str(str2num(datestr(out.SAVE_TIME,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
                        out.TEMP.write_out_init = 1;
                    elseif out.PARA.write_file_mode == 2
                        out.SAVE_TIME = min(tile.FORCING.PARA.end_time, datenum([out.PARA.save_date num2str(str2num(datestr(out.SAVE_TIME,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
                        out.TEMP.write_out_final = 1;
                    end
                end

            end
        end

        function out = out2file_CG(out, tile)

            out = out2file_CG@OUT_BASE_OBSERVABLES(out, tile);
            
            out.STATVAR.timestamp = [];
            out.STATVAR.thawdepth = [];

        end
       
 
        
    end
end
