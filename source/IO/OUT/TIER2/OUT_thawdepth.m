classdef OUT_thawdepth < OUT_BASE
    
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
            out.STATVAR.thawdepth = [];

            out.TEMP.keyword = 'thawdepth';
            out.TEMP.output_var = 'CG_thawdepth';

        end


        function out = store_OUT(out, tile)

            if tile.t >= out.OUTPUT_TIME

                disp([datestr(tile.t)])

                out = state2out(out, tile);
                out.OUTPUT_TIME = out.OUTPUT_TIME + out.PARA.output_timestep;

                if tile.t >= out.SAVE_TIME
                    out = out2file_CG(out, tile);
                    out.STATVAR.timestamp = [];
                    out.STATVAR.thawdepth = [];
                end

            end
        end
        
        function out = state2out(out, tile)
                result = 0;
                CURRENT = tile.TOP.NEXT;
                while ~is_ground_surface(CURRENT) && ~(strcmp(class(CURRENT), 'Bottom'))
                    CURRENT = CURRENT.NEXT;
                end
                i=1;
                while CURRENT.STATVAR.T(i) >= 0
                    result = result + CURRENT.STATVAR.layerThick(i);
                    i=i+1;
                end
                out.STATVAR.thawdepth = [out.STATVAR.thawdepth; result];
                out.STATVAR.timestamp = [out.STATVAR.timestamp; tile.t];
        end
       
 
        
    end
end
