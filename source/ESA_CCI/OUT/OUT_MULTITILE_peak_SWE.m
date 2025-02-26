classdef OUT_MULTITILE_peak_SWE < matlab.mixin.Copyable


    properties
        CONST
        PARA
        STATVAR
        TEMP
        OUTPUT_TIME
        
    end
    
    methods
        
        
        function out = provide_PARA(out)
            out.PARA.out_folder = [];
            out.PARA.output_timestep = [];

        end
        
        function out = provide_CONST(out)
        end
        
        function out = provide_STATVAR(out)

        end
        
        
        function out = finalize_init(out, tile)
                      
            if ~(exist([out.PARA.out_folder tile.RUN_INFO.PARA.run_name])==7)
                mkdir([out.PARA.out_folder tile.RUN_INFO.PARA.run_name]);
            end
            out.PARA.out_folder = [out.PARA.out_folder tile.RUN_INFO.PARA.run_name '/'];

            out.TEMP.peak_SWE = 0;
            out.OUTPUT_TIME = tile.FORCING.PARA.start_time + out.PARA.output_timestep;
            out.STATVAR.SWE = [];
            out.STATVAR.SWE_std = [];
            out.STATVAR.timestamp = [];
            out.STATVAR.da_param = {};

            % out.SAVE_TIME = Inf;
        end
        
        function out = store_OUT(out, tile) 
            if tile.t>=out.OUTPUT_TIME
                disp(datestr(tile.t))
                SWE = reshape(tile.SUBSURFACE_CLASS.STATVAR.SWE, tile.PARA.number_of_realizations, tile.ENSEMBLE.PARA.subgrid_ensemble_size, tile.ENSEMBLE.PARA.grid_ensemble_size);
                SWE=squeeze(mean(SWE,2));
                out.TEMP.peak_SWE = max(out.TEMP.peak_SWE, SWE);

                out.OUTPUT_TIME = tile.t + out.PARA.output_timestep;

                % if tile.t >= out.SAVE_TIME
                %     out.TEMP.peak_SWE = 0;
                %     out.SAVE_TIME = Inf;
                % 
                % end
            end
        end

        function out = save_out_ensemble(out, tile, weights, best_parameters)
            SWE = sum(out.TEMP.peak_SWE .* weights, 2);
            SWE_std = sqrt(sum((out.TEMP.peak_SWE - SWE).^2 .* weights,2));
            out.STATVAR.SWE = [out.STATVAR.SWE SWE];
            out.STATVAR.SWE_std = [out.STATVAR.SWE_std SWE_std];
            out.STATVAR.timestamp = [out.STATVAR.timestamp; tile.t];
            out.STATVAR.da_param = [out.STATVAR.da_param; best_parameters];
            result = out.STATVAR;
            save([out.PARA.out_folder 'out_SWE' num2str(tile.PARA.range(1)) '_' num2str(tile.PARA.range(end)) '.mat'], 'result')
            out = reset_out(out, tile);
        end

        function out = reset_out(out, tile)
            out.TEMP.peak_SWE = 0;
            out.OUTPUT_TIME = tile.t + out.PARA.output_timestep;
            % out.SAVE_TIME = Inf;
        end


    end
end

