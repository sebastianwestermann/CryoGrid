%========================================================================
% CryoGrid OUT class OUT_do_nothing
% does not store any output, generally used during model spin-up
% S. Westermann, Jan 2021
%========================================================================


classdef OUT_do_nothing < OUT_BASE
 

    properties
		
	end
    
    
    methods
		
        
        function out = provide_PARA(out)         
            out.PARA.display_timestep = [];
        end

		
		function out = finalize_init(out, tile)

            if ~isempty(out.PARA.display_timestep)
                out.OUTPUT_TIME = tile.FORCING.PARA.start_time + out.PARA.display_timestep;
            else
                out.OUTPUT_TIME = tile.FORCING.PARA.end_time +1; %set to time that is never reached
            end
        end
        
            
        function out = store_OUT(out, tile)           
            if tile.t >= out.OUTPUT_TIME 
                disp([datestr(tile.t)])
                out.OUTPUT_TIME = out.OUTPUT_TIME + out.PARA.display_timestep;
            end
            
        end

                %-------------param file generation-----
        function out = param_file_info(out)
            out = provide_PARA(out);

            out.PARA.STATVAR = [];
            out.PARA.options = [];
            out.PARA.class_category = 'OUT';
           
            out.PARA.default_value.display_timestep = {5};
            out.PARA.comment.display_timestep = {'timestep that model progress is displayed [days]'};

        end

    end
end