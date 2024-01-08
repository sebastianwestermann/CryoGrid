%========================================================================
% CryoGrid FORCING processing class
%
%
% Authors:
% S. Westermann, December 2022
%
%========================================================================

classdef initialize_TEMP < matlab.mixin.Copyable 
    
    properties
        CONST
        PARA
        STATVAR
        TEMP
    end
    
    methods
        function proc = provide_PARA(proc)
            proc.PARA.is_slope = [];  %rainfall fraction assumed in sumulations (rainfall from the forcing data file is multiplied by this parameter)
        end
        
        
        function proc = provide_CONST(proc)

        end
        
        
        function proc = provide_STATVAR(proc)
            
        end
        
        
        function proc = finalize_init(proc, tile)

        end
        
        
        function forcing = process(proc, forcing, tile)
            
            forcing.TEMP.snowfall=0;
            forcing.TEMP.rainfall=0;
            forcing.TEMP.Lin=0;
            forcing.TEMP.Sin=0;
            forcing.TEMP.Tair=0;
            forcing.TEMP.wind=0;
            forcing.TEMP.q=0;
            forcing.TEMP.p=0;
            
            if proc.PARA.is_slope
                forcing.TEMP.Sin_dif = 0;
                forcing.TEMP.Sin_dir = 0;
                forcing.TEMP.sunElevation = 0;
                forcing.TEMP.azimuth = 0;
            end
        end
        
        
                %-------------param file generation-----
%         function post_proc = param_file_info(post_proc)
%             post_proc = provide_PARA(post_proc);
% 
%             post_proc.PARA.STATVAR = [];
%             post_proc.PARA.class_category = 'FORCING POST_PROCESSING';
%             post_proc.PARA.options = [];
%             
%             post_proc.PARA.eliminate_fraction = [];
%             post_proc.PARA.survive_fraction = [];
%                         
%             post_proc.PARA.default_value.window_size = {7};
%             post_proc.PARA.comment.window_size = {'window size in days within which precipitation is reallocated'};
%             
%             post_proc.PARA.default_value.eliminate_fraction = {0.5};
%             post_proc.PARA.comment.eliminate_fraction = {'fraction of smallest precipitation events (= timestamps with precipitation) that is reallocated to larger events'};
%             
%             post_proc.PARA.default_value.survive_fraction = {0.5};  
%             post_proc.PARA.comment.survive_fraction = {'fraction of largest precipitation events (= timestamps with precipitation) that the small events are reallocated to'};
%             
%         end
        
    end
    
end

