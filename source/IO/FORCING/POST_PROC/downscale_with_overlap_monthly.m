%========================================================================
% CryoGrid FORCING post-processing class 

%interfaces with files downloaded with "download_CMIP" from the Copernicus
%data portal
% Authors:
% S. Westermann, July 2023
%
%========================================================================

classdef downscale_with_overlap_monthly < bias_correct_BASE 
    
    properties
        
    end
    
    methods
        function post_proc = provide_PARA(post_proc)
            
            post_proc.PARA.forcing_class = []; 
            post_proc.PARA.focing_class_index = [];
            post_proc.PARA.start_overlap_period = [];
            post_proc.PARA.end_overlap_period = [];
            
            post_proc.PARA.all_snow_T = [];
            post_proc.PARA.all_rain_T = [];
             
        end
        
        
        function post_proc = provide_CONST(post_proc)
            post_proc.CONST.L_f = []; 
            post_proc.CONST.sigma = [];
            post_proc.CONST.day_sec = [];
            post_proc.CONST.Tmfw = [];
        end
        
        
        function post_proc = provide_STATVAR(post_proc)
            
        end
        
        
        function post_proc = finalize_init(post_proc, tile)
            
        end
        
        
        function forcing = post_process(post_proc, forcing, tile)
            
           
        end
        
        
%                 %-------------param file generation-----
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