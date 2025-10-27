%========================================================================
% CryoGrid FORCING processing class
%
%
% Authors:
% S. Westermann, December 2022
%
%========================================================================

classdef set_start_end_time_nc < matlab.mixin.Copyable 
    
    properties
        CONST
        PARA
        STATVAR
        TEMP
    end
    
    methods
        function proc = provide_PARA(proc)
            proc.PARA.start_time = []; % start time of the simulations (must be within the range of data in forcing file)
            proc.PARA.end_time = [];   % end time of the simulations (must be within the range of data in forcing file)
        end
        
        
        function proc = provide_CONST(proc)

        end
        
        
        function proc = provide_STATVAR(proc)
            
        end
        
        
        function proc = finalize_init(proc, tile)

        end
        
        
        function forcing = process(proc, forcing, tile)
            
            forcing.PARA.start_time = datenum(proc.PARA.start_time(1,1), proc.PARA.start_time(2,1), proc.PARA.start_time(3,1));
            forcing.PARA.end_time = datenum(proc.PARA.end_time(1,1), proc.PARA.end_time(2,1),proc.PARA.end_time(3,1));
            forcing.TEMP.current_month = month(forcing.PARA.start_time);
            forcing.TEMP.current_quarter = floor((forcing.TEMP.current_month-1)./3)+1;
            forcing.TEMP.current_year = year(forcing.PARA.start_time);

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

