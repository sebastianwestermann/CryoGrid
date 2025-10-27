%========================================================================
% CryoGrid FORCING post-processing class condense_precip
%
% The class changes the time distribution of precipitation (both rain- and
% snowfall) by moving the precipitation from small events to large events.
% 
% It is recommended to compare the resulting precipitation statistics to
% measurements of other data sources.
%
% Authors:
% S. Westermann, December 2022
%
%========================================================================

classdef set_min_RH < FORCING_base 
    
    properties
        
    end
    
    methods
        function proc = provide_PARA(proc)
            
            proc.PARA.min_RH = [];  

        end
        
        
        function proc = provide_CONST(proc)

        end
        
        
        function proc = provide_STATVAR(proc)
            
        end
        
        
        function proc = finalize_init(proc, tile)

        end
        
        
        function forcing = process(proc, forcing, tile)
            
            A1=17.62; B1=243.12; C1=611.2;
            Magnus=@(tc) double(tc>=0) .* C1.*exp(A1.*tc./(B1+tc)) + double(tc<0) .* C1.*exp(22.46.*tc./(272.61+tc)); % A version of the Magnus formula with the AERK parameters.
            RH = forcing.DATA.q ./0.622 ./ (Magnus(forcing.DATA.Tair)./forcing.DATA.p);
            
            RH=max(proc.PARA.min_RH, RH);
            forcing.DATA.q = 0.622 .* RH .* (Magnus(forcing.DATA.Tair)./forcing.DATA.p);
        
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