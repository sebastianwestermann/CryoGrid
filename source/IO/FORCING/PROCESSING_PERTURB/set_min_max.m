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

classdef set_min_max < FORCING_base 
    
    properties
        
    end
    
    methods
        function proc = provide_PARA(proc)
            
            proc.PARA.variables = [];  
            proc.PARA.minimum = [];
            proc.PARA.maximum = [];
        end
        
        
        function proc = provide_CONST(proc)

        end
        
        
        function proc = provide_STATVAR(proc)
            
        end
        
        
        function proc = finalize_init(proc, tile)
            if isempty(proc.PARA.minimum) || size(proc.PARA.minimum,1) < size(proc.PARA.variables,1)
                proc.PARA.minimum = repmat(-Inf, size(proc.PARA.variables,1),1);
            else
                proc.PARA.minimum(isnan(proc.PARA.minimum))=-Inf;
            end
            if isempty(proc.PARA.maximum) || size(proc.PARA.maximum,1) < size(proc.PARA.variables,1)
                proc.PARA.maximum = repmat(Inf, size(proc.PARA.variables,1),1);
            else
                proc.PARA.maximum(isnan(proc.PARA.maximum)) = Inf;
            end
        end
        
        
        function forcing = process(proc, forcing, tile)
            
            for i=1:size(proc.PARA.variables,1)
               forcing.DATA.(proc.PARA.variables{i,1}) = max(min(forcing.DATA.(proc.PARA.variables{i,1}), proc.PARA.maximum(i,1)), proc.PARA.minimum(i,1));
            end
        
        end
        
        
        %----------------------------------------------
        %perturb
        function proc = preprocess(proc, forcing, tile)
            
        end
        
        function forcing = perturb_forcing(proc, forcing, tile)
            
            for i=1:size(proc.PARA.variables,1)
               forcing.TEMP.(proc.PARA.variables{i,1}) = max(min(forcing.TEMP.(proc.PARA.variables{i,1}), proc.PARA.maximum(i,1)), proc.PARA.minimum(i,1));
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