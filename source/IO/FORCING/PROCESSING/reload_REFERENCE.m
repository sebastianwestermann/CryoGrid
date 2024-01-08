%========================================================================
% CryoGrid FORCING processing class
%
%
% Authors:
% S. Westermann, December 2022
%
%========================================================================

classdef reload_REFERENCE < matlab.mixin.Copyable %makes the TRANSFORM object
    
    properties
        CONST
        PARA
        STATVAR
        TEMP
    end
    
    methods
        function proc = provide_PARA(proc)
            proc.PARA.reference_forcing_class = [];
            proc.PARA.reference_forcing_class_index = [];
            proc.PARA.offset_from_GMT_reference = [];
        end
        
        
        function proc = provide_CONST(proc)

        end
        
        
        function proc = provide_STATVAR(proc)
            
        end
        
        
        function proc = finalize_init(proc, tile)

        end
        
        
        function forcing = process(proc, forcing, tile)
            if ~isempty(proc.PARA.reference_forcing_class) && sum(isnan(proc.PARA.reference_forcing_class))==0
                reference_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(proc.PARA.reference_forcing_class){proc.PARA.reference_forcing_class_index,1});
                reference_class = finalize_init(reference_class, tile);
                reference_class.DATA.timeForcing = reference_class.DATA.timeForcing - proc.PARA.offset_from_GMT_reference ./ 24;
                forcing.REFERENCE = reference_class;
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

