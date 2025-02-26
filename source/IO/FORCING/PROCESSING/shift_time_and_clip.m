%========================================================================
% CryoGrid FORCING processing class
%
%
% Authors:
% S. Westermann, December 2022
%
%========================================================================

classdef shift_time_and_clip < process_BASE 
    
    properties

    end
    
    methods
        function proc = provide_PARA(proc)
            proc.PARA.start_time_slice = [];
            proc.PARA.end_time_slice = [];
        end
        
        
        function proc = provide_CONST(proc)

        end
        
        
        function proc = provide_STATVAR(proc)
            
        end
        
        
        function proc = finalize_init(proc, tile)
            if ~isempty(proc.PARA.start_time_slice) && sum(isnan(proc.PARA.start_time_slice)) == 0
                proc.PARA.start_time_slice = datenum(proc.PARA.start_time_slice(1), proc.PARA.start_time_slice(2), proc.PARA.start_time_slice(3));
            else
                proc.PARA.start_time_slice = [];
            end
            if ~isempty(proc.PARA.end_time_slice) && sum(isnan(proc.PARA.end_time_slice)) == 0
                proc.PARA.end_time = datenum(proc.PARA.end_time_slice(1), proc.PARA.end_time_slice(2), proc.PARA.end_time_slice(3));
            else
                proc.PARA.end_time_slice = [];
            end
            if isempty(proc.PARA.start_time_slice) && isempty(proc.PARA.end_time_slice)
                proc.PARA.start_time_slice = forcing.DATA.timeForcing(1,1);
            end

        end
        
        
        function forcing = process(proc, forcing, tile)
            number_of_days = forcing.PARA.end_time - forcing.PARA.start_time; %100 - 0
            if ~isempty(proc.PARA.start_time_slice)
                proc.PARA.end_time_slice = proc.PARA.start_time_slice + number_of_days; %1850 +100 =1950
            else
                proc.PARA.start_time_slice = proc.PARA.end_time_slice - number_of_days; 
            end
            time_offset = forcing.PARA.start_time - proc.PARA.start_time_slice; % 0 - 1850 = -1850 

            forcing.DATA.timeForcing = forcing.DATA.timeForcing + time_offset; % 1850- 1850

            range = find(forcing.DATA.timeForcing>=forcing.PARA.start_time & forcing.DATA.timeForcing <= forcing.PARA.end_time);
            variable_list = fieldnames(forcing.DATA);
            for i=1:size(variable_list,1)
                forcing.DATA.(variable_list{i,1}) = forcing.DATA.(variable_list{i,1})(range,1);
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

