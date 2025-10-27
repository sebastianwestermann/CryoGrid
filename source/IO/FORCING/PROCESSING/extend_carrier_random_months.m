%========================================================================
% CryoGrid FORCING processing class 
%
% Authors:
% S. Westermann, December 2023
%
%========================================================================

classdef extend_carrier_random_months < process_BASE
    

    methods
        function proc = provide_PARA(proc)
            proc.PARA.selection_period_before = []; %vector with two entries, startyear and endyear
            proc.PARA.selection_period_after = [];
            
        end
        
        
        function proc = provide_CONST(proc)

        end
        
        
        function proc = provide_STATVAR(proc)
            
        end
        
        
        function proc = finalize_init(proc, tile)

        end
        
        
        function forcing = process(proc, forcing, tile)
            timestep = forcing.CARRIER.DATA.timeForcing(2,1)-forcing.CARRIER.DATA.timeForcing(1,1);
            values_per_day = floor(1./timestep);
            
            if ~any(strcmp(fieldnames(forcing.PARA), 'start_time')) || isempty(forcing.PARA.start_time) || isnan(forcing.PARA.start_time)
                start_date = datenum(year(forcing.REFERENCE.DATA.timeForcing(1,1)),1,1);
            else
                start_date = datenum(year(forcing.PARA.start_time),1,1);
                
            end
            if ~any(strcmp(fieldnames(forcing.PARA), 'end_time')) || isempty(forcing.PARA.end_time) || isnan(forcing.PARA.end_time)
                end_date = datenum(year(forcing.REFERENCE.DATA.timeForcing(end,1))+1,1,1);
            else
                end_date = datenum(year(forcing.PARA.end_time)+1,1,1);
            end
            
            %BEFORE
            variables = fieldnames(forcing.CARRIER.DATA);
            for i=1:size(variables,1)
                forcing_before.(variables{i,1}) = [];
            end
            forcing_before.new_timeForcing = [];
            if start_date < forcing.CARRIER.DATA.timeForcing(1,1)
                y = year(start_date);
                m = month(start_date);
                
                while datenum(y,m,1) < min(end_date, forcing.CARRIER.DATA.timeForcing(1,1))
                    timeForcing = datenum(y,m,1):timestep:datenum(y,m+1,1)-timestep;
                    rand_year = proc.PARA.selection_period_before(1) + floor(rand(1).*(proc.PARA.selection_period_before(2)-proc.PARA.selection_period_before(1)));
                    timeForcing = [datenum(y,m,1):timestep:datenum(y,m+1,1)-timestep]';
                    start_id = find(abs(forcing.CARRIER.DATA.timeForcing(:,1)-datenum(rand_year,m,1)) < timestep/100);
                    end_id = start_id + size(timeForcing,1)-1;
                    if end_id > size(forcing.CARRIER.DATA.timeForcing,1)
                        end_id = size(forcing.CARRIER.DATA.timeForcing,1);
                        start_id = end_id - size(timeForcing,1) + 1;
                    end
                    
                    for i=1:size(variables,1)
                        if any(strcmp(fieldnames(forcing.CARRIER.DATA), variables{i,1}))
                            forcing_before.(variables{i,1}) = [forcing_before.(variables{i,1}); forcing.CARRIER.DATA.(variables{i,1})(start_id:end_id,1)];
                        end
                    end
                    forcing_before.new_timeForcing = [forcing_before.new_timeForcing; timeForcing];
                    
                    m=m+1;
                    if m==13
                        y=y+1;
                        m=1;
                    end
                end
            end
            
            
            %AFTER
            for i=1:size(variables,1)
                forcing_after.(variables{i,1}) = [];
            end
            forcing_after.new_timeForcing = [];
            
            if forcing.CARRIER.DATA.timeForcing(end,1) < end_date
                y = year(forcing.CARRIER.DATA.timeForcing(end,1));
                m = month(forcing.CARRIER.DATA.timeForcing(end,1));
                
                while datenum(y,m,1) < end_date
                    timeForcing = datenum(y,m,1):timestep:datenum(y,m+1,1)-timestep;
                    rand_year = proc.PARA.selection_period_after(1) + floor(rand(1).*(proc.PARA.selection_period_after(2)-proc.PARA.selection_period_after(1)));
                    timeForcing = [datenum(y,m,1):timestep:datenum(y,m+1,1)-timestep]';
                    start_id = find(abs(forcing.CARRIER.DATA.timeForcing(:,1)-datenum(rand_year,m,1)) < timestep/100);
                    end_id = start_id + size(timeForcing,1)-1;
                    if end_id > size(forcing.CARRIER.DATA.timeForcing,1)
                        end_id = size(forcing.CARRIER.DATA.timeForcing,1);
                        start_id = end_id - size(timeForcing,1) + 1;
                    end
                    
                    for i=1:size(variables,1)
                        if any(strcmp(fieldnames(forcing.CARRIER.DATA), variables{i,1}))
                            forcing_after.(variables{i,1}) = [forcing_after.(variables{i,1}); forcing.CARRIER.DATA.(variables{i,1})(start_id:end_id,1)];
                        end
                    end
                    forcing_after.new_timeForcing = [forcing_after.new_timeForcing; timeForcing];
                    
                    m=m+1;
                    if m==13
                        y=y+1;
                        m=1;
                    end
                end
            end
            forcing_before.true_timeForcing = forcing_before.timeForcing;
            forcing_before.timeForcing = forcing_before.new_timeForcing;
            forcing_after.true_timeForcing = forcing_after.timeForcing;
            forcing_after.timeForcing = forcing_after.new_timeForcing;
            
            %concatenate arrays
            if end_date + timestep < forcing.CARRIER.DATA.timeForcing(1,1)
                for i=1:size(variables,1)
                    forcing.CARRIER.DATA.(variables{i,1}) = [ forcing_before.(variables{i,1})];
                end
            elseif start_date - timestep > forcing.CARRIER.DATA.timeForcing(end,1)
                for i=1:size(variables,1)
                    forcing.CARRIER.DATA.(variables{i,1}) = [ forcing_after.(variables{i,1})];
                end
            else
                if ~isempty(forcing_before.timeForcing)
                    forcing.CARRIER.TEMP.start_time_original_forcing = forcing.CARRIER.DATA.timeForcing(1,1);
                    forcing.CARRIER.TEMP.end_time_original_forcing = forcing.CARRIER.DATA.timeForcing(end,1);
                    end_pos = find(abs(forcing_before.timeForcing(:,1) - (forcing.CARRIER.DATA.timeForcing(1,1)-timestep)) <timestep/100);
                    forcing.CARRIER.TEMP.start_id_original_forcing = end_pos+1;
                    forcing.CARRIER.DATA.true_timeForcing = [forcing_before.true_timeForcing(1:end_pos,1); forcing.CARRIER.DATA.timeForcing];
                    for i=1:size(variables,1)
                        forcing.CARRIER.DATA.(variables{i,1}) = [forcing_before.(variables{i,1})(1:end_pos,1); forcing.CARRIER.DATA.(variables{i,1})];
                    end
                end
                if ~isempty(forcing_after.timeForcing)
                    start_pos = find(abs(forcing_after.timeForcing(:,1) - (forcing.CARRIER.DATA.timeForcing(end,1) + timestep)) <timestep/100);
                    forcing.CARRIER.TEMP.end_id_original_forcing = size(forcing.CARRIER.DATA.timeForcing,1);
                    forcing.CARRIER.DATA.true_timeForcing = [forcing.CARRIER.DATA.timeForcing; forcing_after.true_timeForcing(start_pos:end,1)];
                    for i=1:size(variables,1)
                        forcing.CARRIER.DATA.(variables{i,1}) = [forcing.CARRIER.DATA.(variables{i,1}); forcing_after.(variables{i,1})(start_pos:end,1)];
                    end
                end
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

