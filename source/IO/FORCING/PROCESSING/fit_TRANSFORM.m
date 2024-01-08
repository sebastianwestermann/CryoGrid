%========================================================================
% CryoGrid FORCING processing class
%
%
% Authors:
% S. Westermann, December 2022
%
%========================================================================

classdef fit_TRANSFORM < matlab.mixin.Copyable %makes the TRANSFORM object
    
    properties
        CONST
        PARA
        STATVAR
        TEMP
    end
    
    methods
        function proc = provide_PARA(proc)
            proc.PARA.start_overlap = [];
            proc.PARA.end_overlap = [];            
            proc.PARA.variables = [];
            proc.PARA.overlap_target_interval = []; %interval for which the list of overalp-pairs is created
            proc.PARA.transform_class = [];
            proc.PARA.transform_class_index = [];
            
%             proc.PARA.display = 0;

        end
        
        
        function proc = provide_CONST(proc)

        end
        
        
        function proc = provide_STATVAR(proc)
            
        end
        
        
        function proc = finalize_init(proc, tile)

        end
        
        
        function forcing = process(proc, forcing, tile)
            reference_class = forcing.REFERENCE;
            carrier_class = forcing.CARRIER;
            
            %find indices for overlap period
            if ~isempty(proc.PARA.start_overlap) && sum(isnan(proc.PARA.start_overlap))==0
                proc.PARA.start_overlap = datenum(proc.PARA.start_overlap(1,1), proc.PARA.start_overlap(2,1), proc.PARA.start_overlap(3,1));
                start_index_carrier = find(carrier_class.DATA.timeForcing(:,1) >= proc.PARA.start_overlap - (carrier_class.DATA.timeForcing(2,1)-carrier_class.DATA.timeForcing(1,1))./2,1); %subtract half a timestep
                start_index_reference = find(reference_class.DATA.timeForcing(:,1) >= proc.PARA.start_overlap - (reference_class.DATA.timeForcing(2,1)-reference_class.DATA.timeForcing(1,1))./2,1); %subtract half a timestep
            else %determine which class starts first
                if carrier_class.DATA.timeForcing(1,1) - (carrier_class.DATA.timeForcing(2,1)-carrier_class.DATA.timeForcing(1,1))./2 > reference_class.DATA.timeForcing(1,1) - (reference_class.DATA.timeForcing(2,1)-reference_class.DATA.timeForcing(1,1))./2
                    start_index_carrier = 1;
                    start_index_reference = find(reference_class.DATA.timeForcing(:,1) >= carrier_class.DATA.timeForcing(1,1) - (carrier_class.DATA.timeForcing(2,1)-carrier_class.DATA.timeForcing(1,1))./2 - ...
                        (reference_class.DATA.timeForcing(2,1)-reference_class.DATA.timeForcing(1,1))./2,1); %subtract half a timestep
                    
                else
                    start_index_carrier = find(carrier_class.DATA.timeForcing(:,1) >= reference_class.DATA.timeForcing(1,1) - (carrier_class.DATA.timeForcing(2,1)-carrier_class.DATA.timeForcing(1,1))./2 - ...
                        (reference_class.DATA.timeForcing(2,1)-reference_class.DATA.timeForcing(1,1))./2,1); %subtract half a timestep
                    start_index_reference = 1;
                end
                
            end
            if ~isempty(proc.PARA.end_overlap) && sum(isnan(proc.PARA.end_overlap))==0
                proc.PARA.end_overlap = datenum(proc.PARA.end_overlap(1,1), proc.PARA.end_overlap(2,1), proc.PARA.end_overlap(3,1));
                end_index_carrier = find(carrier_class.DATA.timeForcing(:,1) < proc.PARA.end_overlap - (carrier_class.DATA.timeForcing(2,1)-carrier_class.DATA.timeForcing(1,1))./2,1, 'last'); %subtract half a timestep
                end_index_reference = find(reference_class.DATA.timeForcing(:,1) < proc.PARA.end_overlap - (reference_class.DATA.timeForcing(2,1)-reference_class.DATA.timeForcing(1,1))./2,1, 'last'); %subtract half a timestep
            else
                if carrier_class.DATA.timeForcing(end,1) + (carrier_class.DATA.timeForcing(2,1)-carrier_class.DATA.timeForcing(1,1))./2 > reference_class.DATA.timeForcing(end,1) + (reference_class.DATA.timeForcing(2,1)-reference_class.DATA.timeForcing(1,1))./2
                    end_index_carrier = find(carrier_class.DATA.timeForcing(:,1) < reference_class.DATA.timeForcing(end,1) + (carrier_class.DATA.timeForcing(2,1)-carrier_class.DATA.timeForcing(1,1))./2 + ...
                        (reference_class.DATA.timeForcing(2,1)-reference_class.DATA.timeForcing(1,1))./2,1, 'last'); %subtract half a timestep
                    end_index_reference = size(reference_class.DATA.timeForcing,1);
                else
                    end_index_carrier = size(carrier_class.DATA.timeForcing,1);
                    end_index_reference = find(reference_class.DATA.timeForcing(:,1) < carrier_class.DATA.timeForcing(end,1) + (carrier_class.DATA.timeForcing(2,1)-carrier_class.DATA.timeForcing(1,1))./2 + ...
                        (reference_class.DATA.timeForcing(2,1)-reference_class.DATA.timeForcing(1,1))./2, 1, 'last'); %subtract half a timestep
                end
            end
            
%             overlap.time = {};
%             overlap.values= {}; %first argument is carrier, second is reference
            
            for i=1:size(proc.PARA.variables,1)
                
                overlap_pairs_time = [];
                overlap_pairs = [];
                if strcmp(proc.PARA.overlap_target_interval{i,1}, 'full')
                    
                    for j=start_index_carrier:end_index_carrier
                        pos_ref = find(abs(carrier_class.DATA.timeForcing(j,1) - reference_class.DATA.timeForcing(:,1))<1/24/30);
                        if ~isempty(pos_ref)
                            overlap_pairs_time = [overlap_pairs_time; carrier_class.DATA.timeForcing(j,1)];
                            overlap_pairs = [overlap_pairs; [carrier_class.DATA.(proc.PARA.variables{i,1})(j,1) reference_class.DATA.(proc.PARA.variables{i,1})(pos_ref,1)]];
                        end
                    end
                    
                elseif  strcmp(proc.PARA.overlap_target_interval{i,1}, 'day')
                    
                    for j=floor(carrier_class.DATA.timeForcing(start_index_carrier,1)):floor(carrier_class.DATA.timeForcing(end_index_carrier,1))
                        pos_carrier = find(carrier_class.DATA.timeForcing(:,1) >= j & carrier_class.DATA.timeForcing(:,1) < j+1);
                        pos_ref = find(reference_class.DATA.timeForcing(:,1) >= j & reference_class.DATA.timeForcing(:,1) < j+1);
                        if ~isempty(pos_ref) && ~isempty(pos_carrier)
                            overlap_pairs_time = [overlap_pairs_time; j];
                            overlap_pairs = [overlap_pairs; [mean(carrier_class.DATA.(proc.PARA.variables{i,1})(pos_carrier,1)) mean(reference_class.DATA.(proc.PARA.variables{i,1})(pos_ref,1))]];
                        end
                    end
                    
                    
                elseif strcmp(proc.PARA.overlap_target_interval{i,1}, 'month')
                    %compile overlap pairs
                    
                    y = year(max(carrier_class.DATA.timeForcing(start_index_carrier,1), reference_class.DATA.timeForcing(start_index_reference,1)));
                    m = month(max(carrier_class.DATA.timeForcing(start_index_carrier,1), reference_class.DATA.timeForcing(start_index_reference,1)));
                    end_y=year(min(carrier_class.DATA.timeForcing(end_index_carrier), reference_class.DATA.timeForcing(end_index_reference)));
                    end_m=month(min(carrier_class.DATA.timeForcing(end_index_carrier), reference_class.DATA.timeForcing(end_index_reference)));
                    while datenum(y,m,1) <= datenum(end_y, end_m, 1)
                        overlap_pairs_time = [overlap_pairs_time; datenum(y,m,15)];
                        pos_carrier = find(carrier_class.DATA.timeForcing(:,1)>=datenum(y,m,1) & carrier_class.DATA.timeForcing(:,1)<datenum(y,m+1,1));
                        pos_reference = find(reference_class.DATA.timeForcing(:,1)>=datenum(y,m,1) & reference_class.DATA.timeForcing(:,1)<datenum(y,m+1,1));
                        overlap_pairs = [overlap_pairs; mean(carrier_class.DATA.(proc.PARA.variables{i,1})(pos_carrier,1)) ...
                            mean(reference_class.DATA.(proc.PARA.variables{i,1})(pos_reference,1))];
                        m=m+1;
                        if m==13
                            m=1;
                            y=y+1;
                        end
                    end
%                     reference_before_time = [];
%                     reference_before = [];
%                     y = year(reference_class.DATA.timeForcing(1,1));
%                     m = month(reference_class.DATA.timeForcing(1,1));
%                     end_y = year(max(carrier_class.DATA.timeForcing(start_index_carrier,1), reference_class.DATA.timeForcing(start_index_reference,1)));
%                     end_m = month(max(carrier_class.DATA.timeForcing(start_index_carrier,1), reference_class.DATA.timeForcing(start_index_reference,1)));
%                     while datenum(y,m,1) < datenum(end_y, end_m, 1)
%                         reference_before_time = [reference_before_time; datenum(y,m,15)];
%                         pos_reference = find(reference_class.DATA.timeForcing(:,1)>=datenum(y,m,1) & reference_class.DATA.timeForcing(:,1)<datenum(y,m+1,1));
%                         reference_before = [reference_before; mean(reference_class.DATA.(proc.PARA.variables{i,1})(pos_reference,1))];
%                         m=m+1;
%                         if m==13
%                             m=1;
%                             y=y+1;
%                         end
%                     end
%                     reference_after_time = [];
%                     reference_after = [];
%                     y=year(min(carrier_class.DATA.timeForcing(end_index_carrier), reference_class.DATA.timeForcing(end_index_reference)));
%                     m=month(min(carrier_class.DATA.timeForcing(end_index_carrier), reference_class.DATA.timeForcing(end_index_reference)));
%                     while datenum(y,m,1) < datenum(year(reference_class.DATA.timeForcing(end,1)), month(reference_class.DATA.timeForcing(end,1))+1, 1)
%                         reference_after_time = [reference_after_time; datenum(y,m,15)];
%                         pos_reference = find(reference_class.DATA.timeForcing(:,1)>=datenum(y,m,1) & reference_class.DATA.timeForcing(:,1)<datenum(y,m+1,1));
%                         reference_after = [reference_after; mean(reference_class.DATA.(proc.PARA.variables{i,1})(pos_reference,1))];
%                         m=m+1;
%                         if m==13
%                             m=1;
%                             y=y+1;
%                         end
%                     end
                elseif  strcmp(proc.PARA.overlap_target_interval{i,1}, 'year')
                    
                    for j=year(carrier_class.DATA.timeForcing(start_index_carrier,1)):year(carrier_class.DATA.timeForcing(end_index_carrier,1))
                        pos_carrier = find(year(carrier_class.DATA.timeForcing(:,1)) == j);
                        pos_ref = find(year(reference_class.DATA.timeForcing(:,1)) == j);
                        if ~isempty(pos_ref) && ~isempty(pos_carrier)
                            overlap_pairs_time = [overlap_pairs_time; j];
                            overlap_pairs = [overlap_pairs; [mean(carrier_class.DATA.(proc.PARA.variables{i,1})(pos_carrier,1)) mean(reference_class.DATA.(proc.PARA.variables{i,1})(pos_ref,1))]];
                        end
                    end
                else
                    overlap_pairs_time = NaN;
                    overlap_pairs = NaN;
                end
                
                transform_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(proc.PARA.transform_class{i,1}){proc.PARA.transform_class_index(i,1),1});
                transform_class = finalize_init(transform_class, tile);
                transform_class.PARA.variable = proc.PARA.variables{i,1};
                transform_class.PARA.overlap_target_interval = proc.PARA.overlap_target_interval{i,1};
                transform_class.TEMP.overlap_pairs_time = overlap_pairs_time; 
                transform_class.TEMP.overlap_pairs = overlap_pairs;
                
                transform_class = fit_transform(transform_class, forcing, tile);
                                
                forcing.TRANSFORM{i,1} = transform_class;
                 
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

