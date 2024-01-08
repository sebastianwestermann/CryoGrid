%========================================================================
% CryoGrid FORCING processing class detrend_reference_and_carrier
%
% Authors:
% S. Westermann, December 2023
%
%========================================================================

classdef retrend_forcing_monthly < process_BASE
    

    methods
        function proc = provide_PARA(proc)
            proc.PARA.variables = [];
            proc.PARA.relative_correction = [];
        end
        
        
        function proc = provide_CONST(proc)

        end
        
        
        function proc = provide_STATVAR(proc)
            
        end
        
        
        function proc = finalize_init(proc, tile)

        end
        
        
        function forcing = process(proc, forcing, tile)
            for i=1:size(proc.PARA.variables,1)
                
                c_points = forcing.CARRIER.STATVAR.(proc.PARA.variables{i,1}).variable_points;
                r_points = forcing.REFERENCE.STATVAR.(proc.PARA.variables{i,1}).variable_points;
                
                t_points = r_points;
                t_points(2:3,:) = c_points;
                t_points(1,:) = t_points(2,:) + r_points(1,:) - r_points(2,:);
                t_points(4,:) = t_points(1,:) - r_points(1,:) + r_points(4,:);
                
                year_points = [year(forcing.DATA.timeForcing(1,1)):year(forcing.DATA.timeForcing(end,1))]';
                
                offset=[];
                for m=1:12
                    offset = [offset piece_wise_linear(t_points(:,m), year_points, forcing.REFERENCE.STATVAR.(proc.PARA.variables{i,1}).year_points(:,m))];
                end
                if ~proc.PARA.relative_correction(i,1)
                    offset = offset - repmat(forcing.CARRIER.STATVAR.(proc.PARA.variables{i,1}).variable_points(1,:), size(offset,1),1);
                end

                
                y = year(forcing.DATA.timeForcing(1,1));
                m = month(forcing.DATA.timeForcing(1,1));
                year_count = 1;
                while datenum(y,m,1) <= datenum(year(forcing.DATA.timeForcing(end,1)), month(forcing.DATA.timeForcing(end,1)),1)
                    if ~proc.PARA.relative_correction(i,1)
                        forcing.DATA.(proc.PARA.variables{i,1})(forcing.DATA.timeForcing>=datenum(y,m,1) & forcing.DATA.timeForcing<datenum(y,m+1,1),1) = ...
                            forcing.DATA.(proc.PARA.variables{i,1})(forcing.DATA.timeForcing>=datenum(y,m,1) & forcing.DATA.timeForcing<datenum(y,m+1,1),1) + offset(year_count,m);
                    else
                        forcing.DATA.(proc.PARA.variables{i,1})(forcing.DATA.timeForcing>=datenum(y,m,1) & forcing.DATA.timeForcing<datenum(y,m+1,1),1) = ...
                            forcing.DATA.(proc.PARA.variables{i,1})(forcing.DATA.timeForcing>=datenum(y,m,1) & forcing.DATA.timeForcing<datenum(y,m+1,1),1) .* offset(year_count,m);
                    end
                    m=m+1;
                    if m==13
                        y=y+1;
                        year_count = year_count+1;
                        m=1;
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

