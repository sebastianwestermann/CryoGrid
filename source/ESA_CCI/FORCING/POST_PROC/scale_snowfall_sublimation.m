%========================================================================
% CryoGrid FORCING post-processing 
%
%
% Authors:
% S. Westermann, January 2023
%
%========================================================================

classdef scale_snowfall_sublimation < FORCING_base
    
    properties
        
    end
    
    methods
        function post_proc = provide_PARA(post_proc)
            
            post_proc.PARA.start_date = [];
            
        end
        
        
        function post_proc = provide_CONST(post_proc)

        end
        
        
        function post_proc = provide_STATVAR(post_proc)
            
        end
        
        
        function post_proc = finalize_init(post_proc, tile)

        end
        
        
        function forcing = process(post_proc, forcing, tile)
            % a=squeeze(forcing.DATA.ERA_sublimation_downscaled(2,:,:));
            % plot(a(:));
            % hold on
            total_range = logical(repmat(forcing.DATA.timestamp.*0, tile.PARA.number_of_realizations,1,1));
            snow_fraction_av = zeros(tile.PARA.number_of_realizations,1);
            sublimation_factor_av = zeros(tile.PARA.number_of_realizations,1);
            count = 0;
            for i=1:size(post_proc.PARA.scale_parameters, 1)
                count = count + 1;
                timestamp = repmat(forcing.DATA.timestamp, tile.PARA.number_of_realizations, 1, 1);
                range = (timestamp >= datenum([post_proc.PARA.start_date num2str(post_proc.PARA.scale_parameters{i,1}.year)], 'dd.mm.yyyy') & ...
                    timestamp < datenum([post_proc.PARA.start_date num2str(post_proc.PARA.scale_parameters{i,1}.year+1)], 'dd.mm.yyyy'));
                total_range = total_range + range;
                snow_fraction = post_proc.PARA.scale_parameters{i,1}.snow_fraction;
                snow_fraction(isnan(snow_fraction))=1;
                snow_fraction_av = snow_fraction_av + snow_fraction;
                snow_fraction = repmat(snow_fraction, 1, size(timestamp,2), size(timestamp,3));
                forcing.DATA.ERA_snowfall_downscaled(range) =  forcing.DATA.ERA_snowfall_downscaled(range) .* snow_fraction(range);
                sublimation_factor = post_proc.PARA.scale_parameters{i,1}.sublimation_factor;
                sublimation_factor(isnan(sublimation_factor))=1;
                sublimation_factor_av = sublimation_factor_av + sublimation_factor;
                sublimation_factor = repmat(max(0, sublimation_factor), 1, size(timestamp,2), size(timestamp,3));
                forcing.DATA.ERA_sublimation_downscaled(range) = forcing.DATA.ERA_sublimation_downscaled(range) .* sublimation_factor(range);
            end

            % a=squeeze(forcing.DATA.ERA_sublimation_downscaled(2,:,:));
            % plot(a(:));
            
            snow_fraction_av = snow_fraction_av ./ count;
            snow_fraction_av = repmat(snow_fraction_av, 1, size(timestamp,2), size(timestamp,3));
            sublimation_factor_av = sublimation_factor_av ./ count;
            sublimation_factor_av = repmat(max(0, sublimation_factor_av), 1, size(timestamp,2), size(timestamp,3));
            forcing.DATA.ERA_snowfall_downscaled(~total_range) =  forcing.DATA.ERA_snowfall_downscaled(~total_range) .* snow_fraction_av(~total_range);
            forcing.DATA.ERA_sublimation_downscaled(~total_range) = forcing.DATA.ERA_sublimation_downscaled(~total_range) .* sublimation_factor_av(~total_range);
            
            % a=squeeze(forcing.DATA.ERA_sublimation_downscaled(2,:,:));
            % plot(a(:));

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