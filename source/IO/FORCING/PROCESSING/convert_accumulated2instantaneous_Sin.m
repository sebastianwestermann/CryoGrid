%========================================================================
% CryoGrid FORCING processing class
%
%
% Authors:
% S. Westermann, December 2022
%
%========================================================================

classdef convert_accumulated2instantaneous_Sin < terrain_correct_radiation
    

    
    methods
        function proc = provide_PARA(proc)

        end
        
        
        function proc = provide_CONST(proc)

        end
        
        
        function proc = provide_STATVAR(proc)
            
        end
        
        
        function proc = finalize_init(proc, tile)

        end
        
        
                %correct for the fact that ERA5 stores radiation as accumulated for
        %the 1h BEFORE the timestamp -> solution: caluclate average mu0 and
        %mu0 at the end of the interval = actual timestamp and scale with
        %this: CLEAN SOLUTION -> download hourly data and do this for the
        %timestep before and after, average between the two. Interpolate
        %for Lin, this has exactly the same issue, although it does not
        %interfere with the terrain scaling here!

        
        function forcing = process(proc, forcing, tile)
                        
            number_of_increments = 10;
            averaging_interval = 1; %in hours, 1 hour for ERA5
            timeForcing = repmat(forcing.DATA.timeForcing,1,number_of_increments);
            offset = linspace(-1,0,number_of_increments) .* averaging_interval./24;
            offset = repmat(offset, size(timeForcing,1),1);
            timeForcing = timeForcing + offset; 
            [azimuth,sunElevation] = solargeom(proc, timeForcing, forcing.SPATIAL.STATVAR.latitude, forcing.SPATIAL.STATVAR.longitude);
            sunElevation = 90-rad2deg(sunElevation);
            mu0=max(sind(sunElevation),0); % Trunacte negative values.
            scaling_factor = mu0(:, end) ./ mean(mu0, 2); 
            scaling_factor(isnan(scaling_factor)) = 0;
            forcing.DATA.Sin = forcing.DATA.Sin .* scaling_factor;
        
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

