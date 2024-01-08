%========================================================================
% CryoGrid FORCING processing class
%
%
% Authors:
% S. Westermann, December 2022
%
%========================================================================

classdef split_precip_snow_rain < matlab.mixin.Copyable 
    
    properties
        CONST
        PARA
        STATVAR
        TEMP
    end
    
    methods
        function proc = provide_PARA(proc)
            proc.PARA.all_rain_T = [];     % Temperature above which all precipitation is considered as rain
            proc.PARA.all_snow_T = [];     % Temperature below which all precipitation is considered as snow
        end
        
        
        function proc = provide_CONST(proc)

        end
        
        
        function proc = provide_STATVAR(proc)
            
        end
        
        
        function proc = finalize_init(proc, tile)

        end
        
        
        function forcing = process(proc, forcing, tile)

            % Split total precip into snow and rain depending on air temp.
            Tair = forcing.DATA.Tair;
            T_all_rain = proc.PARA.all_rain_T;
            T_all_snow = proc.PARA.all_snow_T;
            
            forcing.DATA.snowfall = forcing.DATA.precip .* (double(Tair <= T_all_snow)  + ...
                double(Tair > T_all_snow & Tair < T_all_rain) .* (1- (Tair - T_all_snow) ./ max(1e-12, (T_all_rain - T_all_snow))));
            forcing.DATA.rainfall = forcing.DATA.precip .* (double(Tair >= T_all_rain)  + ...
                double(Tair > T_all_snow & Tair < T_all_rain) .* (Tair - T_all_snow) ./ max(1e-12, (T_all_rain - T_all_snow)));
        
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

