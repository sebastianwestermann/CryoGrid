%========================================================================
% CryoGrid FORCING processing class compute_time_average
%
% Brute-force time averaging of forcing data for SEB calculations, no ensemble generation 
% all variables stay the same

% Authors:
% S. Westermann, December 2022
%
%========================================================================

classdef compute_time_average < matlab.mixin.Copyable
    
    properties
        
    end
    
    methods
        function proc = provide_PARA(proc)
            
            proc.PARA.averaging_period = [];  %in days
           
        end
        
        
        function proc = provide_CONST(proc)
            proc.CONST.L_f = []; 
            proc.CONST.sigma = [];
            proc.CONST.day_sec = [];
            proc.CONST.Tmfw = [];
        end
        
        
        function proc = provide_STATVAR(proc)
            
        end
        
        
        function proc = finalize_init(proc, tile)
            
        end
        
        
        function forcing = process(proc, forcing, tile)
            
            data_full = forcing.DATA;
            forcing.DATA = [];
            forcing.DATA.snowfall = [];
            forcing.DATA.rainfall = [];
            forcing.DATA.Tair = [];
            forcing.DATA.wind = [];
            forcing.DATA.q = [];
            forcing.DATA.p = [];
            forcing.DATA.Lin = [];
            forcing.DATA.Sin = [];
            forcing.DATA.timeForcing = [];
            
            for i = data_full.timeForcing(1,1):proc.PARA.averaging_period:data_full.timeForcing(end,1)-proc.PARA.averaging_period
                range = find(data_full.timeForcing>=i & data_full.timeForcing < min(data_full.timeForcing(end,1), i + proc.PARA.averaging_period));
                forcing.DATA.timeForcing = [forcing.DATA.timeForcing; mean(data_full.timeForcing(range,1))];
                forcing.DATA.Tair = [forcing.DATA.Tair; mean(data_full.Tair(range,1))];
                forcing.DATA.wind = [forcing.DATA.wind; mean(data_full.wind(range,1))];
                forcing.DATA.q = [forcing.DATA.q; mean(data_full.q(range,1))];
                forcing.DATA.p = [forcing.DATA.p; mean(data_full.p(range,1))];
                forcing.DATA.Lin = [forcing.DATA.Lin; mean(data_full.Lin(range,1))];
                forcing.DATA.Sin = [forcing.DATA.Sin; mean(data_full.Sin(range,1))];   
                forcing.DATA.snowfall = [forcing.DATA.snowfall; mean(data_full.snowfall(range,1))];
                forcing.DATA.rainfall = [forcing.DATA.rainfall; mean(data_full.rainfall(range,1))];  
            end

            %overwrite target variables in TEMP in FORCING
            forcing.TEMP = [];
            forcing.TEMP.Lin = 0;
            forcing.TEMP.Sin = 0;
            forcing.TEMP.q = 0;
            forcing.TEMP.p = 0;
            forcing.TEMP.Tair = 0;
            forcing.TEMP.wind = 0;
            forcing.TEMP.rainfall = 0;
            forcing.TEMP.snowfall = 0;
        end
        
        
%                 %-------------param file generation-----
%         function proc = param_file_info(proc)
%             proc = provide_PARA(proc);
% 
%             proc.PARA.STATVAR = [];
%             proc.PARA.class_category = 'FORCING POST_PROCESSING';
%             proc.PARA.options = [];
%             
%             proc.PARA.eliminate_fraction = [];
%             proc.PARA.survive_fraction = [];
%                         
%             proc.PARA.default_value.window_size = {7};
%             proc.PARA.comment.window_size = {'window size in days within which precipitation is reallocated'};
%             
%             proc.PARA.default_value.eliminate_fraction = {0.5};
%             proc.PARA.comment.eliminate_fraction = {'fraction of smallest precipitation events (= timestamps with precipitation) that is reallocated to larger events'};
%             
%             proc.PARA.default_value.survive_fraction = {0.5};  
%             proc.PARA.comment.survive_fraction = {'fraction of largest precipitation events (= timestamps with precipitation) that the small events are reallocated to'};
%             
%         end
        
    end
    
end