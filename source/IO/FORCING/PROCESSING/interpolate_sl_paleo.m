%========================================================================
% CryoGrid FORCING processing class
%
%
% Authors:
% S. Westermann, December 2022
%
%========================================================================

classdef interpolate_sl_paleo < matlab.mixin.Copyable 
    
    properties
        CONST
        PARA
        STATVAR
        TEMP
    end
    
    methods
        function proc = provide_PARA(proc)

        end
        
        
        function proc = provide_CONST(proc)

        end
        
        
        function proc = provide_STATVAR(proc)
            
        end
        
        
        function proc = finalize_init(proc, tile)

        end
        
        
        function forcing = process(proc, forcing, tile)
            
            dist_lat = abs(forcing.SPATIAL.STATVAR.latitude - forcing.DATA.latitude);
            dist_lon=abs(forcing.SPATIAL.STATVAR.longitude - forcing.DATA.longitude);
            [dist_lat, ind_lat] = sort(dist_lat);
            [dist_lon, ind_lon] = sort(dist_lon);
            dist_lat=dist_lat(1:2);
            dist_lon=dist_lon(1:2);
            ind_lat = ind_lat(1:2);
            ind_lon = ind_lon(1:2);
            weights_lat = 1 - dist_lat./sum(dist_lat);
            weights_lon = 1 - dist_lon./sum(dist_lon);
            
            data = forcing.DATA;
            forcing.DATA=[];

            variables = {'precip'; 'Tair'; 'Lin'; 'Sin'; 'S_TOA'; 'RH'; 'wind'};
            
            for i=1:size(variables,1)
                if isfield(data, variables{i,1})
                    data.(variables{i,1}) = weights_lon(1) .* (weights_lat(1) .* squeeze(data.(variables{i,1})(ind_lon(1),ind_lat(1),:,:)) + weights_lat(2) .* squeeze(data.(variables{i,1})(ind_lon(1), ind_lat(2),:,:))) + ...
                        weights_lon(2) .* (weights_lat(1) .* squeeze(data.(variables{i,1})(ind_lon(2), ind_lat(1), :,:)) + weights_lat(2) .* squeeze(data.(variables{i,1})(ind_lon(2), ind_lat(2),:,:)));

                    forcing.DATA.(variables{i,1}) = [];
                    forcing.DATA.timeForcing=[];
                    for y=year(forcing.PARA.start_time):year(forcing.PARA.end_time)
                        for m=1:12
                            forcing.DATA.timeForcing = [forcing.DATA.timeForcing; datenum(y,m,15)];
                            index = find(y-data.year>=0,1, 'last');
                            fraction = (data.year(index+1,1)-y)./(data.year(index+1,1) - data.year(index,1));
                            forcing.DATA.(variables{i,1}) = [forcing.DATA.(variables{i,1}); fraction .* squeeze(data.(variables{i,1})(m,index)) + (1-fraction) .* squeeze(data.(variables{i,1})(m,index+1))] ;
                        end
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

