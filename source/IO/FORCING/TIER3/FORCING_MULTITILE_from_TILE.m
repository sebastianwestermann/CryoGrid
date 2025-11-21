%========================================================================
% CryoGrid FORCING class FORCING_MULTITILE_from_TILE

% Authors:
% S. Westermann, January 2024
%
%========================================================================

classdef FORCING_MULTITILE_from_TILE < FORCING_base
    
    properties

    end
    
    methods
        
        function forcing = provide_PARA(forcing)         
            forcing.PARA.forcing_class = [];
            forcing.PARA.forcing_class_index = [];
        end
        
               
        function forcing = finalize_init(forcing, tile)

            vars = {'latitude'; 'longitude'; 'altitude'; 'slope_angle'; 'aspect'; 'skyview_factor'; 'horizon_bins'; 'horizon_angles'};
            parameters = [];
            for i=1:size(vars,1)
                parameters = [parameters tile.PARA.(vars{i,1})];
            end
            % for i=1:size(forcing.PARA.variables,1)
            %     parameters = [parameters tile.(forcing.PARA.variables{i,1})];
            % end
            [~,ia, ic] = unique(parameters,'rows','legacy'); %finds the unique variables
            first_round = 1;
            tile2.RUN_INFO = tile.RUN_INFO;
            for i=1:size(ia,1)
                for j=1:size(vars,1)
                    tile2.PARA.(vars{j,1}) = tile.PARA.(vars{j,1})(i,:);
                end
                forcing_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(forcing.PARA.forcing_class){forcing.PARA.forcing_class_index,1});
                % for j=1:size(forcing.PARA.variables,1)
                %     forcing_class.PARA.(forcing.PARA.variables{j,1}) = tile.PARA.(forcing.PARA.variables{j,1})(i,:);
                % end
                forcing_class = finalize_init(forcing_class, tile2);
                fvars = fieldnames(forcing_class.DATA);
                for ii=1:size(fvars,1)
                    if size(forcing_class.DATA.(fvars{ii,1}),2) == 1 %single column
                        if first_round == 1
                            if ~strcmp(fvars{ii,1}, 'timeForcing') && ~strcmp(fvars{ii,1}, 'full_res_data')
                                forcing.DATA.(fvars{ii,1}) = repmat(forcing_class.DATA.(fvars{ii,1}).*0, 1, size(tile.latitude,2));
                                pos = find(ic == i);
                                for jj=1:length(pos)
                                    forcing.DATA.(fvars{ii,1})(:,pos(jj)) = forcing_class.DATA.(fvars{ii,1});
                                end
                            end
                            first_round = 0;
                        else 
                            if ~strcmp(fvars{ii,1}, 'timeForcing') && ~strcmp(fvars{ii,1}, 'full_res_data')
                                pos = find(ic == i);
                                for jj=1:length(pos)
                                    forcing.DATA.(fvars{ii,1})(:,pos(jj)) = forcing_class.DATA.(fvars{ii,1});
                                end
                            end
                        end
                    else %multi-column
                        if first_round == 1
                            if ~strcmp(fvars{ii,1}, 'timeForcing') && ~strcmp(fvars{ii,1}, 'full_res_data')
                                forcing.DATA.(fvars{ii,1}) = repmat(forcing_class.DATA.(fvars{ii,1}).*0, 1, size(tile.PARA.latitude,2));
                                pos = find(ic == i);
                                for jj=1:length(pos)
                                    forcing.DATA.(fvars{ii,1})(:,pos(jj)) = forcing_class.DATA.(fvars{ii,1})(:,jj);
                                end
                            end
                            first_round = 0;
                        else
                            if ~strcmp(fvars{ii,1}, 'timeForcing') && ~strcmp(fvars{ii,1}, 'full_res_data')
                                pos = find(ic == i);
                                for jj=1:length(pos)
                                    forcing.DATA.(fvars{ii,1})(:,pos(jj)) = forcing_class.DATA.(fvars{ii,1})(:,pos(jj));
                                end
                            end
                        end
                    end
                end
            end
            forcing.DATA.timeForcing = forcing_class.DATA.timeForcing;
            forcing.STATVAR.timestep = (forcing.DATA.timeForcing(end,1) - forcing.DATA.timeForcing(1,1)) ./ (size(forcing.DATA.timeForcing,1)-1);
            forcing.TEMP = forcing_class.TEMP;
            %assign PARAs to parent class
            vars = fieldnames(forcing_class.PARA);
            for i=1:size(vars,1)
                forcing.PARA.(vars{i,1}) = forcing_class.PARA.(vars{i,1});
            end
            if size(forcing.PARA.heatFlux_lb,2) == 1
                tile.PARA.geothermal = repmat(forcing.PARA.heatFlux_lb, 1, size(tile.PARA.latitude,2));
            else
                tile.PARA.geothermal = forcing.PARA.heatFlux_lb;
            end %this is an inconsitency between the two models, should be resolved at some point!

        end




    end
end