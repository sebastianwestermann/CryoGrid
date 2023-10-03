classdef ADJUST_STRATIGRAPHY < BASE
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    
    methods
        function grid = make_GRID_fixed(ground, tile)
            grid.STATVAR.GRID = ground.STATVAR.top_depth_rel2groundSurface + cumsum([0; ground.STATVAR.layerThick]);
            grid.STATVAR.MIDPOINTS = (grid.STATVAR.GRID(2:end,1) + grid.STATVAR.GRID(1:end-1,1))./2;
        end
        
        function grid = make_GRID_variable(ground, tile)
            
        end
        
        function ground = adjust_stratigraphy_waterIce_fixed(ground, tile)
            %for classes with fixed water+ice contents
           save_grid = tile.GRID;
           tile.GRID = make_GRID_fixed(ground, tile);
            for i=1:size(tile.PARA.strat_statvar_class,1)
                strat_statvar_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.strat_statvar_class{i,1}){tile.PARA.strat_statvar_class_index(i,1),1});
                strat_statvar_class = finalize_init(strat_statvar_class, tile);
            end
            
            ground_class_name = str2func(class(ground));
            ground2 = ground_class_name();
            variables = fieldnames(tile.GRID.STATVAR);
            for i = 1:size(variables,1)
                if ~strcmp(variables{i,1}, 'MIDPOINTS') && ~strcmp(variables{i,1}, 'GRID')
                    ground2.STATVAR.(variables{i,1}) = tile.GRID.STATVAR.(variables{i,1});
                end
            end
            ground2.STATVAR.layerThick = ground.STATVAR.layerThick;
            ground2 = convert_units(ground2, tile);
            tile.GRID = save_grid;
            
            variables = fieldnames(ground2.STATVAR);
            for i=1:size(variables,1)
                if ~strcmp(variables{i,1}, 'layerThick') && ~strcmp(variables{i,1}, 'T') && ~strcmp(variables{i,1}, 'area') 
                    ground.STATVAR.(variables{i,1}) = ground2.STATVAR.(variables{i,1});
                end
            end
            
            ground = finalize_init2(ground, tile);
            ground = compute_diagnostic(ground, tile);
        end
        
        function ground = adjust_stratigraphy_saturation_fixed(ground, tile)
            %for classes with variable water+ice contents, adjusts the
            %saturation to be as before
           save_grid = tile.GRID;
           tile.GRID = make_GRID_fixed(ground, tile);
            for i=1:size(tile.PARA.strat_statvar_class,1)
                strat_statvar_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.strat_statvar_class{i,1}){tile.PARA.strat_statvar_class_index(i,1),1});
                strat_statvar_class = finalize_init(strat_statvar_class, tile);
            end
            
            ground_class_name = str2func(class(ground));
            ground2 = ground_class_name();
            variables = fieldnames(tile.GRID.STATVAR);
            for i = 1:size(variables,1)
                if ~strcmp(variables{i,1}, 'MIDPOINTS') && ~strcmp(variables{i,1}, 'GRID')
                    ground2.STATVAR.(variables{i,1}) = tile.GRID.STATVAR.(variables{i,1});
                end
            end
            ground2.STATVAR.layerThick = ground.STATVAR.layerThick;
            ground2 = convert_units(ground2, tile);
            tile.GRID = save_grid;
            
            saturation = ground.STATVAR.waterIce ./ (ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.mineral - ground.STATVAR.organic);
            ground.STATVAR.waterIce = saturation .* (ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.mineral - ground.STATVAR.organic);
            
            variables = fieldnames(ground2.STATVAR);
            for i=1:size(variables,1)
                if ~strcmp(variables{i,1}, 'layerThick') && ~strcmp(variables{i,1}, 'T') && ~strcmp(variables{i,1}, 'area')
                    ground.STATVAR.(variables{i,1}) = ground2.STATVAR.(variables{i,1});
                end
            end
            
            ground = finalize_init2(ground, tile);
            ground = compute_diagnostic(ground, tile);
            
        end
        
        function ground = adjust_stratigraphy_Xice(ground, tile)
            
        end
        
    end
end

