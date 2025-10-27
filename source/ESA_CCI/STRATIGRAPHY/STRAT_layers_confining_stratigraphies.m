

classdef STRAT_layers_confining_stratigraphies < matlab.mixin.Copyable
    
    properties
        
		PARA
        CONST
        STATVAR
    end
    
    methods
        

		function strat = provide_PARA(strat) 

			strat.PARA.ID = [];
            strat.PARA.name = [];
            
			strat.PARA.layers_organic_stratigraphy = [];
            strat.PARA.layers_mineral_stratigraphy = [];
            strat.PARA.depth2bedrock_organic = [];
            strat.PARA.depth2bedrock_mineral = [];
            
            strat.PARA.stratigraphy_fraction = []; %same as subgrid_index 
        end
        
        function strat = provide_CONST(strat)

        end
        
        function strat = provide_STATVAR(strat)

        end 
        
        
        function strat = finalize_init(strat, tile)
			
            
            variables = fieldnames(strat.PARA.layers_organic_stratigraphy);
            depth1 = strat.PARA.layers_organic_stratigraphy.depth;
            depth2bedrock_1 = depth1(end,:);
            depth1 = [depth1; Inf];
            depth2 = strat.PARA.layers_mineral_stratigraphy.depth;
            depth2bedrock_2 = depth2(end,:);
            depth2 = [depth2; Inf];
            depth2bedrock =  strat.PARA.stratigraphy_fraction.* strat.PARA.depth2bedrock_organic + strat.PARA.depth2bedrock_mineral .* (1-strat.PARA.stratigraphy_fraction);

            
            for i=1:size(variables,1)
                if ~strcmp(variables{i,1}, 'depth')
                    store_var1 = tile.GRID.STATVAR.GRID(1:end-1,1) .*0;
                    store_var2 = store_var1;
                    for j=1:size(depth1,1)-1
                        range = tile.GRID.STATVAR.MIDPOINTS > depth1(j,1) & tile.GRID.STATVAR.MIDPOINTS <= depth1(j+1,1) & tile.GRID.STATVAR.MIDPOINTS < depth2bedrock_1;
                        store_var1(range,1) = strat.PARA.layers_organic_stratigraphy.(variables{i,1})(j,1);
                    end
                  %  range = tile.GRID.STATVAR.MIDPOINTS >= depth2bedrock_1;
                   % store_var1(range,1) = strat.PARA.layers_organic_stratigraphy.(variables{i,1})(end,1);                    
                    for j=1:size(depth2,1)-1
                        range = tile.GRID.STATVAR.MIDPOINTS > depth2(j,1) & tile.GRID.STATVAR.MIDPOINTS <= depth2(j+1,1)  & tile.GRID.STATVAR.MIDPOINTS < depth2bedrock_2;
                        store_var2(range,1) = strat.PARA.layers_mineral_stratigraphy.(variables{i,1})(j,1);
                    end
                 %   range = tile.GRID.STATVAR.MIDPOINTS >= depth2bedrock_2;
                  %  store_var2(range,1) = strat.PARA.layers_mineral_stratigraphy.(variables{i,1})(end,1);   

                    strat.STATVAR.(variables{i,1}) = store_var1.* strat.PARA.stratigraphy_fraction + store_var2 .* (1-strat.PARA.stratigraphy_fraction);
                    
                    range = (tile.GRID.STATVAR.MIDPOINTS >= depth2bedrock);
                    strat.STATVAR.(variables{i,1})(range,1) = strat.PARA.layers_mineral_stratigraphy.(variables{i,1})(end,1);   
                end
            end
            
        end
        
            

        
    end
    
end