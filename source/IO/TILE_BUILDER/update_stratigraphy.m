%========================================================================
% CryoGrid TLE_BUILDER class update_stratigraphy

% S. Westermann, Jan 2025
%========================================================================


classdef update_stratigraphy
    
    properties
        TILE
    end
    
    methods
        function build_tile(builder)
            builder.TILE = build_tile_update_stratigraphy(builder.TILE);            
        end
        
    end
end

