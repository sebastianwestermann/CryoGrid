classdef FORCING_seb_mat < FORCING_base & READ_FORCING_mat
    
    properties
        
    end
    
    methods
        
        function forcing = finalize_init(forcing, tile)
            forcing = read_mat(forcing, tile)
        end
        
    end
    
end