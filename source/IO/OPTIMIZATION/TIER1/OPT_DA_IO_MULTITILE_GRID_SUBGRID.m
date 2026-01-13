%IO functions for DA for MULTITILEs


classdef DA_IO_MULTITILE_GRID_SUBGRID < DA_IO_MULTITILE 

    
    properties

    end
    
    methods


        function res = get_size_of_modeled_obs(da_IO, da, tile)
            res = tile.ENSEMBLE.PARA.grid_ensemble_size;
        end

    end

end

