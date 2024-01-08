%========================================================================
% CryoGrid FORCING class FORCING_perturb_at_runtime

% Authors:
% S. Westermann, January 2024
%
%========================================================================

classdef FORCING_perturb_at_runtime < FORCING_base
    
    properties
        FORCING_CLASS
        PERTURB_CLASS
    end
    
    methods
        
        function forcing = provide_PARA(forcing)         

            forcing.PARA.forcing_class = [];
            forcing.PARA.forcing_class_index = [];
            forcing.PARA.perturb_forcing_class = [];   %perturbs forcing at runtime
            forcing.PARA.perturb_forcing_class_index = [];
        end
        
               
        function forcing = finalize_init(forcing, tile)
            
            forcing_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(forcing.PARA.forcing_class){forcing.PARA.forcing_class_index,1});
            forcing_class = finalize_init(forcing_class, tile);
            forcing.FORCING_CLASS = forcing_class;
                        
            %optional post-processing with dedicated classes
            if ~isempty(forcing.PARA.perturb_forcing_class) && sum(isnan(forcing.PARA.perturb_forcing_class_index)) == 0
                forcing.PERTURB_CLASS = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(forcing.PARA.perturb_forcing_class){forcing.PARA.perturb_forcing_class_index,1});
                forcing.PERTURB_CLASS = finalize_init(forcing.PERTURB_CLASS, tile);
                forcing.PERTURB_CLASS = preprocess(forcing.PERTURB_CLASS, forcing.FORCING_CLASS, tile);
            end
            forcing.PARA = forcing.FORCING_CLASS.PARA;
        end
        
        
        function forcing = interpolate_forcing(forcing, tile)
            forcing.FORCING_CLASS = interpolate_forcing(forcing.FORCING_CLASS, tile);
            forcing.TEMP = forcing.FORCING_CLASS.TEMP;
            forcing = perturb_forcing(forcing.PERTURB_CLASS, forcing, tile);
        end


    end
end