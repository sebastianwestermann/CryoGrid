%========================================================================
% CryoGrid TIER1 library class for functions for lateral fluxes of salt
% contains push and pull functions used in lateral salt flux classes, e.g. LAT3D_SALT 
% S. Westermann, Nov 2025
%========================================================================


classdef SALT_FLUXES_LATERAL < BASE


    methods

        
        %read information from GROUND class and send it to the LATERAL class   
        function ground = lateral3D_pull_salt_simple(ground, lateral)         
            if isempty(lateral.PARENT.STATVAR.depths_salt)
                depths = ground.STATVAR.upperPos - cumsum([0; ground.STATVAR.layerThick]);
            else
                depths = ground.STATVAR.upperPos - cumsum(ground.STATVAR.layerThick);
            end
            lateral.PARENT.STATVAR.depths_salt = [lateral.PARENT.STATVAR.depths_salt; depths];
            lateral.PARENT.STATVAR.diffusivity_salt = [lateral.PARENT.STATVAR.diffusivity_salt; ground.STATVAR.diffusivitySalt];
            lateral.PARENT.STATVAR.salt_c_brine = [lateral.PARENT.STATVAR.salt_c_brine; ground.STATVAR.salt_c_brine];
        end

         % add lateral heat flux to STATVAR energy 
        function ground = lateral3D_push_salt_simple(ground, lateral)
            ground.STATVAR.saltConc = ground.STATVAR.saltConc + lateral.PARENT.STATVAR.salt_flux(1:size(ground.STATVAR.saltConc,1),1);
            lateral.PARENT.STATVAR.salt_flux(1:size(ground.STATVAR.saltConc,1),:) = [];
        end

    end
end

