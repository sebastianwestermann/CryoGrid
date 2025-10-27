%========================================================================
% CryoGrid TIER1 INTERACTION (IA) class for functions related to salt
% diffusion
% NOTE: at this pooint, only zero flux boundary conditions are implemented
% S. Westermann, October 2020
%========================================================================

classdef IA_SALT < IA_BASE
    
    methods
        
        function get_boundary_condition_ZEROFLUX_SALT_NEXT_m(ia_heat_water) %coupling between classes without (PREVIOUS) and with (NEXT) salt balance
            ia_heat_water.NEXT.TEMP.F_ub_salt = 0;
            ia_heat_water.NEXT.TEMP.d_salt(1) = ia_heat_water.NEXT.TEMP.d_salt(1) + 0;
        end
        
        function get_boundary_condition_ZEROFLUX_SALT_PREVIOUS_m(ia_heat_water) %coupling between classes with (PREVIOUS) and without (NEXT) salt balance
            ia_heat_water.PREVIOUS.TEMP.F_lb_salt = 0;
            ia_heat_water.PREVIOUS.TEMP.d_salt(end) = ia_heat_water.PREVIOUS.TEMP.d_salt(end) + 0;
        end
  
        function get_boundary_condition_SALT_m(ia_heat_water)
            stratigraphy1 = ia_heat_water.PREVIOUS;
            stratigraphy2 = ia_heat_water.NEXT;
            flux = (stratigraphy1.STATVAR.salt_c_brine(end) - stratigraphy2.STATVAR.salt_c_brine(1)) .* stratigraphy1.STATVAR.diffusivitySalt(end) .* stratigraphy2.STATVAR.diffusivitySalt(1) ./...
                (stratigraphy1.STATVAR.diffusivitySalt(end).* stratigraphy2.STATVAR.layerThick(1)./2 + stratigraphy2.STATVAR.diffusivitySalt(1).* stratigraphy1.STATVAR.layerThick(end)./2 );
            flux = flux .* min(stratigraphy1.STATVAR.area(end), stratigraphy2.STATVAR.area(1));
            
            stratigraphy1.TEMP.d_salt(end) = stratigraphy1.TEMP.d_salt(end) - flux;
            stratigraphy2.TEMP.d_salt(1) = stratigraphy2.TEMP.d_salt(1) + flux;
        end

        function get_boundary_condition_SALT_LAKE_m(ia_heat_water)
            stratigraphy1 = ia_heat_water.PREVIOUS;
            stratigraphy2 = ia_heat_water.NEXT;
            flux = (stratigraphy1.STATVAR.salt_c_brine(end) - stratigraphy2.STATVAR.salt_c_brine(1))  .* stratigraphy2.STATVAR.diffusivitySalt(1) ./ (stratigraphy2.STATVAR.layerThick(1) ./ 2);
            flux = flux .* min(stratigraphy1.STATVAR.area(end), stratigraphy2.STATVAR.area(1));

            stratigraphy1.TEMP.d_salt(end) = stratigraphy1.TEMP.d_salt(end) - flux;
            stratigraphy2.TEMP.d_salt(1) = stratigraphy2.TEMP.d_salt(1) + flux;
        end

        function get_boundary_condition_SALTbuoyancy_LAKE_m(ia_heat_water)
            stratigraphy1 = ia_heat_water.PREVIOUS;
            stratigraphy2 = ia_heat_water.NEXT;
            flux_diffusion = (stratigraphy1.STATVAR.salt_c_brine(end) - stratigraphy2.STATVAR.salt_c_brine(1,1))  .* stratigraphy2.STATVAR.diffusivitySalt(1) ./ stratigraphy2.STATVAR.layerThick(1);
            flux_buoyancy = double(stratigraphy1.STATVAR.density_water(end,1) > stratigraphy2.STATVAR.density_water(1,1)) .* (stratigraphy1.STATVAR.density_water(end,1)-stratigraphy2.STATVAR.density_water(1,1)) .* ...
                (stratigraphy1.STATVAR.salt_c_brine(end,1) - stratigraphy2.STATVAR.salt_c_brine(1,1)) .* stratigraphy2.STATVAR.diffusivity_buoyancy(1,1) ./ stratigraphy2.STATVAR.layerThick(1) ;
            flux = flux_diffusion + flux_buoyancy;
            flux = flux .* min(stratigraphy1.STATVAR.area(end), stratigraphy2.STATVAR.area(1));

            stratigraphy1.TEMP.d_salt(end) = stratigraphy1.TEMP.d_salt(end) - flux;
            stratigraphy2.TEMP.d_salt(1) = stratigraphy2.TEMP.d_salt(1) + flux;
        end        

    end
end


           
        
      