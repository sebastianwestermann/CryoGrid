%========================================================================
% CryoGrid INTERACTION (IA) class for heat conduction between a LAKE and a
% normal GROUND class
% S. Westermann, October 2020
%========================================================================

classdef IA_HEAT11_SALTbuoyancy11_LAKE <  IA_HEAT & IA_SALT
    
    methods
        
        function get_boundary_condition_m(ia_heat, tile)
            get_boundary_condition_HEAT_LAKE_m(ia_heat);
            get_boundary_condition_SALTbuoyancy_LAKE_m(ia_heat);
        end
           
    end
end