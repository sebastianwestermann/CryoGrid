%========================================================================
% CryoGrid INTERACTION (IA) class for heat conduction and water fluxes between a SNOW class
% and a GROUND class with Richards equation water
% contains function for SNOW CHILD phase 
% S. Westermann, October 2020
%========================================================================

classdef IA_HEAT11_WATER11_RichardsEq_SNOW < IA_WATER & IA_HEAT 
    
    methods
        
        function get_boundary_condition_m(ia_heat_water, tile)
            get_boundary_condition_HEAT_m(ia_heat_water);
            get_boundary_condition_RichardsEq_SNOW_m(ia_heat_water);
        end
        
        %SNOW
        function get_IA_CHILD_boundary_condition_u(ia_heat_water, tile)
            get_boundary_condition_HEAT_IA_CHILD(ia_heat_water);
            get_boundary_condition_RichardsEq_SNOW_m(ia_heat_water);
        end

        function remove_excessWater_CHILD(ia_heat_water) %move excessWater from SNOW to water and Xwater from PARENT, 0 energy transfer since meltwater must be zero
            %should only be used for pressure classes
            
             space_left = max(0,ia_heat_water.NEXT.STATVAR.layerThick(1) .* ia_heat_water.NEXT.STATVAR.area(1) - ia_heat_water.NEXT.STATVAR.mineral(1) ...
                 - ia_heat_water.NEXT.STATVAR.organic(1) - ia_heat_water.NEXT.STATVAR.waterIce(1)); 
             water_in = min(space_left, ia_heat_water.PREVIOUS.STATVAR.excessWater);
             Xwater_in = max(0, ia_heat_water.PREVIOUS.STATVAR.excessWater - water_in);
            
             ia_heat_water.NEXT.STATVAR.waterIce(1) = ia_heat_water.NEXT.STATVAR.waterIce(1) + water_in;
             ia_heat_water.NEXT.STATVAR.water(1) = ia_heat_water.NEXT.STATVAR.water(1) + water_in;

             ia_heat_water.NEXT.STATVAR.SurfaceWaterIce = ia_heat_water.NEXT.STATVAR.SurfaceWaterIce + Xwater_in;
             ia_heat_water.NEXT.STATVAR.SurfaceWater = ia_heat_water.NEXT.STATVAR.SurfaceWater + Xwater_in;
             ia_heat_water.NEXT.STATVAR.SurfaceLayerThick = ia_heat_water.NEXT.STATVAR.SurfaceLayerThick + Xwater_in ./ ia_heat_water.NEXT.STATVAR.SurfaceArea;
             ia_heat_water.PREVIOUS.STATVAR.excessWater = 0;
             
             %nice to have variable
             ia_heat_water.NEXT.STATVAR.water_in_from_snow = ia_heat_water.NEXT.STATVAR.water_in_from_snow + water_in;
             ia_heat_water.NEXT.STATVAR.water_in_from_snow_surfaceWater = ia_heat_water.NEXT.STATVAR.water_in_from_snow_surfaceWater + Xwater_in;

        end
    end
end