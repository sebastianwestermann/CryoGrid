classdef UNIT_CONVERSION_standard < matlab.mixin.Copyable

    methods

        function ground = convert_normal(unit_converter, ground, tile)
            ground.PARA.heatFlux_lb = tile.FORCING.PARA.heatFlux_lb;
            ground.PARA.airT_height = tile.FORCING.PARA.airT_height;
            ground.STATVAR.area = tile.PARA.area + ground.STATVAR.layerThick .* 0;

            ground.STATVAR.waterIce = ground.STATVAR.waterIce .* ground.STATVAR.layerThick .* ground.STATVAR.area;
            ground.STATVAR.mineral = ground.STATVAR.mineral .* ground.STATVAR.layerThick .* ground.STATVAR.area;
            ground.STATVAR.organic = ground.STATVAR.organic .* ground.STATVAR.layerThick .* ground.STATVAR.area;

        end
        
        function ground = convert_Xice(unit_converter, ground, tile)
            ground.PARA.heatFlux_lb = tile.FORCING.PARA.heatFlux_lb;
            ground.PARA.airT_height = tile.FORCING.PARA.airT_height;
            ground.STATVAR.area = tile.PARA.area + ground.STATVAR.layerThick .* 0;
            ground.STATVAR.Xice(ground.STATVAR.T>0) = 0;

            ground.STATVAR.waterIce = ground.STATVAR.waterIce .* ground.STATVAR.layerThick .* ground.STATVAR.area ./ (1 + ground.STATVAR.Xice);
            ground.STATVAR.mineral = ground.STATVAR.mineral .* ground.STATVAR.layerThick .* ground.STATVAR.area ./ (1 + ground.STATVAR.Xice);
            ground.STATVAR.organic = ground.STATVAR.organic .* ground.STATVAR.layerThick .* ground.STATVAR.area ./ (1 + ground.STATVAR.Xice);
            ground.STATVAR.XwaterIce = ground.STATVAR.Xice./ (1 + ground.STATVAR.Xice) .* ground.STATVAR.layerThick .* ground.STATVAR.area;
        end
        
        function ground = convert_Xice_GLACIER(unit_converter, ground, tile)
            ground.PARA.heatFlux_lb = tile.FORCING.PARA.heatFlux_lb;
            ground.PARA.airT_height = tile.FORCING.PARA.airT_height;
            ground.STATVAR.area = tile.PARA.area + ground.STATVAR.layerThick .* 0;

            ground.STATVAR.layerThick_glacier = ground.STATVAR.glacier_fraction .* ground.STATVAR.layerThick;
            ground.STATVAR.layerThick_sediment = ground.STATVAR.layerThick - ground.STATVAR.layerThick_glacier;

            ground.STATVAR.waterIce = ground.STATVAR.ice_fraction .* ground.STATVAR.layerThick_glacier .* ground.STATVAR.area + ...
                (1 - ground.STATVAR.mineral - ground.STATVAR.organic ) .*  ground.STATVAR.layerThick_sediment .* ground.STATVAR.area;
            ground.STATVAR.mineral = ground.STATVAR.mineral .* ground.STATVAR.layerThick_sediment .* ground.STATVAR.area;
            ground.STATVAR.organic = ground.STATVAR.organic .* ground.STATVAR.layerThick_sediment .* ground.STATVAR.area;              
        end

        function ground = convert_normal_pressure(unit_converter, ground, tile)
            ground.PARA.heatFlux_lb = tile.FORCING.PARA.heatFlux_lb;
            ground.PARA.airT_height = tile.FORCING.PARA.airT_height;
            ground.STATVAR.area = tile.PARA.area + ground.STATVAR.layerThick .* 0;
                        
            ground.STATVAR.initial_layerThick = ground.STATVAR.layerThick;
            ground.STATVAR.initial_porosity = 1 - (ground.STATVAR.mineral + ground.STATVAR.organic);
            ground.STATVAR.initial_voidRatio = ground.STATVAR.initial_porosity ./ (1 - ground.STATVAR.initial_porosity);
            
            %ground.STATVAR.waterIce = ground.STATVAR.waterIce .* ground.STATVAR.layerThick .* ground.STATVAR.area; % not calculated as saturation is given instead of waterIce
            ground.STATVAR.mineral = ground.STATVAR.mineral .* ground.STATVAR.layerThick .* ground.STATVAR.area;
            ground.STATVAR.organic = ground.STATVAR.organic .* ground.STATVAR.layerThick .* ground.STATVAR.area;
            ground.STATVAR.waterIce = ground.STATVAR.saturation .* ground.STATVAR.initial_porosity .* ground.STATVAR.initial_layerThick .* ground.STATVAR.area;
        end


        function ground = convert_normal_ubT(unit_converter, ground, tile)
            ground.PARA.heatFlux_lb = tile.FORCING.PARA.heatFlux_lb;
            ground.STATVAR.area = tile.PARA.area + ground.STATVAR.layerThick .* 0;

            ground.STATVAR.waterIce = ground.STATVAR.waterIce .* ground.STATVAR.layerThick .* ground.STATVAR.area;
            ground.STATVAR.mineral = ground.STATVAR.mineral .* ground.STATVAR.layerThick .* ground.STATVAR.area;
            ground.STATVAR.organic = ground.STATVAR.organic .* ground.STATVAR.layerThick .* ground.STATVAR.area;

        end


        function ground = convert_normal_ubT_lbT(unit_converter, ground, tile)
            ground.STATVAR.area = tile.PARA.area + ground.STATVAR.layerThick .* 0;

            ground.STATVAR.waterIce = ground.STATVAR.waterIce .* ground.STATVAR.layerThick .* ground.STATVAR.area;
            ground.STATVAR.mineral = ground.STATVAR.mineral .* ground.STATVAR.layerThick .* ground.STATVAR.area;
            ground.STATVAR.organic = ground.STATVAR.organic .* ground.STATVAR.layerThick .* ground.STATVAR.area;

        end        
    end
end

