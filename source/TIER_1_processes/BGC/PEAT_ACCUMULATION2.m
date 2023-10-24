%========================================================================
% CryoGrid TIER1 library class PEAT_ACCUMULATION, containing functions related to peat accumulation in BGC_Frolking_peat
% S. Westermann, November 2021
%========================================================================

classdef PEAT_ACCUMULATION2 < BASE
    
    
    methods

        
%         function peat = peat_accumulation_Frolking_newCell(peat)
%             
%             %new_peat = (peat.STATVAR.peat_moss + peat.STATVAR.peat_shrub + peat.STATVAR.peat_graminoid) ./ peat.PARA.number_of_growing_days; % total peat, in g/cm2!!!
%             
%             %assume constant NPP
%             %new_peat = peat.STATVAR.annual_NPP .* peat.TEMP.GPP_acc./70e3 .* peat.PARA.BGC_timestep;
%             %new_peat = peat.STATVAR.annual_NPP .* peat.TEMP.GPP_acc./100e3 .* peat.PARA.BGC_timestep;
%             new_peat = peat.STATVAR.annual_NPP .* peat.TEMP.GPP_acc./100e3;
%             
%             peat.STATVAR.total_peat = [sum(new_peat,2); peat.STATVAR.total_peat];
%             peat.STATVAR.total_peat_PFT = [new_peat; peat.STATVAR.total_peat_PFT];
%             peat.STATVAR.totalpeatC_originalMass = [new_peat; peat.STATVAR.totalpeatC_originalMass]; %;  peat.STATVAR.totalpeatC_originalMass_old];
%             
%             peat.STATVAR.layerThick = (peat.STATVAR.total_peat./peat.PARA.bulkDensity);
%             
%             peat.TEMP.d_layerThick = [sum(new_peat,2) ./ peat.PARA.bulkDensity; peat.TEMP.d_layerThick];
%             peat.TEMP.d_organic = [ sum(new_peat,2) ./ peat.CONST.organicDensity ; peat.TEMP.d_organic];
%             
%         end
        
        function peat = peat_accumulation_Frolking(peat)
            
%              new_peat = peat.STATVAR.annual_NPP .* peat.TEMP.GPP_acc./100e3; 
             
             peat.TEMP.d_peat_accumulate_PFT = peat.STATVAR.annual_NPP .* peat.TEMP.GPP_acc./100e3; 
             
             %must come later in advance_prognostic()
%              peat.STATVAR.total_peat(1,1) = peat.STATVAR.total_peat(1,1) + sum(new_peat,2);
%              peat.STATVAR.total_peat_PFT(1,:) = peat.STATVAR.total_peat_PFT(1,:) + new_peat;
%              peat.STATVAR.total_peat_PFT_originalMass(1,:) = peat.STATVAR.total_peat_PFT_originalMass(1,:) + new_peat; %;  peat.STATVAR.totalpeatC_originalMass_old];
%              
%              peat.STATVAR.layerThick = (peat.STATVAR.total_peat./peat.PARA.bulkDensity); %CONVERSION MUST BE WRONG, SEE DECOMPOSITION.* peat.CONST.mtocm; % g/cm2 ./ cm3/g, then convert with m to cm
%              
%              peat.TEMP.d_layerThick(1,1) = peat.TEMP.d_layerThick(1,1) + sum(new_peat,2) ./ peat.PARA.bulkDensity;
%              peat.TEMP.d_organic(1,1) = peat.TEMP.d_organic(1,1) + sum(new_peat,2) ./ peat.CONST.organicDensity;
             
        end
        
        function peat = update_bulk_density(peat)
            
            %peat.STATVAR.bulkDensity = (peat.STATVAR.organic> 0).*peat.PARA.minbulkDensity+(peat.PARA.diffbulkDensity .*(1/(1+exp(-(40*(1-peat.STATVAR.massRemain(1))-34))))) + (peat.STATVAR.organic <= 0).* peat.PARA.mineral_bulkDensity;
            peat.STATVAR.bulkDensity = peat.PARA.min_bulkDensity + (peat.PARA.max_bulkDensity - peat.PARA.min_bulkDensity) .* (1./(1+exp(-(40*(1-peat.STATVAR.total_peat ./ peat.STATVAR.total_peat_originalMass)-34))));
            peat.STATVAR.bulkDensity(isnan(peat.STATVAR.bulkDensity)) = peat.PARA.min_bulkDensity;
             
        end
        
        function peat = define_NPP_variables_Frolking(peat)
            peat.PARA.NPP_max = [0.85; 0.85; 1.13; 0.56; 0.09; 0.19; 0.19; 0.56; 0.19; 0.19; ...
                0.19; 0.09]';
            
            peat.PARA.PFT_name = cell(12, 1);
            peat.PARA.PFT_name{1,1} = 'grass';
            peat.PARA.PFT_name{2,1} = 'Minerotrophicf orb';
            peat.PARA.PFT_name{3,1} = 'Minerotrophics edge';
            peat.PARA.PFT_name{4,1} = 'Minerotrophic shrub';
            peat.PARA.PFT_name{5,1} = 'Ombrotrophic forb';
            peat.PARA.PFT_name{6,1} = 'Ombrotrophic sedge';
            peat.PARA.PFT_name{7,1} = 'Ombrotrophic shrub';
            peat.PARA.PFT_name{8,1} = 'Brown moss';
            peat.PARA.PFT_name{9,1} = 'Hollow Sphagnum';
            peat.PARA.PFT_name{10,1} = 'Lawn Sphagnum';
            peat.PARA.PFT_name{11,1} = 'Hummock Sphagnum';
            peat.PARA.PFT_name{12,1} = 'feathermoss';
            
            peat.PARA.above_ground_fraction = [0.5; 0.5; 0.2; 0.5; 0.5; 0.2; 0.5; 1; 1; 1; 1; ...
                1]';
            
            peat.PARA.initial_decomposability = [0.20000000000001925; 0.40000000427289167; 0.30000000003861493; ...
                0.25099880287291643; 0.30000000003861493; 0.20000000000001925; ...
                0.20000000000001925; 0.1; 0.1; 0.066666666666666666; ...
                0.05191461747673335; 0.1]' ./ 365;
            
            initial_decomposability_k0 = [0.32; 0.88; 0.57; 0.44; 0.57; 0.32; 0.32; ...
                0.13; 0.13; 0.08; 0.06; 0.13]';
            
            peat.PARA.optimum_peat_depth = ...
                [0.01 1 1;
                0.3 1 1;
                0.1 2 2;
                1 2 2;
                4 2 19;
                4 2 19;
                4 2 19;
                0.1 1.5 1.5;
                2 1 19;
                2 1 19;
                2 1 19;
                4 6 19]';
            
            peat.PARA.optimum_water_table_depth = ...   % [0.4 0.4 0.4; 
                [0.4 0.4 2; 
                0.1 0.3 0.3; 
                0.1 0.4 0.4; 
                0.2 0.2 1; 
                0.2 0.2 0.2;
                0.2 0.3 0.3;
                0.3 0.3 1;
                0.01 0.2 0.05;
                0.01 0.2 0.05;
                0.1 0.3 0.4;
                0.2 0.1 0.5;
                0.4 0.4 0.6]';

        end
        
        function peat = get_annual_NPP_Frolking(peat)
            
            water_table_depth = mean(peat.STATVAR.water_table_depth,1);
            peat_depth = mean(peat.STATVAR.peat_depth,1);
            
            a = double(water_table_depth <= peat.PARA.optimum_water_table_depth(1,:)) .* ((water_table_depth - peat.PARA.optimum_water_table_depth(1,:)) ./ peat.PARA.optimum_water_table_depth(2,:)).^2 + ...
                double(water_table_depth > peat.PARA.optimum_water_table_depth(1,:)) .* ((water_table_depth - peat.PARA.optimum_water_table_depth(1,:)) ./ peat.PARA.optimum_water_table_depth(3,:)).^2;
            
            a= a + double(peat_depth <= peat.PARA.optimum_peat_depth(1,:)) .* ((peat_depth - peat.PARA.optimum_peat_depth(1,:)) ./ peat.PARA.optimum_peat_depth(2,:)).^2 + ...
                double(peat_depth > peat.PARA.optimum_peat_depth(1,:)) .* ((peat_depth - peat.PARA.optimum_peat_depth(1,:)) ./ peat.PARA.optimum_peat_depth(3,:)).^2;
            
            peat.STATVAR.annual_NPP = peat.PARA.NPP_max .* exp(-a);

        end
              
    end
end