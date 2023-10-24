%========================================================================
% CryoGrid TIER1 library class PEAT_DECOMPOSE, containing functions related to peat decompostion in BGC_Frolking_peat
% S. Westermann, November 2021
%========================================================================

classdef PEAT_DECOMPOSE2 < BASE
    
    
    methods
        
        function peat = temp_modifier(peat)
           peat.TEMP.tempModifier = double(peat.STATVAR.T >= 0) .* (peat.PARA.Q10.^(peat.STATVAR.T ./ 10));
        end
        
        function peat = water_modifier(peat)
            peat.TEMP.waterModifier = peat.STATVAR.vol_water .*0;
%             range = peat.STATVAR.vol_water > peat.PARA.fieldCapacity;
%             peat.TEMP.waterModifier(range,1) = 1.0-(1.0-0.025).* ((peat.STATVAR.vol_water(range,1) - peat.PARA.fieldCapacity)./(1-peat.PARA.fieldCapacity)).^3.0;
%             range = peat.STATVAR.vol_water >= 0.00 & peat.STATVAR.vol_water <= peat.PARA.fieldCapacity; %CHECK THIS!, should be 0.01 accoriding to Chaudhary 2017?
%             peat.TEMP.waterModifier(range,1) = 1.0-((peat.PARA.fieldCapacity - peat.STATVAR.vol_water(range,1))./peat.PARA.fieldCapacity).^5.0;%4.88);//IWRO 5 to 2//4.82
            
            saturated = peat.STATVAR.vol_water >0.99;
            depth_below_waterTable = cumsum(peat.STATVAR.layerThick .* double(saturated)) - peat.STATVAR.layerThick .* double(saturated) ./ 2;
            c1 = 2.31;
            vol_water_opt = 0.45;
            min_factor_anoxic = 0.001;
            max_factor_anoxic = 0.3;
            e_folding_depth = 0.3;
            
            peat.TEMP.waterModifier = double(~saturated) .* max(0, 1 - c1.*(peat.STATVAR.vol_water - vol_water_opt).^2) + double(saturated) .* ...
            (min_factor_anoxic + (max_factor_anoxic - min_factor_anoxic) .* exp(-depth_below_waterTable ./ e_folding_depth));
            peat.TEMP.waterModifier(isempty(peat.TEMP.waterModifier)) = 0;
        end
        
        
        function peat = peat_decompose_Frolking(peat)
            
            peat.TEMP.d_peat_decompose_PFT = repmat(peat.PARA.initial_decomposability, size(peat.STATVAR.total_peat_PFT,1),1).*(peat.STATVAR.total_peat_PFT./peat.STATVAR.total_peat_PFT_originalMass);  %0.05

            peat.TEMP.d_peat_decompose_PFT(isnan(peat.TEMP.d_peat_decompose_PFT)) = 0; %if totalpeatC_originalMass = 0 
            
            peat.TEMP.d_peat_decompose_PFT  = peat.TEMP.d_peat_decompose_PFT .* repmat(peat.TEMP.tempModifier,1,size(peat.PARA.initial_decomposability,2)) .* ...
                repmat(peat.TEMP.waterModifier,1,size(peat.PARA.initial_decomposability,2)); %.* peat.PARA.BGC_timestep; 
             %             peat.TEMP.d_organic = peat.TEMP.d_organic - sum(peat.TEMP.d_peat_decompose_PFT .* peat.STATVAR.total_peat_PFT, 2) ./ peat.CONST.organicDensity;

            
            %this needs to come later in advance_prognostic
%             peat.TEMP.d_layerThick = peat.TEMP.d_layerThick - sum(peat.TEMP.d_peat_decompose_PFT.*peat.STATVAR.total_peat_PFT,2) ./ peat.PARA.bulkDensity; % .* peat.CONST.mtocm; %add changes in bulk density
%             
%             
%             peat.STATVAR.total_peat_PFT = peat.STATVAR.total_peat_PFT - peat.TEMP.d_peat_decompose_PFT.*peat.STATVAR.total_peat_PFT;
%             peat.STATVAR.total_peat = sum(peat.STATVAR.total_peat_PFT,2);
%             peat.STATVAR.layerThick = (peat.STATVAR.total_peat./peat.PARA.bulkDensity); %Original must be wrong??? NPP in kg/m2, buld density in kg/m3  .* peat.CONST.mtocm; % g/cm2 ./ cm3/g, then convert with m to cm

        end
        

    end
end