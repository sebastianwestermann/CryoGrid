peat.PARA.isvascular = [1; 1; 0; 1; 1; 0; 1; 0; 0; 0; 0; 0];
peat.PARA.issedge =    [0; 0; 1; 0; 0; 1; 0; 0; 0; 0; 0; 0];
peat.PARA.above_ground_fraction_NPP = [0.5; 0.5; 0.2; 0.5; 0.5; 0.2; 0.5; 1; 1; 1; 1; 1];
peat.PARA.damping_factor_sedges = 5.365;


%depths = [0.005:0.01:2]';
%water_table_depth = -0.1;
% water_table_depth = mean(peat.STATVAR.water_table_depth,1);

depths = cumsum(peat.STATVAR.layerThick) - peat.STATVAR.layerThick./2;

NPP_distribution_vascular = (depths <= max(0.05, min(0.2, mean(peat.STATVAR.water_table_depth,1);)));
NPP_distribution_vascular = NPP_distribution_vascular ./ sum(NPP_distribution_vascular,1);

NPP_distribution_sedge = peat.PARA.damping_factor_sedges.* exp(-peat.PARA.damping_factor_sedges.*depths);
NPP_distribution_sedge(NPP_distribution_sedge<0.05) = 0;
NPP_distribution_sedge = NPP_distribution_sedge ./ sum(NPP_distribution_sedge,1);

NPP_distribution = repmat(depths .* 0, 1,12);
NPP_distribution(1,:) = peat.PARA.above_ground_fraction_NPP';
NPP_distribution(:,peat.PARA.isvascular'==1) = NPP_distribution(:,peat.PARA.isvascular'==1) + repmat(1-peat.PARA.above_ground_fraction_NPP(peat.PARA.isvascular==1)', size(NPP_distribution,1),1) ...
    .* repmat(NPP_distribution_vascular,1, sum(peat.PARA.isvascular));
NPP_distribution(:,peat.PARA.issedge'==1) = NPP_distribution(:,peat.PARA.issedge'==1) + repmat(1-peat.PARA.above_ground_fraction_NPP(peat.PARA.issedge==1)', size(NPP_distribution,1),1) ...
    .* repmat(NPP_distribution_sedge,1, sum(peat.PARA.issedge));






