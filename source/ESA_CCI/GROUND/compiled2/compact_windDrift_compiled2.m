function dD_dt = compact_windDrift_compiled2(T, ice_snow, layerThick_snow, snow_mat1, wind_compaction_timescale, g)


T = min(0,T);

eta_0 = 7.62237e6;
a_eta = 0.1;
b_eta = 0.023;
c_eta = 250;

rho_ice = 920;
rho_max = 350;

rho = ice_snow ./ max(1e-20, layerThick_snow) .* rho_ice;

stress = g .* rho_ice .* (cumsum(ice_snow) - ice_snow./2);

%SW added 4 which corresponds to old coarse-grained snow, compacts less!
eta = 4.* eta_0 .*  rho ./ c_eta .* exp(-a_eta .* T + b_eta .* rho);

dD_dt =  - stress ./max(1e-10,eta) ; %relative compaction

dD_dt = dD_dt - snow_mat1 .* rho_ice .* ice_snow ./ max(1e-8, rho).^2.* (rho_max - min(rho, rho_max)) ./ repmat(wind_compaction_timescale .*24.*3600, 4, 1) ./ max(1e-14, layerThick_snow);


