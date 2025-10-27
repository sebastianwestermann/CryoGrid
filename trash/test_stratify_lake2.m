clear all
layerThick = zeros(60,1)+0.1;
elevations = -cumsum([0;layerThick]);
elevations_midpoint = -(cumsum(layerThick)-layerThick./2);
density_water0 = 1e3;
density_air = 1.293; %[kg m−3];
% T=[4:13]';
% T(2) = 4;
% T(6:end) = [8:-1:4]';

T=linspace(0,4,size(layerThick,1))';


timestep=100; %[sec]

wind= 10;

g=9.81;
coeff_mix_conv = 0.125;
coeff_wind_stir = 0.23;
coeff_mix_turb = 0.51;
coef_wind_drag = 0.0013;

energy_avail_mix = 0; %comes from last timestep 

density_water = density_H2O(T); 

save_T=[];
save_rho=[];


for count=1:4000
    
%mixing from top

mixed_layer_bottom_cell = 1;
moment0 = 0;
moment1 = 0;
energy_mixing = 0;
density_mixed_layer = density_water(mixed_layer_bottom_cell,1);

while mixed_layer_bottom_cell<size(T,1) && density_mixed_layer >= density_water(mixed_layer_bottom_cell+1,1)
        
    density_x_lt = density_water(mixed_layer_bottom_cell,1) .* layerThick(mixed_layer_bottom_cell,1); 
    moment0 = moment0 + density_x_lt;
    moment1 = moment1 + density_x_lt .* elevations_midpoint(mixed_layer_bottom_cell,1);
    
    mixed_layer_bottom_cell = mixed_layer_bottom_cell +1;
    T_mixed_layer = mean(T(1:mixed_layer_bottom_cell,1));
    density_mixed_layer = density_H2O(T_mixed_layer);

end

%last layer
density_x_lt = density_water(mixed_layer_bottom_cell,1) .* layerThick(mixed_layer_bottom_cell,1);
moment0 = moment0 + density_x_lt;
moment1 = moment1 + density_x_lt .* elevations_midpoint(mixed_layer_bottom_cell,1);

%mix the layer
T_mixed_layer = mean(T(1:mixed_layer_bottom_cell,1));
T(1:mixed_layer_bottom_cell,1) = T_mixed_layer;
density_mixed_layer = density_H2O(T_mixed_layer);
density_water(1:mixed_layer_bottom_cell,1) = density_mixed_layer;

energy_convection = max(0, 0.5 .* g .* coeff_mix_conv ./density_mixed_layer .* (moment1 - moment0 .* 0.5.*(elevations(1,1)+elevations(mixed_layer_bottom_cell+1,1)))); 
u_star_sqared = 0.5 .* coef_wind_drag./ density_mixed_layer .* wind.^2;
energy_wind_stir = 0.5 .* coeff_wind_stir  .* u_star_sqared.^1.5 .* timestep;
energy_total_stir =  energy_convection + energy_wind_stir;
energy_avail_mix = energy_avail_mix + energy_total_stir;

q_sqr = max(1e-10, 2 .* energy_total_stir ./ coeff_mix_conv ./ timestep) .^(2/3);
if mixed_layer_bottom_cell<size(T,1)
    density_deviation = g.*(density_water(mixed_layer_bottom_cell+1,1) - density_mixed_layer)./(0.5.*(density_mixed_layer + density_water(mixed_layer_bottom_cell+1,1)));
    energy_uplift_cell_below_ml  = 0.5 .* (density_deviation .* sum(layerThick(1:mixed_layer_bottom_cell,1)) + coeff_mix_turb .* q_sqr) .* layerThick(mixed_layer_bottom_cell+1,1);
end


%here the condition needs to come
while mixed_layer_bottom_cell<size(T,1) && energy_avail_mix >= energy_uplift_cell_below_ml
    energy_avail_mix = energy_avail_mix - energy_uplift_cell_below_ml;
    mixed_layer_bottom_cell = mixed_layer_bottom_cell + 1;

    %next cell
    T_mixed_layer = mean(T(1:mixed_layer_bottom_cell,1));
    T(1:mixed_layer_bottom_cell,1) = T_mixed_layer;
    density_mixed_layer = density_H2O(T_mixed_layer);
    density_water(1:mixed_layer_bottom_cell,1) = density_mixed_layer;

    if mixed_layer_bottom_cell <size(T,1)
        density_deviation = g.*(density_water(mixed_layer_bottom_cell+1,1) - density_mixed_layer )./(0.5.*(density_mixed_layer + density_water(mixed_layer_bottom_cell+1,1)));
        energy_uplift_cell_below_ml = 0.5 .* (density_deviation .* sum(layerThick(1:mixed_layer_bottom_cell,1)) + coeff_mix_turb .* q_sqr) .* layerThick(mixed_layer_bottom_cell+1,1);
    end

end


%%
%mixing from the bottom


%bottom layer mixing
for mixed_layer_bottom_cell = size(T,1):-1:2
    moment0 = 0;
    moment1 = 0;
    energy_mixing = 0;
    density_mixed_layer = density_water(mixed_layer_bottom_cell,1);
    
    mixed_layer_top_cell = mixed_layer_bottom_cell;
    
    while mixed_layer_top_cell>1 && density_mixed_layer < density_water(mixed_layer_top_cell-1,1)
        
        %     density_x_lt = density_water(mixed_layer_bottom_cell,1) .* layerThick(mixed_layer_bottom_cell,1);
        %     moment0 = moment0 + density_x_lt;
        %     moment1 = moment1 + density_x_lt .* elevations_midpoint(mixed_layer_bottom_cell,1);
        
        mixed_layer_top_cell = mixed_layer_top_cell - 1;
        T_mixed_layer = mean(T(mixed_layer_top_cell:mixed_layer_bottom_cell,1));
        density_mixed_layer = density_H2O(T_mixed_layer);
        
    end
    
    % %last layer
    % density_x_lt = density_water(mixed_layer_bottom_cell,1) .* layerThick(mixed_layer_bottom_cell,1);
    % moment0 = moment0 + density_x_lt;
    % moment1 = moment1 + density_x_lt .* elevations_midpoint(mixed_layer_bottom_cell,1);
    
    %mix the layer
    T_mixed_layer = mean(T(mixed_layer_top_cell:mixed_layer_bottom_cell,1));
    T(mixed_layer_bottom_cell:mixed_layer_bottom_cell,1) = T_mixed_layer;
    density_mixed_layer = density_H2O(T_mixed_layer);
    density_water(mixed_layer_top_cell:mixed_layer_bottom_cell,1) = density_mixed_layer;
    
    mixed_layer_bottom_cell = mixed_layer_top_cell;
end

save_T = [save_T T];
save_rho=[save_rho density_water];

end


% w_star3 = g./density_mixed_layer ./timestep.* (moment1 - moment0 .* 0.5.*(elevations(1,1)+elevations(mixed_layer_bottom_cell+1,1)));  
% u_star2 = density_air./density_mixed_layer .* C_D .* wind.^2 ;
% 
% energy_mixing = energy_mixing + C_K .*(w_star3 + C_W .* u_star2.^1.5) .* timestep
% energy_uplift_cell_below_ml = real(density_deviation .* sum(layerThick(1:mixed_layer_bottom_cell,1)) + C_T .* (w_star3 + C_W .* u_star2.^1.5).^(2/3)) .* layerThick(mixed_layer_bottom_cell+1,1)
% energy_uplift_cell_below_ml = max(0,energy_uplift_cell_below_ml);


% %%
% %here the condition needs to come
% while mixed_layer_bottom_cell<size(T,1) && energy_mixing>energy_uplift_cell_below_ml
%     energy_mixing = energy_mixing - energy_uplift_cell_below_ml;
%     mixed_layer_bottom_cell = mixed_layer_bottom_cell +1;
% 
%     %next cell
%     T_mixed_layer = mean(T(1:mixed_layer_bottom_cell,1));
%     T(1:mixed_layer_bottom_cell,1) = T_mixed_layer;
%     density_mixed_layer = density_H2O(T_mixed_layer);
%     density_water(1:mixed_layer_bottom_cell,1) = density_mixed_layer;
%     
% %     moment0 = 0;
% %     moment1 = 0;
% %     for i=1:mixed_layer_bottom_cell
% %         density_x_lt = density_water(i,1) .* layerThick(i,1);
% %         moment0 = moment0 + density_x_lt;
% %         moment1 = moment1 + density_x_lt .* elevations_midpoint(i,1);
% %     end
% %     
% %     %recalculte w and ustar
% %     w_star3 = g./density_mixed_layer ./timestep.* (moment1 - moment0 .* 0.5.*(elevations(1,1)+elevations(mixed_layer_bottom_cell+1,1)));
% %     u_star2 = density_air./density_mixed_layer .* C_D .* wind.^2 ;
% %     energy_mixing = energy_mixing + C_K .*(w_star3 + C_W .* u_star2.^1.5) .* timestep;
%     density_deviation = g.*(density_mixed_layer - density_water(mixed_layer_bottom_cell+1,1))./(0.5.*(density_mixed_layer + density_water(mixed_layer_bottom_cell+1,1)));
%     energy_uplift_cell_below_ml = real(density_deviation .* sum(layerThick(1:mixed_layer_bottom_cell,1)) + C_T .* (w_star3 + C_W .* u_star2.^1.5).^(2/3)) .* layerThick(mixed_layer_bottom_cell+1,1);
%     energy_uplift_cell_below_ml = max(0,energy_uplift_cell_below_ml);
% end
% energy_mixing
% energy_uplift_cell_below_ml
% 




function out = density_H2O(T)

out = (999.83952 + 16.945176 .* T - 7.9870401e-3 .* T.^2 - 46.170461e-6 .* T.^3) ./ (1 + 16.897850e-3 .* T);

end