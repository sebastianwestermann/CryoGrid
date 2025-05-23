%========================================================================
% CryoGrid TIER1 library class, functions related to salt fluxes 
% NOTE: also contains functions related to the  freeze curve representation "fcSimple"
% S. Westermann, T. Ingeman-Nielsen, J. Scheer, October 2020
%========================================================================

classdef SALT < BASE
    
    methods
        %------boundary conditions-----
        
        function ground = get_boundary_condition_u_ZERO_SALT(ground)
            ground.TEMP.F_ub_salt = 0; %zero flux boundary condiition 
            ground.TEMP.d_salt(1) = ground.TEMP.d_salt(1) + ground.TEMP.F_ub_salt .* ground.STATVAR.area(1);
        end
        
        function ground = get_boundary_condition_l_ZERO_SALT(ground)
            ground.TEMP.F_lb_salt = 0; %zero flux boundary condition
            ground.TEMP.d_salt(end) = ground.TEMP.d_salt(end) + ground.TEMP.F_lb_salt .* ground.STATVAR.area(end);
        end
        

        %-----derivatives----------
        
        function ground = get_derivative_salt(ground)
            fluxes = (ground.STATVAR.salt_c_brine(1:end-1) - ground.STATVAR.salt_c_brine(2:end)) .* ground.STATVAR.diffusivitySalt(1:end-1) .* ground.STATVAR.diffusivitySalt(2:end) ./...
                (ground.STATVAR.diffusivitySalt(1:end-1).* ground.STATVAR.layerThick(2:end)./2 + ground.STATVAR.diffusivitySalt(2:end).* ground.STATVAR.layerThick(1:end-1)./2 );
            fluxes(isnan(fluxes)) = 0; % in case diffusivitySalt is 0 for both cells 
            %unit mol/m2, so flux through 1m2 unit cross section between cells!
            
            d_salt=ground.STATVAR.energy.*0;
            d_salt(1) =  - fluxes(1);
            d_salt(2:end-1) = fluxes(1:end-1) - fluxes(2:end);  
            d_salt(end) =  + fluxes(end);
            d_salt = d_salt.*ground.STATVAR.area;  %multiply by area, in [mol/sec] 
            
            ground.TEMP.d_salt = ground.TEMP.d_salt + d_salt;
        end

        function ground = get_derivative_salt_buoancy(ground)
            fluxes_diffusion = (ground.STATVAR.salt_c_brine(1:end-1) - ground.STATVAR.salt_c_brine(2:end)) .* ground.STATVAR.diffusivitySalt(1:end-1) .* ground.STATVAR.diffusivitySalt(2:end) ./...
                (ground.STATVAR.diffusivitySalt(1:end-1).* ground.STATVAR.layerThick(2:end)./2 +  ground.STATVAR.diffusivitySalt(2:end).* ground.STATVAR.layerThick(1:end-1)./2 );
            fluxes_diffusion(isnan(fluxes_diffusion)) = 0; % in case diffusivitySalt is 0 for both cells 
            %unit mol/m2, so flux through 1m2 unit cross section between cells!

            fluxes_buoyancy = double(ground.STATVAR.density_water(1:end-1,1) > ground.STATVAR.density_water(2:end,1)) .* (ground.STATVAR.density_water(1:end-1,1) - ground.STATVAR.density_water(2:end,1)) .* ...
                (ground.STATVAR.salt_c_brine(1:end-1) - ground.STATVAR.salt_c_brine(2:end)) .* ground.STATVAR.diffusivity_buoyancy(1:end-1) .* ground.STATVAR.diffusivity_buoyancy(2:end) ./...
                (ground.STATVAR.diffusivity_buoyancy(1:end-1).* ground.STATVAR.layerThick(2:end)./2 +  ground.STATVAR.diffusivity_buoyancy(2:end).* ground.STATVAR.layerThick(1:end-1)./2 );
            %Eq. 15 in https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2010JC006527

            d_salt=ground.STATVAR.energy.*0;
            d_salt(1) =  - fluxes_diffusion(1) - fluxes_buoyancy(1,1);
            d_salt(2:end-1) = fluxes_diffusion(1:end-1) - fluxes_diffusion(2:end) + fluxes_buoyancy(1:end-1) - fluxes_buoyancy(2:end) ;  
            d_salt(end) =  fluxes_diffusion(end) + fluxes_buoyancy(end);
            d_salt = d_salt.*ground.STATVAR.area;  %multiply by area, in [mol/sec] 
            
            ground.TEMP.d_salt = ground.TEMP.d_salt + d_salt;
        end
        
        
        %-----------timesteps----------
        
        function timestep = get_timestep_salt(ground) %not used at this stage, modify if necessary!
            timestep1 = 10 ./ (max(abs(ground.TEMP.d_salt) ./ ground.STATVAR.layerThick./ ground.STATVAR.area));
            timestep1(isnan(timestep1))  = ground.PARA.dt_max;
            timestep2 = nanmin(double(ground.TEMP.d_salt<0) .* (-ground.STATVAR.saltConc ./ground.TEMP.d_salt) + double(ground.TEMP.d_salt>=0) .* ground.PARA.dt_max);
            timestep2(isnan(timestep2))  = ground.PARA.dt_max;
            timestep = min(timestep1, timestep2); 
        end

       
        
        %----diagnostic step---------
        
        function ground = get_T_water_salt_fcSimple_Xice(ground)
            
            L_f = ground.CONST.L_f;
            c_w = ground.CONST.c_w;
            c_i = ground.CONST.c_i;
            c_o = ground.CONST.c_o;
            c_m = ground.CONST.c_m;
            R = ground.CONST.R;
            Tmfw = ground.CONST.Tmfw;
            
            deltaT = ground.STATVAR.deltaT;
            
            freeWaterIce = 0; %[m]  possible excess ice amount (i.e. free ice that melts at exactly 0 degree)
            
            waterIce = ground.STATVAR.waterIce ./ ground.STATVAR.area ./ ground.STATVAR.layerThick;
            mineral = ground.STATVAR.mineral ./ ground.STATVAR.area ./ ground.STATVAR.layerThick;
            organic = ground.STATVAR.organic ./ ground.STATVAR.area ./ ground.STATVAR.layerThick;
            energy = ground.STATVAR.energy ./ ground.STATVAR.area ./ ground.STATVAR.layerThick ;
            N = ground.STATVAR.saltConc./ ground.STATVAR.area ./ground.STATVAR.layerThick; %concentration per grid cell volume [mol/m3] CHECK, this should be mol ions per m3, not mol NaCl; so molarity of sea wter must be doubled
            
            A = 1 + c_i.*deltaT.*freeWaterIce./(waterIce.*L_f);            
            A1 = c_w.*waterIce + c_m.*mineral + c_o.*organic;            
            A2 = c_i.*waterIce + c_m.*mineral + c_o.*organic;            
            A3 = (c_w-c_i).*waterIce;            
            A4 = waterIce .* L_f;            
            B = -L_f./(R.* Tmfw.^2);
            
            %quadratic equation in Tm, the onset of freezing temperature
            %a*Tm^2 + b*Tm + c = 0
            
            %zero-th order terms
            c = - N .* A2 .* deltaT - N .* A4;
            
            %first order terms
            b = - N .* A3 + B .*(energy ./ A + L_f .*freeWaterIce ./ A + A2 .* deltaT + A4);
            
            %second-order terms
            a = B.* (-c_i.*freeWaterIce./A - (c_i .* deltaT .* freeWaterIce .*A1) ./ (A .*A4) - A2);
            
            %Tm_2 = (-b + sqrt(b.^2 - 4.*a.*c)) ./ (2.*a);
            Tm_1 = (-b - sqrt(b.^2 - 4.*a.*c)) ./ (2.*a); %this is the right branch of the two solutions
            
            thresh1 = 0;
            thresh2 = - L_f.*freeWaterIce;
            thresh3 = - L_f .* freeWaterIce + (c_w.*waterIce + c_m.*mineral + c_o.*organic + c_i .*freeWaterIce) .* (-R.* Tmfw.^2 ./L_f).* N;
            
            T = double(energy >= thresh1) .* energy ./(c_w .*waterIce + c_w.*freeWaterIce + c_m.*mineral + c_o.*organic) + ...
                double(energy < thresh2 & energy > thresh3) .* (energy + L_f.*freeWaterIce) ./  (c_w .*waterIce + c_i.*freeWaterIce + c_m.*mineral + c_o.*organic) + ...
                double(energy <= thresh3) .* (Tm_1 - deltaT + deltaT.* (energy./A - c_i.*freeWaterIce./A .*Tm_1 - c_i.*deltaT.*freeWaterIce./(A.*A4) .* A1.*Tm_1 + L_f.*freeWaterIce ./ A ...
                - (A2 .*(Tm_1-deltaT) - A4)) ./ (A3.*Tm_1 + A2.*deltaT + A4));
            
            ice = double(energy <= thresh3) .* waterIce .*(1-(energy./A - c_i.*freeWaterIce./A .*Tm_1 - c_i.*deltaT.*freeWaterIce./(A.*A4) .* A1.*Tm_1 + L_f.*freeWaterIce ./ A ...
                - (A2 .*(Tm_1-deltaT) - A4)) ./ (A3.*Tm_1 + A2.*deltaT + A4));
            water = waterIce - ice;
            
            freeIce = double(energy >= thresh2 & energy < thresh1) .* energy ./thresh2 .* freeWaterIce + double(energy < thresh2) .* freeWaterIce;
            freeWater = freeWaterIce - freeIce;
            
            ground.STATVAR.T = T;
            ground.STATVAR.ice = ice .* ground.STATVAR.layerThick .*  ground.STATVAR.area;
            ground.STATVAR.water = water .* ground.STATVAR.area .* ground.STATVAR.layerThick;
            
            salt_c_brine = N ./ water;
            salt_c_brine(isnan(salt_c_brine))=0;
            ground.STATVAR.salt_c_brine = salt_c_brine;
        end
        
        function ground = get_T_water_salt_freeW(ground)
            
            L_f = ground.CONST.L_f;
            c_w = ground.CONST.c_w;
            c_i = ground.CONST.c_i;
            c_o = ground.CONST.c_o;
            c_m = ground.CONST.c_m;
            R = ground.CONST.R;
            Tmfw = ground.CONST.Tmfw;
            
            waterIce = ground.STATVAR.waterIce ./ ground.STATVAR.area ./ ground.STATVAR.layerThick;
            mineral = ground.STATVAR.mineral ./ ground.STATVAR.area ./ ground.STATVAR.layerThick;
            organic = ground.STATVAR.organic ./ ground.STATVAR.area ./ ground.STATVAR.layerThick;
            energy = ground.STATVAR.energy ./ ground.STATVAR.area ./ ground.STATVAR.layerThick ;
            n = ground.STATVAR.saltConc./ ground.STATVAR.area ./ground.STATVAR.layerThick; %concentration per grid cell volume [mol/m3] CHECK, this should be mol ions per m3, not mol NaCl; so molarity of sea wter must be doubled
            
            A1 = c_w.*waterIce + c_m.*mineral + c_o.*organic;            
            A2 = c_i.*waterIce + c_m.*mineral + c_o.*organic;            
            A3 = c_w - c_i;            
            A4 = waterIce .* L_f;            
            B = R.* Tmfw.^2;
            
            %quadratic equation in Tm, the onset of freezing temperature
            %a*Tm^2 + b*Tm + c = 0
            
            %zero-th order terms
            c = - n .* B;
            
            %first order terms
            b = - n .* A3 .* B ./ L_f - energy - A4;
            
            %second-order terms
            a = A2;
            
            Tm_2 = (-b + sqrt(b.^2 - 4.*a.*c)) ./ (2.*a);
            Tm_1 = (-b - sqrt(b.^2 - 4.*a.*c)) ./ (2.*a); %this is the right branch of the two solutions
            
            E_f = - A1 .* n ./ waterIce .* R.* Tmfw.^2 ./ L_f; %freezing T
            
            T = double(energy >= E_f) .* energy ./(c_w .*waterIce + c_m.*mineral + c_o.*organic) + ...
                double(energy < E_f) .* Tm_1;
            
            water = double(energy >= E_f) .* waterIce + double(energy < E_f) .* - n ./ T .* R .* Tmfw.^2 ./ L_f;
            water(isnan(water)) = waterIce(isnan(water));
            ice = waterIce - water;
            
            ground.STATVAR.T = T;
            ground.STATVAR.ice = ice .* ground.STATVAR.layerThick .*  ground.STATVAR.area;
            ground.STATVAR.water = water .* ground.STATVAR.area .* ground.STATVAR.layerThick;
            
            salt_c_brine = n ./ water;
            salt_c_brine(isnan(salt_c_brine))=0;
            ground.STATVAR.salt_c_brine = salt_c_brine;
        end
        
        
        
        function ground = get_E_water_salt_FreezeDepress_Xice(ground) %used during initialization to calculate initial state for energy, water, salt concetrations
            
            L_f = ground.CONST.L_f;
            c_w = ground.CONST.c_w;
            c_i = ground.CONST.c_i;
            c_o = ground.CONST.c_o;
            c_m = ground.CONST.c_m;
            R = ground.CONST.R;
            Tmfw = ground.CONST.Tmfw;
            freeWaterIce = 0; %[m]  possible excess ice amount (i.e. free ice that melts at exactly 0 degree), set to 0 here          

            deltaT = ground.STATVAR.deltaT;

            mineral= ground.STATVAR.mineral ./ (ground.STATVAR.layerThick .* ground.STATVAR.area);
            organic = ground.STATVAR.organic ./ (ground.STATVAR.layerThick .* ground.STATVAR.area);
            waterIce = ground.STATVAR.waterIce ./ (ground.STATVAR.layerThick .* ground.STATVAR.area);

%             waterIce = ground.STATVAR.waterIce; % in volumetric fractions as provided by the initialization 
%             mineral = ground.STATVAR.mineral;
%             organic = ground.STATVAR.organic;
            area = ground.STATVAR.area;
            T = ground.STATVAR.T;
            
            N = ground.STATVAR.saltConc.*waterIce; 
            
            A = 1 + c_i.*deltaT.*freeWaterIce./(waterIce.*L_f);
            A1 = c_w.*waterIce + c_m.*mineral + c_o.*organic;
            A2 = c_i.*waterIce + c_m.*mineral + c_o.*organic;
            A3 = (c_w-c_i).*waterIce;
            A4 = waterIce .* L_f;
            B = -L_f./(R.* Tmfw.^2);
            
            %quadratic equation in Tm, the onset of freezing T
            %a*Tm^2 + b*Tm + c = 0
            
            %zero-th order terms
            c = N.*R.*Tmfw.^2./L_f.*deltaT;  
            
            %first order terms
            b=T+deltaT;
            
            %second-order terms
            a=-1;
            
            Tm_1 = (-b + sqrt(b.^2 - 4.*a.*c)) ./ (2.*a); %this is the right branch!!
            %Tm_2 = (-b - sqrt(b.^2 - 4.*a.*c)) ./ (2.*a);
            
            thresh1 = 0;
            thresh2 = (-R.* Tmfw.^2 ./L_f).* N;
            
            E1 = A2 .* (Tm_1-deltaT) - A4;
            E2 = A1 .*Tm_1;
            
            thirdTerm= (N./(-L_f./(R.*Tmfw.^2).*Tm_1) .*(E2-E1) + E1);
            thirdTerm(isnan(thirdTerm) | thirdTerm==Inf)=0;
            
            energy = double(T>= thresh1) .* (c_w .*waterIce + c_w.*freeWaterIce + c_m.*mineral + c_o.*organic) .*T + ...
                double(T < thresh1 & T >= thresh2) .* ((c_w .*waterIce + c_i.*freeWaterIce + c_m.*mineral + c_o.*organic) .*T - freeWaterIce .*L_f) + ...
                double(T < thresh2) .* thirdTerm;
            
            ice = double(T <= thresh2) .* waterIce .* (1-(energy - E1)./(E2-E1));
            water = waterIce - ice;
            
            ground.STATVAR.energy = energy .* ground.STATVAR.layerThick .* ground.STATVAR.area;
            ground.STATVAR.water = water .*  ground.STATVAR.layerThick .* ground.STATVAR.area;
            ground.STATVAR.ice = ice .* ground.STATVAR.layerThick .* ground.STATVAR.area;
            

            %ground.STATVAR.waterIce = ground.STATVAR.waterIce .* ground.STATVAR.layerThick .* ground.STATVAR.area;
            %ground.STATVAR.mineral = ground.STATVAR.mineral .* ground.STATVAR.layerThick .* ground.STATVAR.area;
            %ground.STATVAR.organic = ground.STATVAR.organic .* ground.STATVAR.layerThick .* ground.STATVAR.area;

            ground.STATVAR.air = ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.waterIce - ground.STATVAR.mineral - ground.STATVAR.organic;
            
            ground.STATVAR.saltConc = N .* ground.STATVAR.layerThick .* ground.STATVAR.area;  % [mol]
            
            salt_c_brine = N ./ water;  % [mol/m3]
            salt_c_brine(isnan(salt_c_brine))=0;
            ground.STATVAR.salt_c_brine = salt_c_brine;
        end
        
        
        function ground = get_E_water_salt_freeW(ground) %used during initialization to calculate initial state for energy, water, salt concetrations
            
            L_f = ground.CONST.L_f;
            c_w = ground.CONST.c_w;
            c_i = ground.CONST.c_i;
            c_o = ground.CONST.c_o;
            c_m = ground.CONST.c_m;
            R = ground.CONST.R;
            Tmfw = ground.CONST.Tmfw;
            
            mineral= ground.STATVAR.mineral ./ (ground.STATVAR.layerThick .* ground.STATVAR.area);
            organic = ground.STATVAR.organic ./ (ground.STATVAR.layerThick .* ground.STATVAR.area);
            waterIce = ground.STATVAR.waterIce ./ (ground.STATVAR.layerThick .* ground.STATVAR.area);

            T = ground.STATVAR.T;
            
            n = ground.STATVAR.saltConc.*waterIce; %salt concentration per grid cell
            
            A1 = c_w.*waterIce + c_m.*mineral + c_o.*organic;
            A2 = c_i.*waterIce + c_m.*mineral + c_o.*organic;
            
            T_f = - n ./ waterIce .* R.* Tmfw.^2 ./ L_f;
            E_f =  A1 .* T_f; %freezing T
            
            water = double(T >= T_f) .* waterIce + double(T < T_f) .* - n ./ T .* R.* Tmfw.^2 ./ L_f;
            ice = waterIce - water;
            energy = double(T >= T_f) .* A1 .* T + double(T < T_f) .*  (E_f + (T - T_f) .* A2 - L_f .* ice);

            ground.STATVAR.energy = energy .* ground.STATVAR.layerThick .* ground.STATVAR.area;
            ground.STATVAR.water = water .*  ground.STATVAR.layerThick .* ground.STATVAR.area;
            ground.STATVAR.ice = ice .* ground.STATVAR.layerThick .* ground.STATVAR.area;

            ground.STATVAR.air = ground.STATVAR.layerThick .* ground.STATVAR.area - ground.STATVAR.waterIce - ground.STATVAR.mineral - ground.STATVAR.organic;
            
            ground.STATVAR.saltConc = n .* ground.STATVAR.layerThick .* ground.STATVAR.area;  % [mol]
            
            salt_c_brine = n ./ water;  % [mol/m3]
            salt_c_brine(isnan(salt_c_brine))=0;
            ground.STATVAR.salt_c_brine = salt_c_brine;
        end
        
        
        
        %---diffusivity salt--------------
        
        function ground = diffusivity_salt(ground)
            
            water = ground.STATVAR.water./ground.STATVAR.layerThick ./ ground.STATVAR.area;
            
            D0 = ((6.06 + 9.60)/2  + max(ground.STATVAR.T, 0) .* (0.297  + 0.438)/2) .* 1e-10; %from Boudreau, B., 1997, Diagenetic Models and their implementation, Springer, Berlin.
            %average between values for Na+ and Cl-
            
            ground.STATVAR.diffusivitySalt = D0 .* water ./ground.PARA.tortuosity;
        end

        function ground = diffusivity_salt_buoyancy(ground)
            
            water = ground.STATVAR.water./ground.STATVAR.layerThick ./ ground.STATVAR.area;
            
            D0 = ((6.06 + 9.60)/2  + max(ground.STATVAR.T, 0) .* (0.297  + 0.438)/2) .* 1e-10; %from Boudreau, B., 1997, Diagenetic Models and their implementation, Springer, Berlin.
            %average between values for Na+ and Cl-
            ground.STATVAR.diffusivitySalt = D0 .* water ./ ground.PARA.tortuosity;

            mixing_length = 7e-5; %[m] no idea where this comes from in the first place, from here: https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2010JC006527
            viscosity_water = 1.79e-3; %[Pa sec] coarse estimate for 0 degreeC, make T-dpendent, but should not play a big role
            ground.STATVAR.diffusivity_buoyancy = ground.CONST.g ./ viscosity_water .* mixing_length .* 3e-8 .* water.^3 ./ ground.PARA.tortuosity_buoyancy;
        end
        
        function ground = diffusivity_salt_sea_ice(ground)
            
            water = ground.STATVAR.water./ground.STATVAR.layerThick ./ ground.STATVAR.area;
            
            D0 = ((6.06 + 9.60)/2  + max(ground.STATVAR.T, 0) .* (0.297  + 0.438)/2) .* 1e-10; %from Boudreau, B., 1997, Diagenetic Models and their implementation, Springer, Berlin.
            %average between values for Na+ and Cl-
            
            ground.STATVAR.diffusivitySalt = D0 .* water ;
        end

        function ground = diffusivity_salt_sea_ice_buoyancy(ground)
            
            water = ground.STATVAR.water./ground.STATVAR.layerThick ./ ground.STATVAR.area;
            
            D0 = ((6.06 + 9.60)/2  + max(ground.STATVAR.T, 0) .* (0.297  + 0.438)/2) .* 1e-10; %from Boudreau, B., 1997, Diagenetic Models and their implementation, Springer, Berlin.
            %average between values for Na+ and Cl-
            
            ground.STATVAR.diffusivitySalt = D0 .* water ;

            mixing_length = 7e-5; %[m] no idea where this comes from in the first place, from here: https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2010JC006527
            viscosity_water = 1.79e-3; %[Pa sec] coarse estimate for 0 degreeC, make T-dpendent, but should not play a big role
            
            ground.STATVAR.diffusivity_buoyancy = ground.CONST.g ./ viscosity_water .* mixing_length .* 3e-8 .* water.^3 ;

        end

        function ground = water_density_saline2(ground)

            T = max(0,ground.STATVAR.T);
            %required salt_conc in g salt/kg water 
%             salt_conc = salt_conc .* 1e-3 .* (23+35.5)/2; %salt_conc in the model in mol/m3 = 1e-3 mol/l; assuming NaCl means half of the atoms have atomic mass 23 and half atomic mall 35.5
            salt_conc = ground.STATVAR.salt_c_brine .* 1e-3 .* 35.5./1.1243; %salt_conc in the model in mol/m3 = 1e-3 mol/l; assuming an average ion composition of sea water, 35g/l correpsonds to 1.1243 mol/l (density of water set to 1000 for simplicity)
            
            c1=999.842594;
            c2=6.793952e-2;
            c3=9.095290e-3;
            c4=1.001685e-4;
            c5=1.120083e-6;
            c6=6.536332e-9;
            d1=8.24493e-1;
            d2=4.0899e-3;
            d3=7.6438e-5;
            d4=8.2467e-7;
            d5=5.3875e-9;
            d6=5.72466e-3;
            d7=1.0227e-4;
            d8=1.6546e-6;
            d9=4.8314e-4;
            
            t2 = T.^2;
            t3 = t2.*T;
            t4 = t3.*T;
            t5 = t4.*T;
            s2 = salt_conc.*salt_conc;
            s32 = salt_conc.^1.5;
            
            term1  =  c1;
            term2  =  c2 .* T;
            term3  = -c3 .* t2;
            term4  =  c4 .* t3;
            term5  = -c5 .* t4;
            term6  =  c6 .* t5;
            term7  =  d1;
            term8  = -d2 .* T;
            term9  =  d3 .* t2;
            term10  = -d4 .* t3;
            term11 =  d5 .* t4;
            term12 = -d6;
            term13 =  d7 .* T;
            term14 = -d8 .* t2;
            term15 =  d9;
            
            dpure  =  term6  + term5  + term4  + term2 + term3 + term1;
            csal1  = (term11 + term10  + term9  + term8 + term7) .* salt_conc;
            csal32 = (term14 + term13 + term12) .* s32;
            csal2  =  term15 .* s2;
            
            ground.STATVAR.density_water = min(1214, dpure + csal1 + csal32 + csal2); %maximum density is about 1214 kg/m3 (about 26 wt%), eutectic point reached thereafter
        end
        
    end
end

