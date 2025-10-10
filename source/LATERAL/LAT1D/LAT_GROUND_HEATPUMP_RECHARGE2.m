%========================================================================
% CryoGrid LATERAL_IA class LAT_GROUND_HEATPUMP 

% S. Westermann, Aug 2025
%========================================================================


classdef LAT_GROUND_HEATPUMP_RECHARGE2 < BASE_LATERAL

    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        

        function lateral = provide_CONST(lateral)
            lateral.CONST.day_sec = [];
        end
        
        function lateral = provide_PARA(lateral)
            lateral.PARA.target_T_inside = [];
            lateral.PARA.heating_power_at_20_0 = []; 
            lateral.PARA.COP = ... %COP for different vaporization Ts, propane with condensation T of 40 degreeC, so that the house can be heated with 35 degreeC water
                [-40 2.3;
                -30 2.65;
                -20 3.12;
                -10 3.78;
                 0 4.79;
                10 6.48;
                20 9.88;
                30 20.84];
            lateral.PARA.COP2 = ... %COP for different vaporization Ts, propane with condensation T of 20 degreeC, so that the ground can be heated with 15 degreeC water
                [-40 2.96;
                -30 3.7;
                -20 4.68;
                -10 6.31;
                 0 9.59;
                10 19.49;
                20 100;
                30 100];
            %from https://tlk-energy.de/phasendiagramme/druck-enthalpie
        end
        
        function lateral = provide_STATVAR(lateral)

        end

        function lateral = finalize_init(lateral, tile)
            lateral.TEMP.R_house = lateral.PARA.target_T_inside ./ lateral.PARA.heating_power_at_20_0;        

            rng(1243)
            lateral.TEMP.electricity_price_table = exp(0.75*randn(367,1));
            lateral.TEMP.electricity_price_table = lateral.TEMP.electricity_price_table ./ mean(lateral.TEMP.electricity_price_table);
            lateral.TEMP.average_electricity_price = mean(lateral.TEMP.electricity_price_table,1);
        end

        %----time integration------
        
        %only push function needed
        function lateral = push(lateral, tile)
            
            CURRENT = lateral.PARENT.TOP.NEXT;
            lateral.STATVAR.Tair = tile.FORCING.TEMP.Tair;
            lateral.STATVAR.electricity_price = lateral.TEMP.electricity_price_table(max(1,round(tile.t-datenum(year(tile.t),1,1))),1);
            while ~(strcmp(class(CURRENT), 'Bottom'))
                CURRENT = lateral_push_groundHeatPump_recharge2(CURRENT, lateral);
                CURRENT = compute_diagnostic(CURRENT, tile);
                CURRENT = CURRENT.NEXT;
            end
            
        end
        
        
        function lateral = set_ACTIVE(lateral, i, t)
            lateral.PARENT.ACTIVE(i,1) = 1;
        end
        
        function lateral = get_derivatives(lateral, tile)
            
        end
        
        function lateral = pull(lateral, tile)
            
        end
        
        
        %-------------param file generation-----
        function ground = param_file_info(ground)
            ground = param_file_info@BASE_LATERAL(ground);
            
            ground.PARA.class_category = 'LATERAL_IA';
            
            ground.PARA.options = [];
            ground.PARA.STATVAR = [];
            
        end

        
    end
    
end


