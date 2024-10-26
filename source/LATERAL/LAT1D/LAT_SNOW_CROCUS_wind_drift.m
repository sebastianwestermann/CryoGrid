%========================================================================
% CryoGrid LATERAL_IA class LAT_SNOW_CROCUS_wind_drift 
% simulates lateral wind drift of snow from a 1D realization

% NOTE: works only with the SNOW classes SNOW_crocus_... and SNOW_crocus2_... 
% S. Westermann, Dec 2023
%========================================================================


classdef LAT_SNOW_CROCUS_wind_drift < BASE_LATERAL

    
    methods

        %----mandatory functions---------------
        %----initialization--------------------

        
        function lateral = provide_PARA(lateral)
            lateral.PARA.exposure = []; %exposure positive = ablation; exposure negative = deposition 
            lateral.PARA.n_drift = []; %5 .* swe_per_cell = 0.1;
            lateral.PARA.snow_holding_height = [];
            lateral.PARA.ia_time_increment = []; %0.05; %must be a multiple of the time increment of the main lateral class
        end
        
        function lateral = provide_CONST(lateral)
            lateral.CONST.day_sec = [];
        end
        
        
        function lateral = provide_STATVAR(lateral)
            
        end

        
        function lateral = finalize_init(lateral, tile)
            
        end
        
        %----time integration------------
        
        function lateral = push(lateral, tile)
            
            CURRENT = lateral.PARENT.TOP.NEXT;
            while ~(strcmp(class(CURRENT), 'Bottom'))
                if strcmp(class(CURRENT), 'SNOW_crocus_bucketW_seb') || strcmp(class(CURRENT), 'SNOW_crocus2_bucketW_seb')
                    CURRENT = lateral_pull_snow_wind_drift(CURRENT, lateral);
                    CURRENT = compute_diagnostic(CURRENT, tile);
                end
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
            
%             ground.PARA.class_category = 'LATERAL_IA';
%             
%             ground.PARA.options = [];
%             ground.PARA.STATVAR = [];
%             
%             ground.PARA.default_value.n_drift = {0.1};
%             ground.PARA.comment.n_drift ={'factor translating CROCUS driftability index to snow removed per unit time'};
%             
%             ground.PARA.default_value.weighting_factor_snow_dump = {0.5};
%             ground.PARA.comment.weighting_factor_snow_dump ={'strength of snow dump, i.e. more snow removed when higher, must be equal for all tiles!!'};
%             
%             ground.PARA.default_value.weighting_factor = {1};
%             ground.PARA.comment.weighting_factor ={'weighting factior own tile'};
%             
%             ground.PARA.default_value.snow_holding_height = {0.05};
%             ground.PARA.comment.snow_holding_height ={'snow not removed when below this height [m]'};
%             
%             ground.PARA.default_value.ia_time_increment = {0.25};
%             ground.PARA.comment.ia_time_increment ={'time step [days], must be multiple of of LATERAL class timestep'};
        end
    end
    
end


