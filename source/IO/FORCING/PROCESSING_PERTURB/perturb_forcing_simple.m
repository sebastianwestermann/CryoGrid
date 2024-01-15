classdef perturb_forcing_simple < matlab.mixin.Copyable

    
    properties
        PARA
        CONST
        STATVAR
    end
    
    methods
        function proc = provide_PARA(proc)
            proc.PARA.Tair_bias = [];
            proc.PARA.Sin_rel_error = [];
            proc.PARA.snow_fraction = [];
            proc.PARA.rain_fraction = [];
            proc.PARA.all_rain_T = [];
            proc.PARA.all_snow_T = [];
        end
        
        function proc = provide_CONST(proc)
            proc.CONST.sigma = [];
        end
        
        function proc = provide_STATVAR(proc)
            
        end
        
        function proc = finalize_init(proc, tile)

        end
        
        function forcing = process(proc, forcing, tile)
            %fill in with below for forcing.DATA correction!
        end
        
        %----------------------------------------------
        %proc
        function proc = preprocess(proc, forcing, tile)
            
        end
        
        function forcing = perturb_forcing(proc, forcing, tile)
            sky_emissivity = forcing.TEMP.Lin ./ (forcing.TEMP.Tair+273.15).^4 ./ proc.CONST.sigma;
            forcing.TEMP.Tair = forcing.TEMP.Tair + proc.PARA.Tair_bias;
            forcing.TEMP.Lin = sky_emissivity .* proc.CONST.sigma .* (forcing.TEMP.Tair+273.15).^4;
            forcing.TEMP.Sin = forcing.TEMP.Sin .* (1 + proc.PARA.Sin_rel_error);
            
            total_precip = forcing.TEMP.rainfall + forcing.TEMP.snowfall;
            forcing.TEMP.rainfall = total_precip .* (double(forcing.TEMP.Tair >= proc.PARA.all_rain_T) + ...
                double(forcing.TEMP.Tair > proc.PARA.all_snow_T & forcing.TEMP.Tair < proc.PARA.all_rain_T) .* ...
                (forcing.TEMP.Tair - proc.PARA.all_snow_T) ./ max(1e-12, (proc.PARA.all_rain_T - proc.PARA.all_snow_T)));
            forcing.TEMP.snowfall = total_precip - forcing.TEMP.rainfall;
            forcing.TEMP.snowfall = forcing.TEMP.snowfall .* proc.PARA.snow_fraction;
            forcing.TEMP.rainfall = forcing.TEMP.rainfall .* proc.PARA.rain_fraction;

        end

    end
end

