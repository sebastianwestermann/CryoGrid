classdef READ_FORCING_NC < READ_FORCING_base
    
    properties
        
    end
    
    methods
        
       function forcing = read_NC_ERA5(forcing, tile) % Consider making seperate functions dependant on how data is stored
            % Read forcing data (ERA5 format) from NetCDF files, and
            % populate forcing.DATA
            
            variables = {'t2m'; 'd2m'; 'u10'; 'v10'; 'ssrd'; 'strd'; 'tp'; 'sp'; 'tisr'};
            for i=1:size(variables,1)
                temp.(variables{i,1}) = squeeze(ncread([forcing.PARA.forcing_path variables{i,1} '.nc'], variables{i,1}));
            end
            temp.time = ncread([forcing.PARA.forcing_path 't2m.nc'], 'time');
            
            temp.info = ncinfo([forcing.PARA.forcing_path 't2m.nc']);
            Index1 = find(contains({temp.info.Variables.Name},'time'));
            Index2 = find(contains({temp.info.Variables(Index1).Attributes.Name},'units'));
            reference = split(temp.info.Variables(Index1).Attributes(Index2).Value);
            reference = reference(3);
            reference_date = split(reference,"-");
            reference_date = str2double(reference_date);
            forcing.DATA.timeForcing = datenum(reference_date(1),reference_date(2),reference_date(3)) + double(temp.time)./24;
            
            forcing.DATA.Tair = temp.t2m - forcing.CONST.Tmfw;
            forcing.DATA.wind = sqrt(temp.u10.^2 + temp.u10.^2);
            forcing.DATA.Sin = temp.ssrd./3600;
            forcing.DATA.Sin = [0;0;forcing.DATA.Sin];
            forcing.DATA.Lin = temp.strd./3600;
            forcing.DATA.Lin = [0;0;forcing.DATA.Lin];
            forcing.DATA.p = temp.sp; %in [Pa]!900.*100 + forcing.DATA.Tair.*0;
            forcing.DATA.q = (double(forcing.DATA.Tair<0).*satPresIce(forcing, temp.d2m) + double(forcing.DATA.Tair>=0).*satPresWater(forcing, temp.d2m))./ forcing.DATA.p;
            forcing.DATA.precip = [0;0; temp.tp]; %append missing first two timesteps
            
            % RBZ Sep22: Might remove Sin TOA and calculate this using ERA5
            % parameterizations -> consistent solar geometry as in source
            forcing.DATA.S_TOA = temp.tisr ./ 3600;
            forcing.DATA.S_TOA = [0;0;forcing.DATA.S_TOA];
        end
         
    end
    
    
end
