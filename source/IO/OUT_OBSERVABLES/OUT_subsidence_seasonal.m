classdef OUT_subsidence_seasonal < matlab.mixin.Copyable
    
    properties
        PARA
        CONST
        TEMP
        STATVAR
        TIMESTAMP
        OUTPUT_TIME
        SAVE_TIME
    end
    
    methods

        function out = provide_PARA(out)
            out.PARA.output_timestep = [];
            out.PARA.save_date = [];
            out.PARA.save_interval = [];
            out.PARA.timestamps = [];
            out.PARA.tag = [];
            out.PARA.tag2 = [];
            out.PARA.add_Xice_subsidence = 0;

            out.PARA.timestamps = [];
        end
        
        function out = provide_CONST(out)
            
        end
        
        function out = provide_STATVAR(out)
            
        end
        
        function out = finalize_init(out, tile)
            
            if ~isempty(out.PARA.tag) && sum(isnan(out.PARA.tag))>0
                out.PARA.tag = [];
            end
                
            out.STATVAR.subsidence = []; 
            out.TIMESTAMP = [];

            if ~isempty(out.PARA.timestamps) %DA mode, defined timestamps written by DA
                out.OUTPUT_TIME = out.PARA.timestamps(find(out.PARA.timestamps(:,1)-tile.PARA.start_time > 0, 1, 'first'), 1);
                out.SAVE_TIME = Inf;
            else
                out.OUTPUT_TIME = tile.PARA.start_time + out.PARA.output_timestep;
                if isempty(out.PARA.save_interval) || isnan(out.PARA.save_interval)
                    out.SAVE_TIME = tile.FORCING.PARA.end_time;
                else
                    out.SAVE_TIME = min(tile.FORCING.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(tile.PARA.start_time,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
                end
            end
            if isempty(out.OUTPUT_TIME)
                out.OUTPUT_TIME = Inf;
            end

            out.TEMP.initial_ground_surface_elevation = 0;
            out.TEMP.current_year = 1e20; % random unrealistic number, will be set to current year after first annual observation
            out.TEMP.initial_poreice = 0; %[m3]
            out.TEMP.initial_saturation =0;
            
            %out.TEMP.initial_pore_ice_heave = 0;
        end
        
        function out = store_OUT(out, tile)

            if tile.t >= out.OUTPUT_TIME

                out.TIMESTAMP = [out.TIMESTAMP tile.t];

                if ~isempty(out.PARA.timestamps)
                    out.OUTPUT_TIME = out.PARA.timestamps(find(out.PARA.timestamps(:,1)-tile.t>0, 1, 'first'), 1);
                else
                    out.OUTPUT_TIME = min(out.SAVE_TIME, out.OUTPUT_TIME + out.PARA.output_timestep);
                end
                if isempty(out.OUTPUT_TIME)
                    out.OUTPUT_TIME = Inf;
                end


                %% calc section
                result = 0;
                CURRENT = tile.TOP.NEXT;

                while ~is_ground_surface(CURRENT) && ~(strcmp(class(CURRENT), 'Bottom'))
                    CURRENT = CURRENT.NEXT;
                end

                if year(tile.t) == out.TEMP.current_year  % for all observations after the first (= 0) of this year
    
                    % subsidence from xice
                    ground_surface_elevation = sum(CURRENT.STATVAR.layerThick) - CURRENT.STATVAR.XwaterIce(1,1) ./ CURRENT.STATVAR.area(1,1);
                    subsidence_xice = ground_surface_elevation - out.TEMP.initial_ground_surface_elevation;
                    
                    % subsidence from pore ice melt
                    poreice = CURRENT.STATVAR.ice;
                    air = (CURRENT.STATVAR.layerThick.*CURRENT.STATVAR.area) - CURRENT.STATVAR.waterIce -CURRENT.STATVAR.mineral -CURRENT.STATVAR.organic - CURRENT.STATVAR.XwaterIce; % air is not updated but htere is a field, calculate new
                    air = max(0,air);
                    saturation = CURRENT.STATVAR.ice ./max(1e-100,CURRENT.STATVAR.ice+air);
                    saturation = max(0,min(1,saturation));
                    subsidence_poreice=0;
                    for cell=1:size(poreice,1)
                        if poreice(cell) > out.TEMP.initial_poreice(cell)
                            subsidence_poreice_cell = (poreice(cell) - out.TEMP.initial_poreice(cell)) * sensitivity_heave(saturation(cell)) * 0.08;
    
                        elseif poreice(cell) < out.TEMP.initial_poreice(cell)
                            subsidence_poreice_cell = (poreice(cell) - out.TEMP.initial_poreice(cell)) * sensitivity_heave(out.TEMP.initial_saturation(cell)) * 0.08;
    
                        else % equal
                            subsidence_poreice_cell = 0;
                        end
                        subsidence_poreice = subsidence_poreice + subsidence_poreice_cell;
                    end
                    
                    % total subsidence for current timestamp
                    if out.PARA.add_Xice_subsidence == 1
                        result = subsidence_xice + subsidence_poreice;
                    else
                        result = subsidence_poreice;
                    end
                    
                    % retired
                    % porespace = CURRENT.STATVAR.ice +  CURRENT.STATVAR.air;
                    % saturation = CURRENT.STATVAR.ice  ./ porespace;
                    % saturation(isnan(saturation)) = 0;
                    % current_pore_ice_heave = sum(porespace .* sensitivity_heave(saturation) ./ CURRENT.STATVAR.area .*0.08);
                    % subsidence_poreice = current_pore_ice_heave - out.TEMP.initial_pore_ice_heave; 
    
                else % year has not been initialized yet, initialize year
                    out.TEMP.current_year = year(tile.t);
                    
                    % initial ground surface elevation for subsidence due to xice melt
                    out.TEMP.initial_ground_surface_elevation = sum(CURRENT.STATVAR.layerThick) - CURRENT.STATVAR.XwaterIce(1,1) ./ CURRENT.STATVAR.area(1,1);
                    
                    % initial heave from pore ice
                    out.TEMP.initial_poreice = CURRENT.STATVAR.ice; %[m3]
                    air = (CURRENT.STATVAR.layerThick.*CURRENT.STATVAR.area) - CURRENT.STATVAR.waterIce -CURRENT.STATVAR.mineral -CURRENT.STATVAR.organic - CURRENT.STATVAR.XwaterIce;
                    air = max(0,air);
                    out.TEMP.initial_saturation = CURRENT.STATVAR.ice ./max(1e-100, CURRENT.STATVAR.ice+air);
                    out.TEMP.initial_saturation = max(0,min(1,out.TEMP.initial_saturation)); 

                    % observation for first timestamp
                    result = 0;

                    %retired code
                    %freeporespace = CURRENT.STATVAR.ice +  CURRENT.STATVAR.air;
                    %saturation = CURRENT.STATVAR.ice  ./ freeporespace;
                    %saturation(isnan(saturation)) = 0;
                    %out.TEMP.initial_pore_ice_heave = sum(CURRENT.STATVAR.ice .* sensitivity_heave(saturation) ./ CURRENT.STATVAR.area .*0.08);
                   

                end
                
                out.STATVAR.subsidence = [out.STATVAR.subsidence; result];
            end
            
            if tile.t>=out.SAVE_TIME
                out.TEMP.tag = ['_' out.PARA.tag '_' out.PARA.tag2 '_'];
                out.TEMP.tag = strrep(out.TEMP.tag, '___', '_');
                out.TEMP.tag = strrep(out.TEMP.tag, '__', '_');

                if ~(exist([tile.PARA.result_path tile.PARA.run_name])==7)
                    mkdir([tile.PARA.result_path tile.PARA.run_name])
                end
                CG_out.subsidence = out.STATVAR.subsidence;
                CG_out.timestamp = out.TIMESTAMP';
                CG_out.identifier = tile.RUN_INFO.PPROVIDER.PARA.identifier;
                
                save([tile.PARA.result_path tile.PARA.run_name '/' tile.PARA.run_name '_seasonalsubsidence' out.TEMP.tag datestr(tile.t,'yyyymmdd') '.mat'], 'CG_out')
                
                % Clear the out structure
                out.STATVAR.subsidence = [];
                out.TIMESTAMP=[];
                if ~isnan(out.PARA.save_interval)
                    out.SAVE_TIME = min(tile.FORCING.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(out.SAVE_TIME,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
                end
            end

            % sensitivity function heave(saturation) for pore ice
            function sens = sensitivity_heave(sat)
                sens = sat; % easiest form, 1-1 scaling with saturation
            end

        end

        function result = move_out2obs(out)
            result = out.STATVAR.subsidence;
        end
            
        
        function out = reset_new_stratigraphy(out, tile)
            
        end
        
    end
end