classdef OUT_subsidence < OUT_BASE_OBSERVABLES
    
    properties

        STATVAR

    end
    
    methods

        function out = provide_PARA(out)
            out.PARA.output_timestep = [];
            out.PARA.save_date = [];
            out.PARA.save_interval = [];
            out.PARA.tag = [];
            out.PARA.tag2 = [];
            out.PARA.seasonal = [];
            out.PARA.add_Xice_subsidence=1; % only relevant for seasonal == 1

            out.PARA.timestamps = [];
            out.PARA.write_file_mode = 0; %empty/default setting: DA mode, do not write any output; % 0: normal behaviour, no DA, get normal output times and run regularly; 1: store after DA; 2: store after each run
        end
        
        
        function out = finalize_init(out, tile)
        
            out = finalize_init@OUT_BASE_OBSERVABLES(out, tile);

            out.STATVAR.subsidence = []; 
            out.STATVAR.timestamp = [];

            out.TEMP.initialized = 0;
            out.TEMP.initial_ground_surface_elevation = 0;
            out.TEMP.current_year = 1e20; % random unrealistic number, will be set to current year after first annual observation
            out.TEMP.initial_poreice = 0; %[m3]
            out.TEMP.initial_saturation =0;

            out.TEMP.keyword = 'subsidence';
            out.TEMP.output_var = 'CG_subsidence';

        end

        function out = store_OUT(out, tile)

            if tile.t >= out.OUTPUT_TIME

                disp([datestr(tile.t)])

                out = state2out(out, tile);

                out = set_output_time(out, tile);
                
                if tile.t >= out.SAVE_TIME && isempty(out.PARA.timestamps)
                    if out.PARA.write_file_mode == 0
                        out = out2file_CG(out, tile);
                    elseif out.PARA.write_file_mode == 1
                        out.SAVE_TIME = min(tile.FORCING.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(out.SAVE_TIME,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
                        out.TEMP.write_out_init = 1;
                    elseif out.PARA.write_file_mode == 2
                        out.SAVE_TIME = min(tile.FORCING.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(out.SAVE_TIME,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
                        out.TEMP.write_out_final = 1;
                    end
                end
                
            end
        end

        function out = out2file_CG(out, tile)

            out = out2file_CG@OUT_BASE_OBSERVABLES(out, tile);
            out.STATVAR.timestamp = [];
            out.STATVAR.subsidence = [];

        end

       
        function out = state2out(out,tile)

            result = 0;
            CURRENT = tile.TOP.NEXT;

            while ~is_ground_surface(CURRENT) && ~(strcmp(class(CURRENT), 'Bottom'))
                CURRENT = CURRENT.NEXT;
            end

            if out.PARA.seasonal == 1
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

                else % year has not been initialized yet, initialize year
                    out.TEMP.current_year = year(tile.t);
                    out.TEMP.initialized = 1;

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
                end

            else % out.PARA.seasonal == 0
                if ~out.TEMP.initialized

                    out.TEMP.initialized = 1;
                    out.TEMP.initial_ground_surface_elevation = sum(CURRENT.STATVAR.layerThick) - CURRENT.STATVAR.XwaterIce(1,1) ./ CURRENT.STATVAR.area(1,1);

                    result = 0;
                else
                    ground_surface_elevation = sum(CURRENT.STATVAR.layerThick) - CURRENT.STATVAR.XwaterIce(1,1) ./ CURRENT.STATVAR.area(1,1);

                    result = ground_surface_elevation - out.TEMP.initial_ground_surface_elevation;
                end

            end
            out.STATVAR.subsidence = [out.STATVAR.subsidence; result];
            %ADDED SW
            out.STATVAR.timestamp = [out.STATVAR.timestamp; tile.t];

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
        

        % function out = store_OUT(out, tile)
        % 
        %     if tile.t >= out.OUTPUT_TIME
        % 
        %         out.TIMESTAMP = [out.TIMESTAMP tile.t];
        % 
        %         if ~isempty(out.PARA.timestamps)
        %             out.OUTPUT_TIME = out.PARA.timestamps(find(out.PARA.timestamps(:,1)-tile.t>0, 1, 'first'), 1);
        %         else
        %             out.OUTPUT_TIME = min(out.SAVE_TIME, out.OUTPUT_TIME + out.PARA.output_timestep);
        %         end
        %         if isempty(out.OUTPUT_TIME)
        %             out.OUTPUT_TIME = Inf;
        %         end
        % 
        % 
        %         %% calc section
        %         result = 0;
        %         CURRENT = tile.TOP.NEXT;
        % 
        %         while ~is_ground_surface(CURRENT) && ~(strcmp(class(CURRENT), 'Bottom'))
        %             CURRENT = CURRENT.NEXT;
        %         end
        % 
        %         if out.PARA.seasonal == 1
        %             if year(tile.t) == out.TEMP.current_year  % for all observations after the first (= 0) of this year
        % 
        %                 % subsidence from xice
        %                 ground_surface_elevation = sum(CURRENT.STATVAR.layerThick) - CURRENT.STATVAR.XwaterIce(1,1) ./ CURRENT.STATVAR.area(1,1);
        %                 subsidence_xice = ground_surface_elevation - out.TEMP.initial_ground_surface_elevation;
        % 
        %                 % subsidence from pore ice melt
        %                 poreice = CURRENT.STATVAR.ice;
        %                 air = (CURRENT.STATVAR.layerThick.*CURRENT.STATVAR.area) - CURRENT.STATVAR.waterIce -CURRENT.STATVAR.mineral -CURRENT.STATVAR.organic - CURRENT.STATVAR.XwaterIce; % air is not updated but htere is a field, calculate new
        %                 air = max(0,air);
        %                 saturation = CURRENT.STATVAR.ice ./max(1e-100,CURRENT.STATVAR.ice+air);
        %                 saturation = max(0,min(1,saturation));
        %                 subsidence_poreice=0;
        %                 for cell=1:size(poreice,1)
        %                     if poreice(cell) > out.TEMP.initial_poreice(cell)
        %                         subsidence_poreice_cell = (poreice(cell) - out.TEMP.initial_poreice(cell)) * sensitivity_heave(saturation(cell)) * 0.08;
        % 
        %                     elseif poreice(cell) < out.TEMP.initial_poreice(cell)
        %                         subsidence_poreice_cell = (poreice(cell) - out.TEMP.initial_poreice(cell)) * sensitivity_heave(out.TEMP.initial_saturation(cell)) * 0.08;
        % 
        %                     else % equal
        %                         subsidence_poreice_cell = 0;
        %                     end
        %                     subsidence_poreice = subsidence_poreice + subsidence_poreice_cell;
        %                 end
        % 
        %                 % total subsidence for current timestamp
        %                 if out.PARA.add_Xice_subsidence == 1
        %                     result = subsidence_xice + subsidence_poreice;
        %                 else
        %                     result = subsidence_poreice;
        %                 end
        % 
        %             else % year has not been initialized yet, initialize year
        %                 out.TEMP.current_year = year(tile.t);
        %                 out.TEMP.initialized = 1;
        % 
        %                 % initial ground surface elevation for subsidence due to xice melt
        %                 out.TEMP.initial_ground_surface_elevation = sum(CURRENT.STATVAR.layerThick) - CURRENT.STATVAR.XwaterIce(1,1) ./ CURRENT.STATVAR.area(1,1);
        % 
        %                 % initial heave from pore ice
        %                 out.TEMP.initial_poreice = CURRENT.STATVAR.ice; %[m3]
        %                 air = (CURRENT.STATVAR.layerThick.*CURRENT.STATVAR.area) - CURRENT.STATVAR.waterIce -CURRENT.STATVAR.mineral -CURRENT.STATVAR.organic - CURRENT.STATVAR.XwaterIce;
        %                 air = max(0,air);
        %                 out.TEMP.initial_saturation = CURRENT.STATVAR.ice ./max(1e-100, CURRENT.STATVAR.ice+air);
        %                 out.TEMP.initial_saturation = max(0,min(1,out.TEMP.initial_saturation)); 
        % 
        %                 % observation for first timestamp
        %                 result = 0;    
        %             end
        % 
        %         else % out.PARA.seasonal == 0
        %             if ~out.TEMP.initialized
        % 
        %                 out.TEMP.initialized = 1;
        %                 out.TEMP.initial_ground_surface_elevation = sum(CURRENT.STATVAR.layerThick) - CURRENT.STATVAR.XwaterIce(1,1) ./ CURRENT.STATVAR.area(1,1);
        % 
        %                 result = 0;
        %             else
        %                 ground_surface_elevation = sum(CURRENT.STATVAR.layerThick) - CURRENT.STATVAR.XwaterIce(1,1) ./ CURRENT.STATVAR.area(1,1);
        % 
        %                 result = ground_surface_elevation - out.TEMP.initial_ground_surface_elevation;
        %             end
        % 
        %         end
        %         out.STATVAR.subsidence = [out.STATVAR.subsidence; result];
        %     end
        % 
        %     if tile.t>=out.SAVE_TIME
        %         out.TEMP.tag = ['_' out.PARA.tag '_' out.PARA.tag2 '_'];
        %         out.TEMP.tag = strrep(out.TEMP.tag, '___', '_');
        %         out.TEMP.tag = strrep(out.TEMP.tag, '__', '_');
        % 
        %         if ~(exist([tile.PARA.result_path tile.PARA.run_name])==7)
        %             mkdir([tile.PARA.result_path tile.PARA.run_name])
        %         end
        %         CG_out.subsidence = out.STATVAR.subsidence;
        %         CG_out.timestamp = out.TIMESTAMP';
        % 
        %         if out.PARA.seasonal == 1
        %             save([tile.PARA.result_path tile.PARA.run_name '/' tile.PARA.run_name '_seasonalsubsidence' out.TEMP.tag datestr(tile.t,'yyyymmdd') '.mat'], 'CG_out')
        %         else % out.PARA.seasonal == 0 / interannual
        %             save([tile.PARA.result_path tile.PARA.run_name '/' tile.PARA.run_name '_interannualsubsidence' out.TEMP.tag datestr(tile.t,'yyyymmdd') '.mat'], 'CG_out')
        %         end
        % 
        % 
        %         % Clear the out structure
        %         out.STATVAR.subsidence = [];
        %         out.TIMESTAMP=[];
        %         if ~isnan(out.PARA.save_interval)
        %             out.SAVE_TIME = min(tile.FORCING.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(out.SAVE_TIME,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
        %         end
        %     end
        % 
        %     % sensitivity function heave(saturation) for pore ice
        %     function sens = sensitivity_heave(sat)
        %         sens = sat; % easiest form, 1-1 scaling with saturation
        %     end
        % 
        % end

    end
end
