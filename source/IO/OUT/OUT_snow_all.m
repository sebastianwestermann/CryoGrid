%========================================================================
% CryoGrid OUT class OUT_snow_all delivering simple snow depths, in addition to
% important snow parameters (density, sphericity, grain size, water
% content, T below snow) - designed as output class for ITCH 2025 SnowMIP
% works together with an Xice class so Xice layers forming at the surface are included in the results 
% S. Westermann, Nov 2025
%========================================================================


classdef OUT_snow_all < matlab.mixin.Copyable
 

    properties
        MISC
        TEMP
        PARA
        OUTPUT_TIME
        SAVE_TIME
        CONST
        result
        timestamp
        
    end
    
    
    methods
        
        %initialization
        
        function out = provide_PARA(out)         

            out.PARA.variables = []; %depth density, albedo, water, sphericity, grain size, dendricity, conductivity
            out.PARA.output_timestep = [];
            out.PARA.target_grid_size = [];
            out.PARA.include_Xice = [];
            out.PARA.regrid_results = [];
            out.PARA.save_date = [];
            out.PARA.save_interval = [];
            out.PARA.tag = [];
        end

        
        function out = provide_CONST(out)

        end

        
        function out = provide_STATVAR(out)

        end

        
        function out = finalize_init(out, tile)

            forcing = tile.FORCING;
            
            % Set the next (first) output time. This is the next (first) time output
            % is collected (in memory) for later storage to disk.
            out.OUTPUT_TIME = forcing.PARA.start_time + out.PARA.output_timestep;
            out.timestamp = [];
            out.PARA.variables = [out.PARA.variables; 'layerThick'];

            for i=1:size(out.PARA.variables,1)
                str2 = ['snow_' out.PARA.variables{i,1}];
                out.result.(str2) = {};
            end

            % Set the next (first) save time. This is the next (first) time all the
            % collected output is saved to disk.
            if isempty(out.PARA.save_interval) || isnan(out.PARA.save_interval) 
                out.SAVE_TIME = forcing.PARA.end_time;
            else
                out.SAVE_TIME = min(forcing.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(forcing.PARA.start_time,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
            end
                        
        end
        
        %---------------time integration-------------
                    
        function out = store_OUT(out, tile)           
            
            if tile.t >= out.OUTPUT_TIME

                t = tile.t;
                disp(datestr(t))
                TOP = tile.TOP;
                BOTTOM = tile.BOTTOM;
                forcing = tile.FORCING;
                run_name = tile.PARA.run_name; %tile.RUN_NUMBER;
                result_path = tile.PARA.result_path;
                timestep = tile.timestep;
                out_tag = out.PARA.tag;
            
                % Store the current state of the model in the out structure.

                out.timestamp = [out.timestamp t];
                
                snow_depth = 0;
                CURRENT = TOP.NEXT;
                if isprop(CURRENT, 'CHILD') && CURRENT.CHILD ~= 0 %snow is CHILD
                    class_name = class(CURRENT.CHILD);
                    if strcmp(class_name(1:4), 'SNOW')
                        for i=1:size(out.PARA.variables,1)
                            str1 = ['get_snow_' out.PARA.variables{i,1} '_CHILD'];
                            str2 = ['snow_' out.PARA.variables{i,1}];
                            a = str2func(str1);
                            out.result.(str2) = [out.result.(str2); a(out, CURRENT, tile)];
                        end
                    end
                else
                    class_name = class(CURRENT);
                    if strcmp(class_name(1:4), 'SNOW') %snow is normal class
                        for i=1:size(out.PARA.variables,1)
                            str1 = ['get_snow_' out.PARA.variables{i,1}];
                            str2 = ['snow_' out.PARA.variables{i,1}];
                            a = str2func(str1);
                            out.result.(str2) = [out.result.(str2); a(out, CURRENT, tile)];
                        end
                    else
                        if out.PARA.include_Xice == 1
                            for i=1:size(out.PARA.variables,1)
                                str1 = ['get_snow_' out.PARA.variables{i,1} '_Xice'];
                                str2 = ['snow_' out.PARA.variables{i,1}];
                                a = str2func(str1);
                                out.result.(str2) = [out.result.(str2); a(out, CURRENT, tile)];
                            end
                        else
                            for i=1:size(out.PARA.variables,1)
                                str2 = ['snow_' out.PARA.variables{i,1}];
                                if strcmp(out.PARA.variables{i,1}, 'GST')
                                    str1 = ['get_snow_' out.PARA.variables{i,1} '_Xice'];
                                    a = str2func(str1);
                                    out.result.(str2) = [out.result.(str2); a(out, CURRENT, tile)];
                                else
                                    out.result.(str2) = [out.result.(str2); NaN];
                                end
                            end
                        end
                    end
                end

                
                % Set the next OUTPUT_TIME
                out.OUTPUT_TIME = min(out.SAVE_TIME, out.OUTPUT_TIME + out.PARA.output_timestep);
                
                if t>=out.SAVE_TIME
                    % It is time to save all the collected model output to disk
                     
                    if ~(exist([result_path run_name])==7)
                        mkdir([result_path run_name])
                    end

                    if out.PARA.regrid_results == 1
                        CG_out = regrid_results_snow(out);
                        CG_out.timestamp = out.timestamp;
                    else
                        CG_out = out.result;
                        CG_out.timestamp = out.timestamp;
                    end
                    if isempty(out_tag) || all(isnan(out_tag))
                        save([result_path run_name '/' run_name '_snow_' datestr(t,'yyyymmdd') '.mat'], 'CG_out')
                    else
                        save([result_path run_name '/' run_name '_snow_' out_tag '_' datestr(t,'yyyymmdd') '.mat'], 'CG_out')
                    end
                    
                    % Clear the out structure
                    out.timestamp = [];
                    for i=1:size(out.PARA.variables,1)
                        str2 = ['snow_' out.PARA.variables{i,1}];
                        out.result.(str2) = {};
                    end

                    if ~isnan(out.PARA.save_interval)
                        % If save_interval is defined, uptate SAVE_TIME for next save opertion 
                        out.SAVE_TIME = min(forcing.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(out.SAVE_TIME,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
                        % If save_interval is not defined, we will save at the very end of the model run
					end
                end
            end
        end
        
        %normal snow
        function res = get_snow_MAAT(out, CURRENT, tile)
            res = tile.FORCING.TEMP.Tair;
        end

        function res = get_snow_GST(out, CURRENT, tile)
            res = CURRENT.NEXT.STATVAR.T(1,1);
        end

        function res = get_snow_depth(out, CURRENT, tile)
            res = sum(CURRENT.STATVAR.layerThick);
            if out.PARA.include_Xice == 1
                res = res + CURRENT.NEXT.STATVAR.Xice(1,1) ./ CURRENT.NEXT.STATVAR.area(1,1);
            end
        end

        function res = get_snow_SWE(out, CURRENT, tile)
            res = sum(CURRENT.STATVAR.ice./ CURRENT.STATVAR.area);
            if out.PARA.include_Xice == 1
                res = res + CURRENT.NEXT.STATVAR.Xice(1,1) ./ CURRENT.NEXT.STATVAR.area(1,1);
            end
        end

        function res = get_snow_albedo(out, CURRENT, tile)
            res = CURRENT.STATVAR.albedo;
        end

        function res = get_snow_layerThick(out, CURRENT, tile)
            res =  CURRENT.STATVAR.layerThick;
            if out.PARA.include_Xice == 1
                res = [res; CURRENT.NEXT.STATVAR.Xice(1,1) ./ CURRENT.NEXT.STATVAR.area(1,1)];
            end
        end

        function res = get_snow_density(out, CURRENT, tile)
            res = CURRENT.STATVAR.ice ./ (CURRENT.STATVAR.layerThick .* CURRENT.STATVAR.area);
            if out.PARA.include_Xice == 1
                res = [res; 1];
            end
        end

        function res = get_snow_water(out, CURRENT, tile)
            res = CURRENT.STATVAR.water ./ (CURRENT.STATVAR.layerThick .* CURRENT.STATVAR.area);
            if out.PARA.include_Xice == 1
                res = [res; 0];
            end
        end

        function res = get_snow_sphericity(out, CURRENT, tile)
            res = CURRENT.STATVAR.s;
            if out.PARA.include_Xice == 1
                res = [res; -1];
            end
        end

        function res = get_snow_grainSize(out, CURRENT, tile)
            res = CURRENT.STATVAR.gs;
            if out.PARA.include_Xice == 1
                res = [res; -1];
            end
        end

        function res = get_snow_dendricity(out, CURRENT, tile)
            res = CURRENT.STATVAR.d;
            if out.PARA.include_Xice == 1
                res = [res; -1];
            end
        end

        function res = get_snow_conductivity(out, CURRENT, tile)
            res = CURRENT.STATVAR.thermCond;
            if out.PARA.include_Xice == 1
                res = [res; 2.2];
            end
        end

        %CHILD functions
        function res = get_snow_MAAT_CHILD(out, CURRENT, tile)
            res = tile.FORCING.TEMP.Tair;
        end

        function res = get_snow_GST_CHILD(out, CURRENT, tile)
            res = CURRENT.STATVAR.T(1,1);
        end

        function res = get_snow_depth_CHILD(out, CURRENT, tile)
            res = CURRENT.CHILD.STATVAR.area .* CURRENT.CHILD.STATVAR.layerThick ./ CURRENT.STATVAR.area(1,1);
            if out.PARA.include_Xice == 1
                res = res + CURRENT.STATVAR.Xice(1,1) ./ CURRENT.STATVAR.area(1,1);
            end
        end

        function res = get_snow_SWE_CHILD(out, CURRENT, tile)
            res = CURRENT.CHILD.STATVAR.ice ./ CURRENT.STATVAR.area(1,1);
            if out.PARA.include_Xice == 1
                res = res + CURRENT.STATVAR.Xice(1,1) ./ CURRENT.STATVAR.area(1,1);
            end
        end

        function res = get_snow_albedo_CHILD(out, CURRENT, tile)
            res = CURRENT.CHILD.STATVAR.albedo;
        end

        function res = get_snow_layerThick_CHILD(out, CURRENT, tile)
            res =  CURRENT.CHILD.STATVAR.area .* CURRENT.CHILD.STATVAR.layerThick ./ CURRENT.STATVAR.area(1,1);
            if out.PARA.include_Xice == 1
                res = [res; CURRENT.STATVAR.Xice(1,1) ./ CURRENT.STATVAR.area(1,1)];
            end
        end

        function res = get_snow_density_CHILD(out, CURRENT, tile)
            res =  CURRENT.CHILD.STATVAR.ice ./ (CURRENT.CHILD.STATVAR.layerThick .* CURRENT.CHILD.STATVAR.area);
            if out.PARA.include_Xice == 1
                res = [res; 1];
            end
        end

        function res = get_snow_water_CHILD(out, CURRENT, tile)
            res =  CURRENT.CHILD.STATVAR.water ./ (CURRENT.CHILD.STATVAR.layerThick .* CURRENT.CHILD.STATVAR.area);
            if out.PARA.include_Xice == 1
                res = [res; 0];
            end
        end

        function res = get_snow_sphericity_CHILD(out, CURRENT, tile)
            res =  CURRENT.CHILD.STATVAR.s;
            if out.PARA.include_Xice == 1
                res = [res; -1];
            end
        end

        function res = get_snow_grainSize_CHILD(out, CURRENT, tile)
            res =  CURRENT.CHILD.STATVAR.gs;
            if out.PARA.include_Xice == 1
                res = [res; -1];
            end
        end

        function res = get_snow_dendricity_CHILD(out, CURRENT, tile)
            res =  CURRENT.CHILD.STATVAR.d;
            if out.PARA.include_Xice == 1
                res = [res; -1];
            end
        end

        function res = get_snow_conductivity_CHILD(out, CURRENT, tile)
            res =  CURRENT.CHILD.STATVAR.thermCond;
            if out.PARA.include_Xice == 1
                res = [res; 2.2];
            end
        end

        %Xice
        function res = get_snow_MAAT_Xice(out, CURRENT, tile)
            res = tile.FORCING.TEMP.Tair;
        end

        function res = get_snow_GST_Xice(out, CURRENT, tile)
            res = CURRENT.STATVAR.T(1,1);
        end

        function res = get_snow_depth_Xice(out, CURRENT, tile)
            res = CURRENT.STATVAR.Xice(1,1) ./ CURRENT.NEXT.STATVAR.area(1,1);
        end

        function res = get_snow_SWE_Xice(out, CURRENT, tile)
            res = CURRENT.STATVAR.Xice(1,1) ./ CURRENT.NEXT.STATVAR.area(1,1);
        end

        function res = get_snow_albedo_Xice(out, CURRENT, tile)
            res = CURRENT.PARA.albedo;
        end

        function res = get_snow_layerThick_Xice(out, CURRENT, tile)
            res = CURRENT.STATVAR.Xice(1,1) ./ CURRENT.STATVAR.area(1,1);
        end

        function res = get_snow_density_Xice(out, CURRENT, tile)
            res = 1;
        end

        function res = get_snow_water_Xice(out, CURRENT, tile)
            res = 0;
        end

        function res = get_snow_sphericity_Xice(out, CURRENT, tile)
            res = -1;
        end

        function res = get_snow_grainSize_Xice(out, CURRENT, tile)
            res = -1;
        end

        function res = get_snow_dendricity_Xice(out, CURRENT, tile)
            res = -1;
        end

        function res = get_snow_conductivity_Xice(out, CURRENT, tile)
            res = 2.2;
        end

        function CG_out = regrid_results_snow(out)
            vars = fieldnames(out.result);
            first_round = 1;
            for i=1:size(vars,1)
                if strcmp(vars{i,1}, 'snow_albedo') || strcmp(vars{i,1}, 'snow_depth') || strcmp(vars{i,1}, 'snow_SWE') || strcmp(vars{i,1}, 'snow_GST') || strcmp(vars{i,1}, 'snow_MAAT')
                    CG_out.(vars{i,1}) = cell2mat(out.result.(vars{i,1}));
                elseif strcmp(vars{i,1}, 'snow_density') || strcmp(vars{i,1}, 'snow_water') || strcmp(vars{i,1}, 'snow_sphericity') || strcmp(vars{i,1}, 'snow_grainSize') || strcmp(vars{i,1}, 'snow_dendricity') || strcmp(vars{i,1}, 'snow_conductivity')
                    if first_round == 1 %establish new grid
                        max_depth = 0;
                        lt_all = out.result.snow_layerThick;
                        for j=1:size(lt_all,1)
                            max_depth = max(max_depth, sum(lt_all{j,1}));
                        end
                        new_grid = [0:out.PARA.target_grid_size:out.PARA.target_grid_size+max_depth]';
                        new_grid = (new_grid(1:end-1,1)+new_grid(2:end,1))./2;
                        CG_out.z = new_grid;
                        first_round = 0;
                    end

                    %regrid
                    CG_out.(vars{i,1}) = [];
                    for j=1:size(lt_all,1)
                        lt = lt_all{j,1};
                        target_var = out.result.(vars{i,1}){j,1};
                        delete_cells = find(lt(:,1)<=1e-4);
                        target_var(delete_cells,:) = [];

                        lt(delete_cells,:) = [];
                        if sum(lt)>0
                            if size(lt,1) == 1
                                lt3=[0; lt];
                                target_var2 = [target_var; target_var];
                            else
                                delta=0.001;
                                lt2=[0; lt(1)-delta; 2.*delta];
                                for ii=2:size(lt,1)-1
                                    lt2=[lt2; lt(ii)-2.*delta; 2.*delta];
                                end
                                lt2 =[lt2; lt(end)-delta];
                                lt3=cumsum(lt2);
                                target_var2 = [];
                                for ii=1:size(target_var,1)
                                    target_var2 = [target_var2; target_var(ii); target_var(ii)];
                                end
                            end
                            lt3 = -(lt3-lt3(end));

                            % lt = [0; lt];
                            % lt2=cumsum(lt);
                            % lt2=[lt2(1); (lt2(1:end-1)+lt2(2:end))/2; lt2(end)];
                            % lt2 = -(lt2-lt2(end));
                            %target_var = [target_var(1); target_var; target_var(end)];

                            CG_out.(vars{i,1}) = [CG_out.(vars{i,1}) interp1(lt3, target_var2, new_grid, 'nearest')];
                        else
                            CG_out.(vars{i,1}) = [CG_out.(vars{i,1}) new_grid.*NaN];
                        end
                    end
                    if strcmp(vars{i,1}, 'snow_density')
                        CG_out.snow_ice =  CG_out.snow_density;
                        CG_out.snow_density = CG_out.snow_density .* 920;
                    elseif strcmp(vars{i,1}, 'snow_sphericity') || strcmp(vars{i,1}, 'snow_grainSize') || strcmp(vars{i,1}, 'snow_dendricity')
                        CG_out.(vars{i,1})(find(CG_out.(vars{i,1})==-1)) = NaN;
                    end
                end
            end
        end


        
        % 
        % %-------------param file generation-----
        % function out = param_file_info(out)
        %     out = provide_PARA(out);
        % 
        %     out.PARA.STATVAR = [];
        %     out.PARA.options = [];
        %     out.PARA.class_category = 'OUT';
        % 
        %     out.PARA.comment.variables = {'select output variables: T, water, ice, waterIce, XwaterIce, Xwater, Xice, saltConc are supported'};
        %     out.PARA.options.variables.name = 'H_LIST';
        %     out.PARA.options.variables.entries_x = {'T' 'water' 'ice' 'waterIce'};
        % 
        %     out.PARA.comment.upper_elevation = {'upper elevation of output domain, in m a.s.l.; must match altitude in TILE!'};
        % 
        %     out.PARA.comment.lower_elevation = {'lower elevation of output domain, in m a.s.l.; must match altitude in TILE!'};
        % 
        %     out.PARA.default_value.target_grid_size = {0.02};
        %     out.PARA.comment.target_grid_size = {'cell size of regular grid that values are interpolated to'};
        % 
        %     out.PARA.default_value.output_timestep = {0.25};
        %     out.PARA.comment.output_timestep = {'timestep of output [days]'};
        % 
        %     out.PARA.default_value.save_date = {'01.09.'};
        %     out.PARA.comment.save_date = {'date (dd.mm.) when output file is written'};
        % 
        %     out.PARA.default_value.save_interval = {1};
        %     out.PARA.comment.save_interval = {'interval of output files [years]'};
        % 
        %     out.PARA.default_value.tag = {''};
        %     out.PARA.comment.tag = {'additional tag added to file name'};
        % end
        

    end
end