

classdef OUT_MULTITILE_snow < matlab.mixin.Copyable

    properties
        
        PARA
        STATVAR
        CONST
        TEMP
        OUTPUT_TIME
        SAVE_TIME
    end
    
    methods

        function out = provide_PARA(out)
            out.PARA.output_timestep = [];
            out.PARA.save_date = [];
            out.PARA.save_interval = [];
            out.PARA.regrid_results = [];
            out.PARA.target_grid_size = [];
            out.PARA.tag = [];
        end
        
        function out = provide_CONST(out)
            out.CONST.day_sec = [];
        end
        
        function out = provide_STATVAR(out)

        end
        
        function out = finalize_init(out, tile)

            out = reset_STATVAR(out); %initializes STATVAR
            
            if isempty(out.PARA.save_interval) || isnan(out.PARA.save_interval) 
                out.SAVE_TIME = tile.FORCING.PARA.end_time;
            else
                if datenum([out.PARA.save_date datestr(tile.FORCING.PARA.start_time,'yyyy')], 'dd.mm.yyyy') <= tile.FORCING.PARA.start_time + 30
                    out.SAVE_TIME = min(tile.FORCING.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(tile.FORCING.PARA.start_time,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
                else
                    out.SAVE_TIME = min(tile.FORCING.PARA.end_time,  datenum([out.PARA.save_date datestr(tile.FORCING.PARA.start_time,'yyyy')], 'dd.mm.yyyy'));
                end
            end
            
            if isempty(out.PARA.output_timestep) || isnan(out.PARA.output_timestep)
                out.OUTPUT_TIME = out.SAVE_TIME;
            else
                out.OUTPUT_TIME = tile.FORCING.PARA.start_time + out.PARA.output_timestep;
            end
            
        end
        
        function out = store_OUT(out, tile)

            ground = tile.SUBSURFACE_CLASS;
            if tile.t >= out.OUTPUT_TIME

                disp(datestr(tile.t))
                snow_ice = ground.STATVAR.ice_snow(1:3,:);
                snow_ice(3,:) =  snow_ice(3,:) + ground.STATVAR.ice_snow(4,:);
                snow_water = ground.STATVAR.water_snow(1:3,:);
                snow_water(3,:) =  snow_water(3,:) + ground.STATVAR.water_snow(4,:);
                snow_water(snow_ice<1e-12) = NaN;
                snow_layerTick =  ground.STATVAR.layerThick_snow(1:3,:);
                snow_layerTick(3,:) =  snow_layerTick(3,:) + ground.STATVAR.layerThick_snow(4,:);
                snow_layerTick(snow_ice<1e-12) = NaN;
                thermCond = ground.STATVAR.thermCond_snow;
                thermCond(ground.STATVAR.ice_snow<1e-12) = 0;
                thermCond(3,:) =  thermCond(3,:) + thermCond(4,:);
                thermCond = thermCond(1:3,:);
                snow_ice(snow_ice<1e-12) = NaN;

                out.TEMP.snow_density = cat(3, out.TEMP.snow_density, snow_ice./snow_layerTick);
                out.TEMP.snow_ice = cat(3,out.TEMP.snow_ice, snow_ice);
                out.TEMP.snow_layerThick = cat(3, out.TEMP.snow_layerThick, snow_layerTick);
                out.TEMP.snow_water = cat(3, out.TEMP.snow_water, snow_water./snow_layerTick);
                out.TEMP.snow_conductivity  = cat(3, out.TEMP.snow_conductivity, thermCond);
                out.TEMP.snow_albedo = [out.TEMP.snow_albedo; ground.STATVAR.albedo];
                out.TEMP.snow_depth = [out.TEMP.snow_depth; sum(ground.STATVAR.layerThick_snow,1)];
                out.TEMP.snow_SWE = [out.TEMP.snow_SWE; sum(ground.STATVAR.ice_snow,1)];
                out.TEMP.snow_GST = [out.TEMP.snow_GST; ground.STATVAR.T(5,:)];
                out.TEMP.snow_MAAT = [out.TEMP.snow_MAAT; tile.FORCING.TEMP.Tair];
                out.TEMP.timestamp = [out.TEMP.timestamp tile.t];

                if tile.t>=out.SAVE_TIME

                    run_name = tile.PARA.run_name;
                    result_path = tile.PARA.result_path;
                    out_tag = out.PARA.tag;
                    t = tile.t;

                    if out.PARA.regrid_results == 1
                        CG_out = regrid_results_snow(out);
                    else
                        CG_out = out.TEMP;
                    end

                    if ~(exist([result_path run_name])==7)
                        mkdir([result_path run_name])
                    end
                    if isempty(out_tag) || all(isnan(out_tag))
                        save([result_path run_name '/' run_name '_' datestr(t,'yyyymmdd') '.mat'], 'CG_out')
                    else
                        save([result_path run_name '/' run_name '_' out_tag '_' datestr(t,'yyyymmdd') '.mat'], 'CG_out')
                    end


                    out = reset_STATVAR(out);
                    out.SAVE_TIME = min(tile.FORCING.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(out.SAVE_TIME,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
                end

                if isempty(out.PARA.output_timestep) || isnan(out.PARA.output_timestep)
                    out.OUTPUT_TIME = out.SAVE_TIME;
                else
                    out.OUTPUT_TIME = min(out.SAVE_TIME, out.OUTPUT_TIME + out.PARA.output_timestep);
                end
            end
        end

        function CG_out = regrid_results_snow(out)
            vars = fieldnames(out.TEMP);
            first_round = 1;
            for i=1:size(vars,1)
                if strcmp(vars{i,1}, 'snow_albedo') || strcmp(vars{i,1}, 'snow_depth') || strcmp(vars{i,1}, 'snow_SWE') || strcmp(vars{i,1}, 'snow_GST') || strcmp(vars{i,1}, 'snow_MAAT') || strcmp(vars{i,1}, 'timestamp')
                    CG_out.(vars{i,1}) = out.TEMP.(vars{i,1});
                elseif strcmp(vars{i,1}, 'snow_density') || strcmp(vars{i,1}, 'snow_water') || strcmp(vars{i,1}, 'snow_conductivity') || strcmp(vars{i,1}, 'snow_ice')
                    if first_round == 1 %establish new grid
                        lt = out.TEMP.snow_layerThick;
                        lt(isnan(lt))=0;
                        lt = sum(out.TEMP.snow_layerThick,1);
                        max_depth = max(lt(:));
                        max_depth = max(max_depth, 3);
                        new_grid = [0:out.PARA.target_grid_size:out.PARA.target_grid_size+max_depth]';
                        new_grid = (new_grid(1:end-1,1)+new_grid(2:end,1))./2;
                        CG_out.z = new_grid;
                        first_round = 0;
                    end

                    %regrid
                    CG_out.(vars{i,1}) = [];
                    for jj=1:size(out.TEMP.snow_layerThick,2) %dimension of ensemble
                        var_interp = [];
                        for j=1:size(out.TEMP.snow_layerThick,3)%dimension of time
                            lt = out.TEMP.snow_layerThick(:, jj, j);
                            target_var = out.TEMP.(vars{i,1})(:, jj, j);
                            delete_cells = find(lt(:,1)<=1e-3 | isnan(lt(:,1)));
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

                                var_interp  = [var_interp interp1(lt3, target_var2, new_grid, 'nearest')];
                            else
                                var_interp  = [var_interp new_grid.*NaN];
                            end
                        end
                        CG_out.(vars{i,1}) = cat(3, CG_out.(vars{i,1}), var_interp); %dimension of ensemble
                        
                    end
                end
            end
        end
            
        function out = reset_STATVAR(out)
            out.TEMP.snow_density = [];
            out.TEMP.snow_ice = [];
            out.TEMP.snow_layerThick = [];
            out.TEMP.snow_water = [];
            out.TEMP.snow_conductivity = []; 
            out.TEMP.snow_albedo = [];
            out.TEMP.snow_depth = [];
            out.TEMP.snow_SWE = [];
            out.TEMP.snow_GST = [];
            out.TEMP.snow_MAAT = [];
            out.TEMP.timestamp = [];
        end

    end
end

