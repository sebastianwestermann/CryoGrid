%========================================================================
% CryoGrid OUT class OUT_snow_all delivering simple snow depths, in addition to
% important snow parameters (density, sphericity, grain size, water
% content, T below snow) - designed as output class for ITCH 2025 SnowMIP
% works together with an Xice class so Xice layers forming at the surface are included in the results 
% S. Westermann, Nov 2025
%========================================================================


classdef OUT_snow_all_OPT < OUT_BASE_OBSERVABLES & OUT_snow_all
 

    properties

    end
    
    
    methods
    
        function out = provide_PARA(out)         

            out.PARA.variables = []; %depth density, albedo, water, sphericity, grain size, dendricity, conductivity
            out.PARA.output_timestep = [];
            out.PARA.target_grid_size = [];
            out.PARA.include_Xice = [];
            out.PARA.regrid_results = [];
            out.PARA.save_date = [];
            out.PARA.save_interval = [];
            out.PARA.tag = [];
            out.PARA.tag2 = [];
            out.PARA.write_file_mode = 0; %empty/default setting: DA mode, do not write any output; % 0: normal behaviour, no DA, get normal output times and run regularly; 1: store after DA; 2: store after each run
            out.PARA.timestamps = [];
        end

        
        function out = finalize_init(out, tile)

            out = finalize_init@OUT_BASE_OBSERVABLES(out, tile);

            out.STATVAR.timestamp = [];
            out.PARA.variables = [out.PARA.variables; 'layerThick'];

            for i=1:size(out.PARA.variables,1)
                str2 = out.PARA.variables{i,1}; %['snow_' out.PARA.variables{i,1}];
                out.STATVAR.(str2) = {};
            end
            out.TEMP.keyword = 'snow';
            out.TEMP.output_var = 'CG_snow';
        end
        
        function out = store_OUT(out, tile)

            if tile.t >= out.OUTPUT_TIME
                
                disp([datestr(tile.t)])

                out = state2out(out, tile);

                out = set_output_time(out, tile);
                
                if tile.t >= out.SAVE_TIME
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

            if out.PARA.regrid_results == 1
                out = regrid_results_snow(out);
            end
            out = out2file_CG@OUT_BASE_OBSERVABLES(out, tile);
            
            out.STATVAR.timestamp = [];
            for i=1:size(out.PARA.variables,1)
                out.STATVAR.(out.PARA.variables{i,1}) = {};
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