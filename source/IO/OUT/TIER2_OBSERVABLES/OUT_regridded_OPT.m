%========================================================================
% CryoGrid OUT class OUT_all
% CryoGrid OUT class defining storage format of the output 
% OUT_all stores identical copies of all GROUND classses (including STATVAR, TEMP, PARA) in the
% stratigraphy for each output timestep, while lateral classes are not stored.
% The user can specify the save date and the save interval (e.g. yearly
% files), as well as the output timestep (e.g. 6 hourly). The output files
% are in Matlab (".mat") format.
% S. Westermann, T. Ingeman-Nielsen, J. Scheer, June 2021
%========================================================================


classdef OUT_regridded_OPT < OUT_BASE_OBSERVABLES & OUT_regridded
 

    properties
    
    end
    
    
    methods
        
        
        function out = provide_PARA(out)         

            out.PARA.variables = [];
            out.PARA.upper_elevation = [];
            out.PARA.lower_elevation = [];
            out.PARA.relative2surface = [];
            out.PARA.target_grid_size = [];
            out.PARA.output_timestep = [];
            out.PARA.save_date = [];
            out.PARA.save_interval = [];
            out.PARA.tag = [];
            out.PARA.tag2 = [];
            out.PARA.write_file_mode = 0; %empty/default setting: DA mode, do not write any output; % 0: normal behaviour, no DA, get normal output times and run regularly; 1: store after DA; 2: store after each run
            out.PARA.timestamps = [];
        end

        
        
        function out = finalize_init(out, tile)

            out = finalize_init@OUT_BASE_OBSERVABLES(out, tile);

            out.TEMP.variables = [out.PARA.variables; 'class_number'];
            for i = 1:size(out.TEMP.variables,1)
               out.STATVAR.(out.TEMP.variables{i,1}) = []; 
            end
            out.STATVAR.timestamp = [];
            out.STATVAR.depths = [];
            
            if ~isempty(out.PARA.relative2surface) && ~isnan(out.PARA.relative2surface) && out.PARA.relative2surface
                out.PARA.upper_elevation = tile.PARA.altitude + out.PARA.upper_elevation;
                out.PARA.lower_elevation = tile.PARA.altitude - out.PARA.lower_elevation;
            end
           
            out.TEMP.output_var = 'CG_ground';
            out.TEMP.keyword = 'ground';

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



        %-------------param file generation-----
        function out = param_file_info(out)
            out = provide_PARA(out);

            out.PARA.STATVAR = [];
            out.PARA.options = [];
            out.PARA.class_category = 'OUT';
            
            out.PARA.comment.variables = {'select output variables: T, water, ice, waterIce, XwaterIce, Xwater, Xice, saltConc are supported'};
            out.PARA.options.variables.name = 'H_LIST';
            out.PARA.options.variables.entries_x = {'T' 'water' 'ice' 'waterIce'};
      
            out.PARA.comment.upper_elevation = {'upper elevation of output domain, in m a.s.l.; must match altitude in TILE!'};
            
            out.PARA.comment.lower_elevation = {'lower elevation of output domain, in m a.s.l.; must match altitude in TILE!'};
            
            out.PARA.default_value.target_grid_size = {0.02};
            out.PARA.comment.target_grid_size = {'cell size of regular grid that values are interpolated to'};
           
            out.PARA.default_value.output_timestep = {0.25};
            out.PARA.comment.output_timestep = {'timestep of output [days]'};

            out.PARA.default_value.save_date = {'01.09.'};
            out.PARA.comment.save_date = {'date (dd.mm.) when output file is written'};
            
            out.PARA.default_value.save_interval = {1};
            out.PARA.comment.save_interval = {'interval of output files [years]'};
            
            out.PARA.default_value.tag = {''};
            out.PARA.comment.tag = {'additional tag added to file name'};
        end
        

    end
end