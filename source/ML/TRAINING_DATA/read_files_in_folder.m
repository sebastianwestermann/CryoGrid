classdef read_files_in_folder < matlab.mixin.Copyable
    
    properties
        PARA
        STATVAR
        CONST
        TEMP
        DATA
    end

    methods
        function features = provide_PARA(features)
            features.PARA.out_folder = []; 
            features.PARA.out_filename = [];
            features.PARA.out_date = [];
            features.PARA.start_year = [];
            features.PARA.end_year = [];

            features.PARA.variables = [];
        end


        function features = provide_CONST(features)

        end

        function features = provide_STATVAR(features)

        end

        function features = finalize_init(features, tile)

            str1 = num2str(features.PARA.out_date(1,1)+100);
            str2 = num2str(features.PARA.out_date(2,1)+100);
            features.PARA.out_str = [str1(2:3) str2(2:3)];

            valid_files = dir([features.PARA.out_folder  features.PARA.out_filename '_*' features.PARA.out_str '.mat']);

            file_index = [];
            year_list=[];
            length_of_filename = size(features.PARA.out_filename,2);
            for i=1:size(valid_files,1)
                next_pos = find(valid_files(i).name(length_of_filename+2:end) == '_', 1, 'first');
                file_index = [file_index; str2num(valid_files(i).name(length_of_filename+2:length_of_filename+2+next_pos-2))];
                year_list = [year_list; str2num(valid_files(i).name(length_of_filename+1+next_pos+1:length_of_filename+1+next_pos+4))];
            end

            features.TEMP.file_index = sort(unique(file_index)); %list of available points.

            load([features.PARA.out_folder valid_files(1).name])
            cell_size = CG_out.depths(1)-CG_out.depths(2);
            ground_surface_cell = find(sum(isnan(CG_out.T),2)==0, 1, 'first'); %bit nonsense if snow stays around!
            
            features.TEMP.target_cells = [];
            for i=1:size(features.PARA.variables,1)
                if strcmp('T', features.PARA.variables{i,1}(1,1))
                    target_depth = str2num(features.PARA.variables{i,1}(2:end-1));
                    features.TEMP.target_cells = [features.TEMP.target_cells; round(ground_surface_cell + target_depth./cell_size)];
                elseif strcmp('ALT', features.PARA.variables{i,1})
                    features.TEMP.get_ALT = 1;
                end
            end

        end

        function [out, data_groups] = generate_target_data(features, tile)
            disp('reading out-files')
            out = [];
            data_groups = [];

            for i=1:size(features.TEMP.file_index,1)
                for y = features.PARA.start_year:features.PARA.end_year
                    load([features.PARA.out_folder  features.PARA.out_filename '_' num2str(features.TEMP.file_index(i,1)) '_' num2str(y)  features.PARA.out_str '.mat'])
                    out_i = [];
                    for j=1:size(features.TEMP.target_cells,1)
                        out_i = [out_i mean(CG_out.T(features.TEMP.target_cells(j,1),:))];
                    end
                    out = [out; out_i];
                end
            end
            features.DATA.out = out;
            data_groups = [1:size(features.TEMP.target_cells,1)];
        end

    end
end

