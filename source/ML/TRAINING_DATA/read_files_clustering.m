classdef read_files_clustering < matlab.mixin.Copyable
    
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

            load([features.PARA.out_folder  features.PARA.out_filename '_' num2str(tile.RUN_INFO.CLUSTER.STATVAR.key_centroid_index(1,1)) '_' num2str(features.PARA.start_year)  features.PARA.out_str '.mat'])
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

        function out = generate_target_data(features, tile)
            disp('reading out-files')
            out = [];
            valid_ind = sort(tile.RUN_INFO.CLUSTER.STATVAR.sample_centroid_index);
            for i=1:size(valid_ind,1)
                for y = features.PARA.start_year:features.PARA.end_year
                    load([features.PARA.out_folder  features.PARA.out_filename '_' num2str(valid_ind(i,1)) '_' num2str(y)  features.PARA.out_str '.mat'])
                    out_i = [];
                    for j=1:size(features.TEMP.target_cells,1)
                        out_i = [out_i mean(CG_out.T(features.TEMP.target_cells(j,1),:))];
                    end
                    out = [out; out_i];
                end
            end
            features.DATA.out = out;
        end

    end
end

