classdef read_files_clustering_GT < matlab.mixin.Copyable
    
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

            features.PARA.depths = [];
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
            for i=1:size(features.PARA.depths,1)
                target_depth = features.PARA.depths(i,1);
                features.TEMP.target_cells = [features.TEMP.target_cells; round(ground_surface_cell + target_depth./cell_size)];
            end

        end

        function [out, data_groups]  = generate_target_data(features, tile)
            disp('reading out-files')
            out = [];
            valid_ind = sort(tile.RUN_INFO.CLUSTER.STATVAR.key_centroid_index);
            for i=1:size(valid_ind,1)
                out_i = [];
                for y = features.PARA.start_year:features.PARA.end_year
                    load([features.PARA.out_folder  features.PARA.out_filename '_' num2str(valid_ind(i,1)) '_' num2str(y)  features.PARA.out_str '.mat'])
                    out_ii = [];
                    for j=1:size(features.TEMP.target_cells,1)
                        out_ii = cat(3, out_ii, mean(CG_out.T(features.TEMP.target_cells(j,1),:)));
                    end
                    out_i = cat(2, out_i, out_ii);
                end
                out = cat(1, out, out_i);
            end
            features.DATA.out = out;

            data_groups = out(1,:,:).*0 + 1;
            for i=1:size(data_groups,3)
                data_groups(1,:,i) = data_groups(1,:,i).*i;
            end
        end

    end
end

