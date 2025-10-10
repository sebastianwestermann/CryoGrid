classdef mask_outfiles < matlab.mixin.Copyable
    
    properties
        PARA
        STATVAR
        CONST
        TEMP
    end

    methods
        function mask = provide_PARA(mask)
            mask.PARA.out_folder = [];
            mask.PARA.out_filename = [];
            mask.PARA.out_date = [];
            mask.PARA.additive = 0; 
        end


        function mask = provide_CONST(mask)

        end

        function mask = provide_STATVAR(mask)

        end

        function mask = finalize_init(mask, tile)
            str1 = num2str(mask.PARA.out_date(1,1)+100);
            str2 = num2str(mask.PARA.out_date(2,1)+100);

            valid_files = dir([mask.PARA.out_folder  mask.PARA.out_filename '_*' str1(2:3) str2(2:3) '.mat']);

            file_index = [];
            year_list=[];
            length_of_filename = size(mask.PARA.out_filename,2);
            for i=1:size(valid_files,1)
                next_pos = find(valid_files(i).name(length_of_filename+2:end) == '_', 1, 'first');
                file_index = [file_index; str2num(valid_files(i).name(length_of_filename+2:length_of_filename+2+next_pos-2))];
                year_list = [year_list; str2num(valid_files(i).name(length_of_filename+1+next_pos+1:length_of_filename+1+next_pos+4))];
            end

            mask.TEMP.file_index = sort(unique(file_index)); %list of available points.
        end

        function valid = apply_mask2(mask, tile, valid)
            valid_vector = zeros(size(tile.RUN_INFO.SPATIAL.STATVAR.latitude,1),1);
            %valid_vector(tile.RUN_INFO.CLUSTER.STATVAR.sample_centroid_index,1) = 1;
            for i=1:size(mask.TEMP.file_index,1)
                valid_vector(find(tile.RUN_INFO.SPATIAL.STATVAR.key(:,1)==mask.TEMP.file_index(i,1)),1) = 1;
            end
            valid_vector = logical(valid_vector);
            if mask.PARA.additive
                valid = valid | valid_vector;
            else
                valid = valid & valid_vector;
            end
            
        end

    end
end

