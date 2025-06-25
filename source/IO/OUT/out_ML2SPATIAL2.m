%========================================================================
% CryoGrid OUT class out_ML2SPATIAL

% CryoGrid OUT class defining storage format of the output 
% S. Westermann, April 2025
%========================================================================


classdef out_ML2SPATIAL2 < matlab.mixin.Copyable
 

    properties
        MEAN
        STD
        TIMESTAMP
        TEMP
        PARA
        CONST
    end
    
    
    methods
        
        
        function out = provide_PARA(out)         
            out.PARA.tag = [];
        end

        
        function out = provide_CONST(out)

        end

        
        function out = provide_STATVAR(out)
    
        end

        
        function out = finalize_init(out, tile)

            out.MEAN.ML_out = [];
            out.STD.ML_out = [];
            if ~any(strcmp(fieldnames(tile.RUN_INFO.SPATIAL.STATVAR), 'ML_std'))
                tile.RUN_INFO.SPATIAL.STATVAR.ML_std = tile.RUN_INFO.SPATIAL.STATVAR.latitude .*0;
            end

        end
        
        %---------------time integration-------------

        function out = store_OUT(out, tile)

            %take mean/std and reshape, so that gridcell is dim1, time dim2 and the different variables dim3
            out.MEAN.ML_out = cat(2, out.MEAN.ML_out, reshape(mean(tile.STATVAR.ML_out,3), size(tile.STATVAR.ML_out,1), 1, size(tile.STATVAR.ML_out,2)));
            out.STD.ML_out = cat(2, out.STD.ML_out,  reshape(std(tile.STATVAR.ML_out, [], 3), size(tile.STATVAR.ML_out,1), 1, size(tile.STATVAR.ML_out,2)));
        end

        function out = write_OUT(out, tile)
            if ~(exist([tile.PARA.result_path tile.PARA.run_name])==7)
                mkdir([tile.PARA.result_path tile.PARA.run_name])
            end
            if tile.RUN_INFO.PARA.number_of_cores>1
                worker_no_str = ['_' num2str(tile.PARA.worker_number)];
            else
                worker_no_str = '';
            end
            if isempty(out.PARA.tag) || all(isnan(out.PARA.tag))
                save([tile.PARA.result_path tile.PARA.run_name '/' tile.PARA.run_name '_' num2str(tile.RUN_INFO.STATVAR.run_info_count) worker_no_str '.mat'], 'out')
            else
                save([tile.PARA.result_path tile.PARA.run_name '/' tile.PARA.run_name '_' num2str(tile.RUN_INFO.STATVAR.run_info_count) worker_no_str '_' out.PARA.tag '.mat'], 'out')
            end

            tile.RUN_INFO.SPATIAL.STATVAR.ML_std(tile.PARA.range,1) = mean(mean(out.STD.ML_out, 2),3);

            disp('writing out files')

            
            if tile.RUN_INFO.PARA.number_of_cores > 1
                spmdBarrier;
                for i=1:tile.RUN_INFO.PARA.number_of_cores
                    if i == tile.PARA.worker_number
                        %own worker, send info to all workers except yourself
                        for j=1:tile.RUN_INFO.PARA.number_of_cores
                            if j~=tile.PARA.worker_number
                                spmdSend(tile.PARA.range, j, 1);
                                spmdSend(mean(mean(out.STD.ML_out, 2),3), j, 2);
                            end
                        end
                    else
                        %receive from sending worker
                        range = spmdReceive(i, 1);
                        tile.RUN_INFO.SPATIAL.STATVAR.ML_std(range,1) = spmdReceive(i, 2);
                    end
                    spmdBarrier;
                end
                spmdBarrier;
            end

        end
      
        %-------------param file generation-----
    end
end