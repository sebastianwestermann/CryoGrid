%========================================================================
% CryoGrid OUT class out_ML2SPATIAL

% CryoGrid OUT class defining storage format of the output 
% S. Westermann, April 2025
%========================================================================


classdef out_ML2SPATIAL < matlab.mixin.Copyable
 

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

           for i=1:size(tile.PARA.variables,1)
               out.MEAN.(tile.PARA.variables{i,1}) = [];
               out.STD.(tile.PARA.variables{i,1}) = [];
               if ~any(strcmp(fieldnames(tile.RUN_INFO.SPATIAL.STATVAR), ['ml_std_' (tile.PARA.variables{i,1})]))
                   tile.RUN_INFO.SPATIAL.STATVAR.(['ml_std_' (tile.PARA.variables{i,1})]) = tile.RUN_INFO.SPATIAL.STATVAR.latitude .*0;
               end
           end



        end
        
        %---------------time integration-------------
                    
        function out = store_OUT(out, tile, i)           
            
            out.MEAN.(tile.PARA.variables{i,1}) = [out.MEAN.(tile.PARA.variables{i,1}) mean(tile.STATVAR.(tile.PARA.variables{i,1}),2)];
            out.STD.(tile.PARA.variables{i,1}) = [out.STD.(tile.PARA.variables{i,1}) std(tile.STATVAR.(tile.PARA.variables{i,1}),[],2)];

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
                save([tile.PARA.result_path tile.PARA.run_name '/' tile.PARA.run_name worker_no_str '.mat'], 'out')
            else
                save([tile.PARA.result_path tile.PARA.run_name '/' tile.PARA.run_name worker_no_str '_' out.PARA.tag '.mat'], 'out')
            end

            for i=1:size(tile.PARA.variables,1)
                tile.RUN_INFO.SPATIAL.STATVAR.(['ml_std_' (tile.PARA.variables{i,1})])(tile.PARA.range,1) = mean(out.STD.(tile.PARA.variables{i,1}),2);
            end

            disp('writing out files')

            
            if tile.RUN_INFO.PARA.number_of_cores > 1
                spmdBarrier;
                for i=1:tile.RUN_INFO.PARA.number_of_cores
                    if i == tile.PARA.worker_number
                        %own worker, send info to all workers except yourself
                        for j=1:tile.RUN_INFO.PARA.number_of_cores
                            if j~=tile.PARA.worker_number
                                spmdSend(tile.PARA.range, j, 1);
                                for k=1:size(tile.PARA.variables,1)
                                    spmdSend(mean(out.STD.(tile.PARA.variables{k,1}),2), j, k+1);
                                end
                            end
                        end
                    else
                        %receive from sending worker 
                        range = spmdReceive(i, 1);
                        for k=1:size(tile.PARA.variables,1)
                            tile.RUN_INFO.SPATIAL.STATVAR.(['ml_std_' (tile.PARA.variables{k,1})])(range,1) = spmdReceive(i, k+1);
                        end
                    end
                    spmdBarrier;
                end
                spmdBarrier;
            end

        end
      
        %-------------param file generation-----
    end
end