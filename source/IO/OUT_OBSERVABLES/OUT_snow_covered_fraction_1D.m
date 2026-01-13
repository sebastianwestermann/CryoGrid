classdef OUT_snow_covered_fraction_1D < matlab.mixin.Copyable
    
    properties
        PARA
        CONST
        TEMP
        STATVAR
        TIMESTAMP
        OUTPUT_TIME
    end
    
    methods

        function out = provide_PARA(out)
            out.PARA.timestamps = [];
            out.PARA.reference_SWE = []; 
            out.PARA.add_Xice = [];
            out.PARA.tag = [];

        end
        
        function out = provide_CONST(out)
            
        end
        
        function out = provide_STATVAR(out)
            
        end
        
        function out = finalize_init(out, tile)
            out.STATVAR.fsca = []; 
            %out.OUTPUT_TIME = tile.RUN_INFO.OPT.STATVAR.observation_times(find(tile.RUN_INFO.OPT.STATVAR.observation_times(:,1)-tile.PARA.start_time>=0, 1, 'first'), 1);
            out.OUTPUT_TIME = out.PARA.timestamps(find(out.PARA.timestamps(:,1)-tile.PARA.start_time>=0, 1, 'first'), 1);

            
        end
        
        function out = store_OUT(out, tile)

            if tile.t >= out.OUTPUT_TIME

                out.TIMESTAMP = [out.TIMESTAMP tile.t];
                % out.OUTPUT_TIME = tile.RUN_INFO.OPT.STATVAR.observation_times(find(tile.RUN_INFO.OPT.STATVAR.observation_times(:,1)-tile.t>0, 1, 'first'), 1);
                out.OUTPUT_TIME = out.PARA.timestamps(find(out.PARA.timestamps(:,1)-tile.t>0, 1, 'first'), 1);

                result = 0;
                CURRENT = tile.TOP.NEXT;
                while ~(strcmp(class(CURRENT), 'Bottom'))
                    class_name = class(CURRENT);
                    if strcmp(class_name(1:4), 'SNOW')
                        result = 1;
                        break
                    end
                    if out.PARA.add_Xice == 1
                        if any(strcmp(properties(CURRENT), 'CHILD'))
                            if isempty(CURRENT.CHILD) || CURRENT.CHILD == 0
                                ice_height = CURRENT.STATVAR.Xice(1,1) ./ CURRENT.STATVAR.area(1,1);
                                result = min(1,max(0, ice_height./ out.PARA.reference_SWE));
                            else
                                ice_height = (CURRENT.CHILD.STATVAR.ice + CURRENT.STATVAR.Xice(1,1)) ./ CURRENT.STATVAR.area(1,1);
                                result = min(1,max(0, ice_height./ out.PARA.reference_SWE));
                            end
                            break
                        end
                    else
                        if any(strcmp(properties(CURRENT), 'CHILD'))
                            ice_height = CURRENT.CHILD.STATVAR.ice ./ CURRENT.STATVAR.area(1,1);
                            result = min(1,max(0, ice_height./ out.PARA.reference_SWE));
                        end
                        break
                    end
                    CURRENT = CURRENT.NEXT;
                end
                
                out.STATVAR.fsca = [out.STATVAR.fsca; result];
            end
            

            %add possibility to save output 

        end

        function result = move_out2obs(out)
            result = out.STATVAR.fsca;
        end
            
        
        function out = reset_new_stratigraphy(out, tile)
            
        end
        
    end
end

