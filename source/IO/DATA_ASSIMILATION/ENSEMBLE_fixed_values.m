classdef ENSEMBLE_fixed_values < ENSEMBLE_general


    properties
%         PARA
%         CONST
%         STATVAR
%         TEMP
    end
    
    methods
        
 function ensemble = provide_PARA(ensemble) 

            ensemble.PARA.variables1_names= [];
            ensemble.PARA.variables1_values= [];
            
            ensemble.PARA.variables2_names= [];
            ensemble.PARA.variables2_values= [];
            
            ensemble.PARA.variables3_names= [];
            ensemble.PARA.variables3_values= [];
            
%             ensemble.PARA.variables4_names = [];
%             ensemble.PARA.variables4_values= [];
%             
%             ensemble.PARA.variables5_names= [];
%             ensemble.PARA.variables5_values= [];
            
            ensemble.PARA.modify_class_name = [];
            ensemble.PARA.modify_class_index = [];
            ensemble.PARA.variable_in_ensemble = [];
            ensemble.PARA.variable_in_class = [];
            ensemble.PARA.position_in_class = [];
        end
        
        function ensemble = provide_CONST(ensemble)

        end
        
        function ensemble = provide_STATVAR(ensemble)

        end 
        
        
        function ensemble = finalize_init(ensemble, tile)
            ensemble.PARA.position_in_class(isnan(ensemble.PARA.position_in_class)) = 1; 
            
            index=[];
            if isempty(ensemble.PARA.variables2_names) || sum(isnan(ensemble.PARA.variables2_values(:)))>0
                for i=1:size(ensemble.PARA.variables1_values,1)
                    index=[index; [i]];
                end
            elseif isempty(ensemble.PARA.variables3_names) || sum(isnan(ensemble.PARA.variables3_values(:)))>0
                for i=1:size(ensemble.PARA.variables1_values,1)
                    for j=1:size(ensemble.PARA.variables2_values,1)
                        index=[index; [i j]];
                    end
                end
            else
                for i=1:size(ensemble.PARA.variables1_values,1)
                    for j=1:size(ensemble.PARA.variables2_values,1)
                        for k=1:size(ensemble.PARA.variables3_values,1)
                            index=[index; [i j k]];
                        end
                    end
                end
            end
            
            index = index(tile.PARA.worker_number,:);

            for i=1:size(ensemble.PARA.variables1_names,1)
                ensemble.STATVAR.(ensemble.PARA.variables1_names{i,1}) = ensemble.PARA.variables1_values(index(1,1), i);
            end
            if ~isempty(ensemble.PARA.variables2_names) && ~(sum(isnan(ensemble.PARA.variables2_values(:)))>0)
                for i=1:size(ensemble.PARA.variables2_names,1)
                    ensemble.STATVAR.(ensemble.PARA.variables2_names{i,1}) = ensemble.PARA.variables2_values(index(1,2), i);
                end
            end
            if ~isempty(ensemble.PARA.variables3_names) && ~(sum(isnan(ensemble.PARA.variables3_values(:)))>0)
                for i=1:size(ensemble.PARA.variables3_names,1)
                    ensemble.STATVAR.(ensemble.PARA.variables3_names{i,1}) = ensemble.PARA.variables3_values(index(1,3), i);
                end
            end
            
            ensemble = set_pointers2classes(ensemble, tile);
            variable_list = [];
            for i=1:size(ensemble.PARA.modify_class_pointer,1)
                variable_list = [variable_list; {class(ensemble.PARA.modify_class_pointer{i,1})}];
            end
            
            %write variables in PROVIDER
            for i=1:size(ensemble.PARA.modify_class_name,1)
                %rewrite in already called classes
                if any(strcmp(variable_list, ensemble.PARA.modify_class_name{i,1}))
                    pos = find(strcmp(variable_list,ensemble.PARA.modify_class_name{i,1}));
                    ensemble.PARA.modify_class_pointer{pos,1}.PARA.(ensemble.PARA.variable_in_class{i,1})(ensemble.PARA.position_in_class(i,1)) = ...
                    ensemble.STATVAR.(ensemble.PARA.variable_in_ensemble{i,1});
                end
               %rewrite in Provider
                    tile.RUN_INFO.PPROVIDER.CLASSES.(ensemble.PARA.modify_class_name{i,1}){ensemble.PARA.modify_class_index(i,1),1}.PARA.(ensemble.PARA.variable_in_class{i,1})(ensemble.PARA.position_in_class(i,1)) = ...
                    ensemble.STATVAR.(ensemble.PARA.variable_in_ensemble{i,1});
%                 
            end
            
           ensemble.PARA.modify_class_pointer = [];
                      
        end
        

    end
end

