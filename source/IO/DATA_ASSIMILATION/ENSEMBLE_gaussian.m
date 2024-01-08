classdef ENSEMBLE_gaussian < matlab.mixin.Copyable


    properties
        PARA
        CONST
        STATVAR
        TEMP
    end
    
    methods
        
 function ensemble = provide_PARA(ensemble) 

            ensemble.PARA.variable_in_ensemble = [];
            ensemble.PARA.center = [];
            ensemble.PARA.width = [];
            ensemble.PARA.lower_bound = [];
            ensemble.PARA.upper_bound = [];
            
            ensemble.PARA.modify_class_name = [];
            ensemble.PARA.modify_class_index = [];
            ensemble.PARA.variable_in_class = [];
            ensemble.PARA.index_V = [];
            ensemble.PARA.index_H = [];
        end
        
        function ensemble = provide_CONST(ensemble)

        end
        
        function ensemble = provide_STATVAR(ensemble)

        end 
        
        
        function ensemble = finalize_init(ensemble, tile)
            if isempty(ensemble.PARA.index_V) || (size(ensemble.PARA.index_V,1)==1 && isnan(ensemble.PARA.index_V))
                ensemble.PARA.index_V = ones(size(ensemble.PARA.modify_class_index,1));
            end
            ensemble.PARA.index_V(isnan(ensemble.PARA.index_V)) = 1;
            if isempty(ensemble.PARA.index_H) || (size(ensemble.PARA.index_H,1)==1 && isnan(ensemble.PARA.index_H))
                ensemble.PARA.index_H = ones(size(ensemble.PARA.modify_class_index,1));
            end
            ensemble.PARA.index_H(isnan(ensemble.PARA.index_H)) = 1;
            
            
            %variables in Gaussian space, use this for resampling
            ensemble.TEMP.value_gaussian = [];
            ensemble.TEMP.mean_gaussian = [];
            ensemble.TEMP.std_gaussian = [];
            ensemble.TEMP.variable_name = {};
            
            for i=1:size(ensemble.PARA.variable_in_ensemble,1)
                ensemble.STATVAR.(ensemble.PARA.variable_in_ensemble{i,1}) = max(ensemble.PARA.lower_bound(i,1), min(ensemble.PARA.upper_bound(i,1), ensemble.PARA.center(i,1) + randn(1,tile.PARA.ensemble_size) .* ensemble.PARA.width(i,1)));
                ensemble.TEMP.value_gaussian = [ensemble.TEMP.value_gaussian; ensemble.STATVAR.(ensemble.PARA.variable_in_ensemble{i,1})(1,tile.PARA.worker_number)];
                ensemble.TEMP.mean_gaussian = [ensemble.TEMP.mean_gaussian; ensemble.PARA.center(i,1)];
                ensemble.TEMP.std_gaussian = [ensemble.TEMP.std_gaussian; ensemble.PARA.width(i,1)];
                ensemble.TEMP.variable_name = [ensemble.TEMP.variable_name; ensemble.PARA.variable_in_ensemble{i,1}];
                
                ensemble.STATVAR.(ensemble.PARA.variable_in_ensemble{i,1}) = ensemble.STATVAR.(ensemble.PARA.variable_in_ensemble{i,1})(1,tile.PARA.worker_number);
            end
            
            %write variables in PROVIDER
            for i=1:size(ensemble.PARA.modify_class_name,1)
                tile.RUN_INFO.PPROVIDER.CLASSES.(ensemble.PARA.modify_class_name{i,1}){ensemble.PARA.modify_class_index(i,1),1}.PARA.(ensemble.PARA.variable_in_class{i,1})(ensemble.PARA.index_V(i,1), ensemble.PARA.index_H(i,1)) = ...
                    ensemble.STATVAR.(ensemble.PARA.variable_in_ensemble{i,1});
            end
                      
        end
        
        %add possibility to restrict the variables only to the ones which should be changed       
        function ensemble = recalculate_ensemble_parameters_after_DA(ensemble, tile, variable_list)
            %recalculate ensemble parameters based on DA weights and
            %previous ensemble-the new variable centers and widths are set
            %by the DA class, plus the variables that are affected by the
            %DA

            for i=1:size(variable_list,1)
                pos = find(strcmp(variable_list{i,1}, ensemble.TEMP.variable_name));
                ensemble.STATVAR.(variable_list{i,1}) = ensemble.TEMP.value_gaussian(pos,1);

            end
            
            %establish pointers to all classes in the stratigraphy
            ensemble = set_pointers2classes(ensemble, tile);
            
            %rewrite PARA in PROVIDER and all classes in the
            %stratigraphy
            for i=1:size(ensemble.PARA.modify_class_name,1)
                if any(strcmp(variable_list, ensemble.PARA.variable_in_ensemble{i,1}))
                    tile.RUN_INFO.PPROVIDER.CLASSES.(ensemble.PARA.modify_class_name{i,1}){ensemble.PARA.modify_class_index(i,1),1}.PARA.(ensemble.PARA.variable_in_class{i,1})(ensemble.PARA.index_V(i,1), ensemble.PARA.index_H(i,1)) = ...
                        max(ensemble.PARA.lower_bound(i,1), min(ensemble.PARA.upper_bound(i,1), ensemble.STATVAR.(ensemble.PARA.variable_in_ensemble{i,1})));
                    for j=1:size(ensemble.PARA.modify_class_pointer{i,1})
                        ensemble.PARA.modify_class_pointer{i,1}(j,1).PARA.(ensemble.PARA.variable_in_class{i,1})(ensemble.PARA.index_V(i,1), ensemble.PARA.index_H(i,1)) = ...
                            max(ensemble.PARA.lower_bound(i,1), min(ensemble.PARA.upper_bound(i,1), ensemble.STATVAR.(ensemble.PARA.variable_in_ensemble{i,1})));
                    end
                end
            end
            
            ensemble.PARA.modify_class_pointer = [];
        end
        
        function ensemble = set_pointers2classes(ensemble, tile)
            for i=1:size(ensemble.PARA.modify_class_name,1)
                ensemble.PARA.modify_class_pointer{i,1} = {};
            end
            
            ensemble = find_classes_in_variable(ensemble, tile, 1);
            
        end
        
        function ensemble = find_classes_in_variable(ensemble, current_class, level)
            class_name = class(current_class);
            if isobject(current_class)
                for i=1:size(ensemble.PARA.modify_class_name, 1)
                    if strcmp(ensemble.PARA.modify_class_name{i,1}, class_name)
                        if  ensemble.PARA.modify_class_index(i,1) == current_class.PARA.class_index
                            dec = 1;
                            for j=1:size(ensemble.PARA.modify_class_pointer{i,1}, 1)
                                if isequal(ensemble.PARA.modify_class_pointer{i,1}(j,1), current_class)
                                    dec = 0;
                                end
                            end
                            if dec
                                ensemble.PARA.modify_class_pointer{i,1} = [ensemble.PARA.modify_class_pointer{i,1}; current_class];
                            end
                        end
                    end
                end
            end
            if isstruct(current_class) || isobject(current_class)
                variables = fieldnames(current_class);
                for i = 1:size(variables,1)
                    if iscell(current_class.(variables{i,1}))
                        if level <= 5
                            for j=1:size(current_class.(variables{i,1}), 1)
                                for k=1:size(current_class.(variables{i,1}), 2)
                                    ensemble = find_classes_in_variable(ensemble, current_class.(variables{i,1}){j,k}, level+1);
                                end
                            end
                        end
                    elseif isstruct(current_class.(variables{i,1}))
                        if level <= 5
                            ensemble = find_classes_in_variable(ensemble, current_class.(variables{i,1}), level+1);
                        end
                    elseif isobject(current_class.(variables{i,1})) && ~strcmp(variables{i,1}, 'RUN_INFO')
                        if level <= 5
                            ensemble = find_classes_in_variable(ensemble, current_class.(variables{i,1}), level+1);
                        end
                    end
                end
            end
        end

    end
end

