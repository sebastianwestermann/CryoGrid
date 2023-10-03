classdef ENSEMBLE_general < matlab.mixin.Copyable


    properties
        PARA
        CONST
        STATVAR
        TEMP
    end
    
    methods
        
 function ensemble = provide_PARA(ensemble) 

            ensemble.PARA.gaussian_variables_name= [];
            ensemble.PARA.gaussian_variables_center= [];
            ensemble.PARA.gaussian_variables_width= [];
            
            ensemble.PARA.boxcar_variables_name = [];
            ensemble.PARA.boxcar_variables_lower_bound = [];
            ensemble.PARA.boxcar_variables_upper_bound = [];
            
            ensemble.PARA.lognormal_variables_name = [];
            ensemble.PARA.lognormal_variables_lower_bound = [];
            ensemble.PARA.lognormal_variables_std = [];
            ensemble.PARA.lognormal_variables_mean = [];
            
            ensemble.PARA.modify_class_name = [];
            ensemble.PARA.modify_class_index = [];
            ensemble.PARA.variable_in_ensemble = [];
            ensemble.PARA.variable_in_class = [];
        end
        
        function ensemble = provide_CONST(ensemble)

        end
        
        function ensemble = provide_STATVAR(ensemble)

        end 
        
        
        function ensemble = finalize_init(ensemble, tile)
            
            %variables in Gaussian space, use this for resampling
            ensemble.TEMP.value_gaussian = [];
            ensemble.TEMP.mean_gaussian = [];
            ensemble.TEMP.std_gaussian = [];
            ensemble.TEMP.variable_type = []; %1 gaussian; 2: boxcar; 3: lognormal
            ensemble.TEMP.variable_name = {};
            
            for i=1:size(ensemble.PARA.gaussian_variables_name,1)
                ensemble.STATVAR.(ensemble.PARA.gaussian_variables_name{i,1}) = ensemble.PARA.gaussian_variables_center(i,1) + randn(1,tile.PARA.ensemble_size) .* ensemble.PARA.gaussian_variables_width(i,1);
                ensemble.TEMP.value_gaussian = [ensemble.TEMP.value_gaussian; ensemble.STATVAR.(ensemble.PARA.gaussian_variables_name{i,1})(1,tile.PARA.worker_number)];
                ensemble.TEMP.mean_gaussian = [ensemble.TEMP.mean_gaussian; ensemble.PARA.gaussian_variables_center(i,1)];
                ensemble.TEMP.std_gaussian = [ensemble.TEMP.std_gaussian; ensemble.PARA.gaussian_variables_width(i,1)];
                ensemble.TEMP.variable_type = [ensemble.TEMP.variable_type; 1]; %1 gaussian; 2: boxcar; 3: lognormal
                ensemble.TEMP.variable_name = [ensemble.TEMP.variable_name; ensemble.PARA.gaussian_variables_name{i,1}];
                
                ensemble.STATVAR.(ensemble.PARA.gaussian_variables_name{i,1}) = ensemble.STATVAR.(ensemble.PARA.gaussian_variables_name{i,1})(1,tile.PARA.worker_number);
            end
            for i=1:size(ensemble.PARA.boxcar_variables_name,1)
                %logit transformation
                prior_median = (ensemble.PARA.boxcar_variables_lower_bound(i,1) +  ensemble.PARA.boxcar_variables_upper_bound(i,1))./2;
                prior_mean = logit(ensemble, prior_median, ensemble.PARA.boxcar_variables_lower_bound(i,1), ensemble.PARA.boxcar_variables_upper_bound(i,1));
                prior_std=1.7; %makes a near perfect boxcar 
                prior_untransformed = prior_mean + prior_std .* randn(1,tile.PARA.ensemble_size);
                
                ensemble.TEMP.value_gaussian = [ensemble.TEMP.value_gaussian; prior_untransformed(1,tile.PARA.worker_number)];
                ensemble.TEMP.mean_gaussian = [ensemble.TEMP.mean_gaussian; prior_mean];
                ensemble.TEMP.std_gaussian = [ensemble.TEMP.std_gaussian; prior_std];
                ensemble.TEMP.variable_type = [ensemble.TEMP.variable_type; 2]; %1 gaussian; 2: boxcar; 3: lognormal
                ensemble.TEMP.variable_name = [ensemble.TEMP.variable_name; ensemble.PARA.boxcar_variables_name{i,1}];
                
                ensemble.STATVAR.(ensemble.PARA.boxcar_variables_name{i,1}) = expit(ensemble, prior_untransformed, ensemble.PARA.boxcar_variables_lower_bound(i,1), ensemble.PARA.boxcar_variables_upper_bound(i,1));
                ensemble.STATVAR.(ensemble.PARA.boxcar_variables_name{i,1}) = ensemble.STATVAR.(ensemble.PARA.boxcar_variables_name{i,1})(1,tile.PARA.worker_number);
            end

            for i=1:size(ensemble.PARA.lognormal_variables_name,1)
                mean_gaussian = log(ensemble.PARA.lognormal_variables_mean(i,1).^2./sqrt(ensemble.PARA.lognormal_variables_mean(i,1).^2 + ensemble.PARA.lognormal_variables_std(i,1).^2));
                std_gaussian = sqrt(log(1 + ensemble.PARA.lognormal_variables_std(i,1).^2 ./ ensemble.PARA.lognormal_variables_mean(i,1).^2)); %tranform mean and std from real to Gaussian space
                value_gaussian = mean_gaussian + std_gaussian .* randn(1,tile.PARA.ensemble_size);

                ensemble.TEMP.value_gaussian = [ensemble.TEMP.value_gaussian; value_gaussian(1,tile.PARA.worker_number)];
                ensemble.TEMP.mean_gaussian = [ensemble.TEMP.mean_gaussian; mean_gaussian];
                ensemble.TEMP.std_gaussian = [ensemble.TEMP.std_gaussian; std_gaussian];
                ensemble.TEMP.variable_type = [ensemble.TEMP.variable_type; 3]; %1 gaussian; 2: boxcar; 3: lognormal
                ensemble.TEMP.variable_name = [ensemble.TEMP.variable_name; ensemble.PARA.lognormal_variables_name{i,1}];
                
                ensemble.STATVAR.(ensemble.PARA.lognormal_variables_name{i,1}) = ensemble.PARA.lognormal_variables_lower_bound(i,1) + exp(value_gaussian);
                ensemble.STATVAR.(ensemble.PARA.lognormal_variables_name{i,1}) = ensemble.STATVAR.(ensemble.PARA.lognormal_variables_name{i,1})(1,tile.PARA.worker_number);
            end
            
            %write variables in PROVIDER
            for i=1:size(ensemble.PARA.modify_class_name,1)
                tile.RUN_INFO.PPROVIDER.CLASSES.(ensemble.PARA.modify_class_name{i,1}){ensemble.PARA.modify_class_index(i,1),1}.PARA.(ensemble.PARA.variable_in_class{i,1}) = ...
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
                if ensemble.TEMP.variable_type(pos,1) == 1
                    ensemble.STATVAR.(variable_list{i,1}) = ensemble.TEMP.value_gaussian(pos,1);
                elseif ensemble.TEMP.variable_type(pos,1) == 2
                    pos_boxcar = find(strcmp(variable_list{i,1}, ensemble.PARA.boxcar_variables_name)); 
                    ensemble.STATVAR.(variable_list{i,1}) = expit(ensemble, ensemble.TEMP.value_gaussian(pos,1), ensemble.PARA.boxcar_variables_lower_bound(pos_boxcar,1), ensemble.PARA.boxcar_variables_upper_bound(pos_boxcar,1));
                elseif ensemble.TEMP.variable_type(pos,1) == 3
                    pos_lognormal = find(strcmp(variable_list{i,1}, ensemble.PARA.lognormal_variables_name));
                    ensemble.STATVAR.(variable_list{i,1}) = ensemble.PARA.lognormal_variables_lower_bound(pos_lognormal,1) + exp(ensemble.TEMP.value_gaussian(pos,1));
                end
            end
            
            %establish pointers to all classes in the stratigraphy
            ensemble = set_pointers2classes(ensemble, tile);
            
            %rewrite PARA in PROVIDER and all classes in the
            %stratigraphy
            for i=1:size(ensemble.PARA.modify_class_name,1)
                if any(strcmp(variable_list, ensemble.PARA.variable_in_ensemble{i,1}))
                    tile.RUN_INFO.PPROVIDER.CLASSES.(ensemble.PARA.modify_class_name{i,1}){ensemble.PARA.modify_class_index(i,1),1}.PARA.(ensemble.PARA.variable_in_class{i,1}) = ...
                        ensemble.STATVAR.(ensemble.PARA.variable_in_ensemble{i,1});
                    for j=1:size(ensemble.PARA.modify_class_pointer{i,1})
                        ensemble.PARA.modify_class_pointer{i,1}(j,1).PARA.(ensemble.PARA.variable_in_class{i,1}) = ensemble.STATVAR.(ensemble.PARA.variable_in_ensemble{i,1});
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
        
        
        
        
        %logit tranformation for boxcar variables
        function res = logit(ensemble, x, a, b) 
            res = log(((x-a)./(b-a))./(1-(x-a)./(b-a))); % Transforms from physical to unbounded (Gaussian) space
        end
        
        function res = expit(ensemble, xt,a,b) 
            res = a+(b-a)./(1+exp(-xt)); % Transforms back from unbounded to physical space
        end
    end
end

