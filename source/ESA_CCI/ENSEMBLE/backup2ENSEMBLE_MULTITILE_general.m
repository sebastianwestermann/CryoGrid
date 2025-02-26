classdef ENSEMBLE_MULTITILE_general < matlab.mixin.Copyable


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
            
            variables = {'gaussian_variables_name'; 'gaussian_variables_center'; 'gaussian_variables_width'; 'boxcar_variables_name'; 'boxcar_variables_lower_bound';...
                'boxcar_variables_upper_bound'; 'lognormal_variables_name'; 'lognormal_variables_lower_bound'; 'lognormal_variables_std'; 'lognormal_variables_mean'};
            
            for i=1:size(variables,1)
                if size(ensemble.PARA.(variables{i,1}),2) == 1
                    ensemble.PARA.(variables{i,1}) = repmat(ensemble.PARA.(variables{i,1}),1,tile.PARA.number_of_realizations);
                    %ensemble.PARA.(variables{i,1}) = repmat(ensemble.PARA.(variables{i,1}),1,tile.PARA.number_of_realizations);
                    %replaced by ensemble_size at all occurrences
                end
            end
           
            %variables in Gaussian space, use this for resampling
            ensemble.TEMP.value_gaussian = [];
            ensemble.TEMP.mean_gaussian = [];
            ensemble.TEMP.std_gaussian = [];
            ensemble.TEMP.variable_type = []; %1 gaussian; 2: boxcar; 3: lognormal
            ensemble.TEMP.variable_name = {};
            
            for i=1:size(ensemble.PARA.gaussian_variables_name,1)
                ensemble.STATVAR.(ensemble.PARA.gaussian_variables_name{i,1}) = [];
                for j=1:size(ensemble.PARA.gaussian_variables_center,2)
                    ensemble.STATVAR.(ensemble.PARA.gaussian_variables_name{i,1}) = [ensemble.STATVAR.(ensemble.PARA.gaussian_variables_name{i,1}); ensemble.PARA.gaussian_variables_center(i,j) + randn(1,tile.PARA.number_of_realizations) .* ensemble.PARA.gaussian_variables_width(i,j)];
                end
                ensemble.STATVAR.(ensemble.PARA.gaussian_variables_name{i,1}) = ensemble.STATVAR.(ensemble.PARA.gaussian_variables_name{i,1})(:)';
                ensemble.TEMP.value_gaussian = [ensemble.TEMP.value_gaussian; ensemble.STATVAR.(ensemble.PARA.gaussian_variables_name{i,1})];
                ensemble.TEMP.mean_gaussian = [ensemble.TEMP.mean_gaussian; ensemble.PARA.gaussian_variables_center(i,:)];
                ensemble.TEMP.std_gaussian = [ensemble.TEMP.std_gaussian; ensemble.PARA.gaussian_variables_width(i,:)];
                ensemble.TEMP.variable_type = [ensemble.TEMP.variable_type; 1]; %1 gaussian; 2: boxcar; 3: lognormal
                ensemble.TEMP.variable_name = [ensemble.TEMP.variable_name; ensemble.PARA.gaussian_variables_name{i,1}];
%                 ensemble.STATVAR.(ensemble.PARA.gaussian_variables_name{i,1}) = ensemble.STATVAR.(ensemble.PARA.gaussian_variables_name{i,1})(1,tile.PARA.worker_number);
            end
            for i=1:size(ensemble.PARA.boxcar_variables_name,1)
                %logit transformation
                prior_median = (ensemble.PARA.boxcar_variables_lower_bound(i,1) +  ensemble.PARA.boxcar_variables_upper_bound(i,1))./2;
                prior_mean = logit(ensemble, prior_median, ensemble.PARA.boxcar_variables_lower_bound(i,1), ensemble.PARA.boxcar_variables_upper_bound(i,1));
                prior_std = repmat(1.7, 1, tile.PARA.number_of_realizations); %makes a near perfect boxcar 
                prior_untransformed = [];
                prior_transformed = [];
                for j=1:size(ensemble.PARA.boxcar_variables_lower_bound,2)
                    random_gaussian = prior_mean(1,j) + prior_std(1,j) .* randn(1,tile.PARA.number_of_realizations);
                    prior_untransformed = [prior_untransformed; random_gaussian];
                    prior_transformed = [prior_transformed; expit(ensemble, prior_untransformed, ensemble.PARA.boxcar_variables_lower_bound(i,j), ensemble.PARA.boxcar_variables_upper_bound(i,j))];
                end
                prior_untransformed = prior_untransformed(:)';
                prior_transformed = prior_transformed(:)';
                ensemble.TEMP.value_gaussian = [ensemble.TEMP.value_gaussian; prior_untransformed];
                ensemble.TEMP.mean_gaussian = [ensemble.TEMP.mean_gaussian; prior_mean];
                ensemble.TEMP.std_gaussian = [ensemble.TEMP.std_gaussian; prior_std];
                ensemble.TEMP.variable_type = [ensemble.TEMP.variable_type; 2]; %1 gaussian; 2: boxcar; 3: lognormal
                ensemble.TEMP.variable_name = [ensemble.TEMP.variable_name; ensemble.PARA.boxcar_variables_name{i,1}];
                ensemble.STATVAR.(ensemble.PARA.boxcar_variables_name{i,1}) = prior_transformed; %expit(ensemble, prior_untransformed, ensemble.PARA.boxcar_variables_lower_bound(i,1), ensemble.PARA.boxcar_variables_upper_bound(i,1));
%                 ensemble.STATVAR.(ensemble.PARA.boxcar_variables_name{i,1}) = ensemble.STATVAR.(ensemble.PARA.boxcar_variables_name{i,1})(1,tile.PARA.worker_number);
            end

            for i=1:size(ensemble.PARA.lognormal_variables_name,1)
                mean_gaussian = log(ensemble.PARA.lognormal_variables_mean(i,1).^2./sqrt(ensemble.PARA.lognormal_variables_mean(i,1).^2 + ensemble.PARA.lognormal_variables_std(i,1).^2));
                std_gaussian = sqrt(log(1 + ensemble.PARA.lognormal_variables_std(i,1).^2 ./ ensemble.PARA.lognormal_variables_mean(i,1).^2)); %tranform mean and std from real to Gaussian space
                value_gaussian = [];
                value_normal_space = [];
                for j=1:size(ensemble.PARA.lognormal_variables_mean,2)
                    random_gaussian = mean_gaussian(1,j) + std_gaussian(1,j) .* randn(1,tile.PARA.number_of_realizations);
                    value_gaussian = [value_gaussian; random_gaussian];
                    value_normal_space = [value_normal_space; ensemble.PARA.lognormal_variables_lower_bound(i,j) + exp(random_gaussian)];
                end
                value_gaussian = value_gaussian(:)';
                value_normal_space = value_normal_space(:)';
                ensemble.TEMP.value_gaussian = [ensemble.TEMP.value_gaussian; value_gaussian];
                ensemble.TEMP.mean_gaussian = [ensemble.TEMP.mean_gaussian; mean_gaussian];
                ensemble.TEMP.std_gaussian = [ensemble.TEMP.std_gaussian; std_gaussian];
                ensemble.TEMP.variable_type = [ensemble.TEMP.variable_type; 3]; %1 gaussian; 2: boxcar; 3: lognormal
                ensemble.TEMP.variable_name = [ensemble.TEMP.variable_name; ensemble.PARA.lognormal_variables_name{i,1}];
                
                ensemble.STATVAR.(ensemble.PARA.lognormal_variables_name{i,1}) = value_normal_space; %ensemble.PARA.lognormal_variables_lower_bound(i,1) + exp(value_gaussian);
%                 ensemble.STATVAR.(ensemble.PARA.lognormal_variables_name{i,1}) = ensemble.STATVAR.(ensemble.PARA.lognormal_variables_name{i,1})(1,tile.PARA.worker_number);
            end
            
            
            %write variables in PROVIDER
            for i=1:size(ensemble.PARA.modify_class_name,1)
                tile.RUN_INFO.PPROVIDER.CLASSES.(ensemble.PARA.modify_class_name{i,1}){ensemble.PARA.modify_class_index(i,1),1}.PARA.(ensemble.PARA.variable_in_class{i,1}) = ...
                    ensemble.STATVAR.(ensemble.PARA.variable_in_ensemble{i,1});
                tile.RUN_INFO.PPROVIDER.CLASSES.(ensemble.PARA.modify_class_name{i,1}){ensemble.PARA.modify_class_index(i,1),1}.PARA.number_of_realizations = tile.PARA.number_of_realizations; 
            end 

%            for i=1:size(ensemble.PARA.gaussian_variables_name,1)
%                ensemble.STATVAR.(ensemble.PARA.gaussian_variables_name{i,1}) = ensemble.PARA.gaussian_variables_center(i,1) + randn(1,tile.PARA.number_of_realizations) .* ensemble.PARA.gaussian_variables_width(i,1);
%            end
%            for i=1:size(ensemble.PARA.boxcar_variables_name,1)
%                ensemble.STATVAR.(ensemble.PARA.boxcar_variables_name{i,1}) = ensemble.PARA.boxcar_variables_center(i,1) + (rand(1,tile.PARA.number_of_realizations)-0.5) .* ensemble.PARA.boxcar_variables_width(i,1);
%            end

           % ensemble.STATVAR = [];
            
            
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
                    ensemble.STATVAR.(variable_list{i,1}) = ensemble.TEMP.value_gaussian(pos,:);
                elseif ensemble.TEMP.variable_type(pos,1) == 2
                    pos_boxcar = find(strcmp(variable_list{i,1}, ensemble.PARA.boxcar_variables_name)); 
                    ensemble.STATVAR.(variable_list{i,1}) = expit(ensemble, ensemble.TEMP.value_gaussian(pos,:), ensemble.PARA.boxcar_variables_lower_bound(pos_boxcar,1), ensemble.PARA.boxcar_variables_upper_bound(pos_boxcar,1));
                elseif ensemble.TEMP.variable_type(pos,1) == 3
                    pos_lognormal = find(strcmp(variable_list{i,1}, ensemble.PARA.lognormal_variables_name));
                    ensemble.STATVAR.(variable_list{i,1}) = ensemble.PARA.lognormal_variables_lower_bound(pos_lognormal,1) + exp(ensemble.TEMP.value_gaussian(pos,:));
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

        
        
        function ensemble = finalize_init2(ensemble, tile)
            
            variables = {'gaussian_variables_name'; 'gaussian_variables_center'; 'gaussian_variables_width'; 'boxcar_variables_name'; 'boxcar_variables_lower_bound';...
                'boxcar_variables_upper_bound'; 'lognormal_variables_name'; 'lognormal_variables_lower_bound'; 'lognormal_variables_std'; 'lognormal_variables_mean'};
            
            for i=1:size(variables,1)
                if size(ensemble.PARA.(variables{i,1}),2) == 1
                    ensemble.PARA.(variables{i,1}) = repmat(ensemble.PARA.(variables{i,1}),1,tile.PARA.number_of_realizations);
                end
            end
           
            %variables in Gaussian space, use this for resampling
            ensemble.TEMP.value_gaussian = [];
            ensemble.TEMP.mean_gaussian = [];
            ensemble.TEMP.std_gaussian = [];
            ensemble.TEMP.variable_type = []; %1 gaussian; 2: boxcar; 3: lognormal
            ensemble.TEMP.variable_name = {};
            
            for i=1:size(ensemble.PARA.gaussian_variables_name,1)
                ensemble.STATVAR.(ensemble.PARA.gaussian_variables_name{i,1}) = [];
                for j=1:size(ensemble.PARA.gaussian_variables_center,2)
                    ensemble.STATVAR.(ensemble.PARA.gaussian_variables_name{i,1}) = [ensemble.STATVAR.(ensemble.PARA.gaussian_variables_name{i,1}); ensemble.PARA.gaussian_variables_center(i,1) + randn(1,tile.PARA.ensemble_size) .* ensemble.PARA.gaussian_variables_width(i,1)];
                end
                ensemble.STATVAR.(ensemble.PARA.gaussian_variables_name{i,1}) = ensemble.STATVAR.(ensemble.PARA.gaussian_variables_name{i,1})(:)';
                ensemble.TEMP.value_gaussian = [ensemble.TEMP.value_gaussian; ensemble.STATVAR.(ensemble.PARA.gaussian_variables_name{i,1})];
                ensemble.TEMP.mean_gaussian = [ensemble.TEMP.mean_gaussian; ensemble.PARA.gaussian_variables_center(i,:)];
                ensemble.TEMP.std_gaussian = [ensemble.TEMP.std_gaussian; ensemble.PARA.gaussian_variables_width(i,:)];
                ensemble.TEMP.variable_type = [ensemble.TEMP.variable_type; 1]; %1 gaussian; 2: boxcar; 3: lognormal
                ensemble.TEMP.variable_name = [ensemble.TEMP.variable_name; ensemble.PARA.gaussian_variables_name{i,1}];
%                 ensemble.STATVAR.(ensemble.PARA.gaussian_variables_name{i,1}) = ensemble.STATVAR.(ensemble.PARA.gaussian_variables_name{i,1})(1,tile.PARA.worker_number);
            end
            for i=1:size(ensemble.PARA.boxcar_variables_name,1)
                %logit transformation
                prior_median = (ensemble.PARA.boxcar_variables_lower_bound(i,:) +  ensemble.PARA.boxcar_variables_upper_bound(i,:))./2;
                prior_mean = logit(ensemble, prior_median, ensemble.PARA.boxcar_variables_lower_bound(i,:), ensemble.PARA.boxcar_variables_upper_bound(i,:));
                prior_std = repmat(1.7, 1, tile.PARA.ensemble_size); %makes a near perfect boxcar 
                prior_untransformed = [];
                prior_transformed = [];
                for j=1:size(ensemble.PARA.boxcar_variables_lower_bound,2)
                    random_gaussian = prior_mean(1,j) + prior_std(1,j) .* randn(1,tile.PARA.ensemble_size);
                    prior_untransformed = [prior_untransformed; random_gaussian];
                    prior_transformed = [prior_transformed; expit(ensemble, prior_untransformed, ensemble.PARA.boxcar_variables_lower_bound(i,j), ensemble.PARA.boxcar_variables_upper_bound(i,j))];
                end
                prior_untransformed = prior_untransformed(:)';
                prior_transformed = prior_transformed(:)';
                ensemble.TEMP.value_gaussian = [ensemble.TEMP.value_gaussian; prior_untransformed];
                ensemble.TEMP.mean_gaussian = [ensemble.TEMP.mean_gaussian; prior_mean];
                ensemble.TEMP.std_gaussian = [ensemble.TEMP.std_gaussian; prior_std];
                ensemble.TEMP.variable_type = [ensemble.TEMP.variable_type; 2]; %1 gaussian; 2: boxcar; 3: lognormal
                ensemble.TEMP.variable_name = [ensemble.TEMP.variable_name; ensemble.PARA.boxcar_variables_name{i,1}];
                ensemble.STATVAR.(ensemble.PARA.boxcar_variables_name{i,1}) = prior_transformed; %expit(ensemble, prior_untransformed, ensemble.PARA.boxcar_variables_lower_bound(i,1), ensemble.PARA.boxcar_variables_upper_bound(i,1));
%                 ensemble.STATVAR.(ensemble.PARA.boxcar_variables_name{i,1}) = ensemble.STATVAR.(ensemble.PARA.boxcar_variables_name{i,1})(1,tile.PARA.worker_number);
            end

            for i=1:size(ensemble.PARA.lognormal_variables_name,1)
                mean_gaussian = log(ensemble.PARA.lognormal_variables_mean(i,:).^2./sqrt(ensemble.PARA.lognormal_variables_mean(i,:).^2 + ensemble.PARA.lognormal_variables_std(i,:).^2));
                std_gaussian = sqrt(log(1 + ensemble.PARA.lognormal_variables_std(i,:).^2 ./ ensemble.PARA.lognormal_variables_mean(i,:).^2)); %tranform mean and std from real to Gaussian space
                value_gaussian = [];
                value_normal_space = [];
                for j=1:size(ensemble.PARA.lognormal_variables_mean,2)
                    random_gaussian = mean_gaussian(1,j) + std_gaussian(1,j) .* randn(1,tile.PARA.ensemble_size);
                    value_gaussian = [value_gaussian; random_gaussian];
                    value_normal_space = [value_normal_space; ensemble.PARA.lognormal_variables_lower_bound(i,j) + exp(random_gaussian)];
                end
                value_gaussian = value_gaussian(:)';
                value_normal_space = value_normal_space(:)';
                ensemble.TEMP.value_gaussian = [ensemble.TEMP.value_gaussian; value_gaussian];
                ensemble.TEMP.mean_gaussian = [ensemble.TEMP.mean_gaussian; mean_gaussian];
                ensemble.TEMP.std_gaussian = [ensemble.TEMP.std_gaussian; std_gaussian];
                ensemble.TEMP.variable_type = [ensemble.TEMP.variable_type; 3]; %1 gaussian; 2: boxcar; 3: lognormal
                ensemble.TEMP.variable_name = [ensemble.TEMP.variable_name; ensemble.PARA.lognormal_variables_name{i,1}];
                
                ensemble.STATVAR.(ensemble.PARA.lognormal_variables_name{i,1}) = value_normal_space; %ensemble.PARA.lognormal_variables_lower_bound(i,1) + exp(value_gaussian);
%                 ensemble.STATVAR.(ensemble.PARA.lognormal_variables_name{i,1}) = ensemble.STATVAR.(ensemble.PARA.lognormal_variables_name{i,1})(1,tile.PARA.worker_number);
            end
        end
        
        function ensemble = write2provider(ensemble, tile)
            %write variables in PROVIDER
            for i=1:size(ensemble.PARA.modify_class_name,1)
                tile.RUN_INFO.PPROVIDER.CLASSES.(ensemble.PARA.modify_class_name{i,1}){ensemble.PARA.modify_class_index(i,1),1}.PARA.(ensemble.PARA.variable_in_class{i,1}) = ...
                    ensemble.STATVAR.(ensemble.PARA.variable_in_ensemble{i,1});
                tile.RUN_INFO.PPROVIDER.CLASSES.(ensemble.PARA.modify_class_name{i,1}){ensemble.PARA.modify_class_index(i,1),1}.PARA.ensemble_size = tile.PARA.ensemble_size; 
            end 
        end
        
        %called by DA
        function result = get_best_fitting_parameters(ensemble, tile, variable_list, weights)
            clipping_threshold = 0.95;

            % [~, pos] = max(weights,[],2);
            not_valid = sum(double(isnan(weights)),2)>0;
            for i=1:size(variable_list,1)
                pos_in_var_list = find(strcmp(ensemble.TEMP.variable_name, variable_list{i,1}));
                if ensemble.TEMP.variable_type(pos_in_var_list,1)==1
                    parameter = reshape(ensemble.TEMP.value_gaussian(pos_in_var_list,:), tile.PARA.ensemble_size, size(ensemble.TEMP.value_gaussian,2) ./ tile.PARA.ensemble_size);
                    
                    result.(variable_list{i,1}) = get_weighted_average_of_parameter(ensemble, parameter, weights, clipping_threshold);
                    %diag(parameter(:,pos));
                elseif ensemble.TEMP.variable_type(pos_in_var_list,1)==2
                    parameter = reshape(ensemble.TEMP.value_gaussian(pos_in_var_list,:), tile.PARA.ensemble_size, size(ensemble.TEMP.value_gaussian,2) ./ tile.PARA.ensemble_size);
                    optim_param = get_weighted_average_of_parameter(ensemble, parameter, weights, clipping_threshold);

                    parameter_transformed=[];
                    for j=1:size(pos,1)
                        parameter_transformed=[parameter_transformed; expit(ensemble, optim_param(j,1), ensemble.PARA.boxcar_variables_lower_bound(pos_in_var_list,j), ensemble.PARA.boxcar_variables_upper_bound(pos_in_var_list,j))];
                    end
                    result.(variable_list{i,1}) = parameter_transformed;
                elseif ensemble.TEMP.variable_type(pos_in_var_list,1)==3
                    parameter = reshape(ensemble.TEMP.value_gaussian(pos_in_var_list,:), tile.PARA.ensemble_size, size(ensemble.TEMP.value_gaussian,2) ./ tile.PARA.ensemble_size);
                    result.(variable_list{i,1}) = exp(get_weighted_average_of_parameter(ensemble, parameter, weights, clipping_threshold));
                end
                result.(variable_list{i,1})(not_valid,1) = NaN;
            end
           
        end

        function optim_param = get_weighted_average_of_parameter(ensemble, parameter, weights, clipping_threshold)

            [weights_sorted,weights_sorted_pos] = sort(weights,2, 'descend');
            weights_sorted_accum = cumsum(weights_sorted,2);
            clip_pos = [];
            for i=1:size(weights_sorted_accum,1)
                pos=find(weights_sorted_accum(i,:)>clipping_threshold, 1);
                if ~isempty(pos)
                    clip_pos=[clip_pos; pos];
                else
                    clip_pos=[clip_pos; size(weights_sorted_accum,2)];
                end
            end

            for i=1:size(weights_sorted,1)
                weights_sorted(i,clip_pos(i,1)+1:end)=0;
            end
            weights_sorted = weights_sorted./repmat(sum(weights_sorted,2), 1, size(weights_sorted,2));
            optim_param = [];

            for i=1:size(weights_sorted,1)
                optim_param = [optim_param; sum(parameter(i, weights_sorted_pos(i,:)).* weights_sorted(i,:)) ];
            end

        end


    end
end

