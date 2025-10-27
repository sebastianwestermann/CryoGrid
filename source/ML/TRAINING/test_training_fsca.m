da = da_store;

target_ML = reshape(da.ENSEMBLE.modeled_obs, size(da.ENSEMBLE.modeled_obs,1), size(da.ENSEMBLE.modeled_obs,2).*size(da.ENSEMBLE.modeled_obs,3))';
feature_data_in_ML = reshape(da.ENSEMBLE.value_gaussian, size(da.ENSEMBLE.value_gaussian,1), size(da.ENSEMBLE.value_gaussian,2).*size(da.ENSEMBLE.value_gaussian,3))';
%%
da.STATVAR.out_mean = mean(target_ML(:),1);
da.STATVAR.out_std = std(target_ML(:),[],1);
da.STATVAR.in_mean = mean(feature_data_in_ML,1);
da.STATVAR.in_std = std(feature_data_in_ML,[],1);
%da.STATVAR.out = (target_ML(:) - da.STATVAR.out_mean) ./ da.STATVAR.out_std;
da.STATVAR.out = target_ML(:);

% test=[ones(size(target_ML,1),5) target_ML zeros(size(target_ML,1),5)];
% score = [];
% for i=1:size(target_ML,2)
% j=i+5;
% score = [score mean(test(:,j-5:j),2) - mean(test(:,j:j+5),2)];
% end
% [~, max_id] = max(score, [], 2);
% da.STATVAR.out_mean = mean(max_id,1);
% da.STATVAR.out_std = std(max_id,[],1);
% 
% da.STATVAR.out = (max_id - da.STATVAR.out_mean) ./ da.STATVAR.out_std;

%da.STATVAR.out = max_id;  %do not normalize so that standard exp function can become active 
da.STATVAR.in = (feature_data_in_ML - da.STATVAR.in_mean) ./ da.STATVAR.in_std;
ml =  neural_net_ensemble();
ml.PARA.activation_functions = {'ReLU'; 'ReLU'; 'ReLU'; 'standard_logistic'};
ml.PARA.number_of_neurons = [15; 15; 15; 1];
ml.PARA.ensemble_size = 100;
ml = finalize_init(ml, da);
training_class = train_EnKA();
training_class.PARA.training_fraction = 1;
training_class.PARA.number_of_iterations = 100;
training_class.PARA.relative_error_term = 1e-3;
da.TEMP.var_ID = 1;
training_class = finalize_init(training_class, da);

ml = train_ML2(training_class, ml, da);

[predicted_ensemble_ML, ml_parameters] =  progapagate_ML(ml, (feature_data_in_ML2 - da.STATVAR.in_mean) ./ da.STATVAR.in_std);

%predicted_ensemble_ML= predicted_ensemble_ML.*da.STATVAR.out_std + da.STATVAR.out_mean;
predicted_ensemble = mean(predicted_ensemble_ML,2);
% plot(predicted_ensemble, 'red')
% hold on
