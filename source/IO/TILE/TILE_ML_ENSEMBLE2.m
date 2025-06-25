%========================================================================
% CryoGrid TILE class TILE_ML_ENSEMBLE
% TILE class designed for machine learning 

% S. Westermann, April 2025
%========================================================================

classdef TILE_ML_ENSEMBLE2 < matlab.mixin.Copyable
    
    properties
        
        PARA
        RUN_INFO
        FORCING
        CONST
        OUT        
        STORE
        TEMP
        STATVAR
        ML
        ML_STORE
        t
        FEATURE_CLASS
        TRANSFORM_ML_out
        TRANSFORM_ML_in
    end
    
    
    methods
        
       

        function tile = provide_PARA(tile)

            tile.PARA.number_of_MLs = [];
            
            tile.PARA.ml_class = [];
            tile.PARA.ml_class_index = [];

            tile.PARA.training_class = [];
            tile.PARA.training_class_index = [];

            tile.PARA.target_data_class = [];
            tile.PARA.target_data_class_index = [];

            tile.PARA.feature_data_class = [];
            tile.PARA.feature_data_class_index = [];

            tile.PARA.in2out_ML_features_class = [];%must be stored
            tile.PARA.in2out_ML_features_class_index = [];

            tile.PARA.in2out_ML_target_class = [];
            tile.PARA.in2out_ML_target_class_index = [];

            tile.PARA.out_class = [];
            tile.PARA.out_class_index = [];

            tile.PARA.read_from_store = 0;
            tile.PARA.strip4store = 0;
            
            tile.PARA.worker_number = 1;
        end

        function tile = provide_CONST(tile)

        end
        
        function tile = provide_STATVAR(tile)

        end

        function tile = finalize_init(tile)

            tile.PARA.run_name =  tile.RUN_INFO.PPROVIDER.PARA.run_name;
            tile.PARA.result_path =  tile.RUN_INFO.PPROVIDER.PARA.result_path;
            
            disp('getting training data')

            rng(1234+tile.PARA.worker_number)
            target_data_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.target_data_class){tile.PARA.target_data_class_index,1}); %read out
            target_data_class = finalize_init(target_data_class, tile);

            feature_data_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.feature_data_class){tile.PARA.feature_data_class_index,1}); %read terrain and forcing
            feature_data_class = finalize_init(feature_data_class, tile);
            [feature_data_in_ML, feature_data_groups] = generate_feature_data_training(feature_data_class, tile);
            tile.FEATURE_CLASS = feature_data_class;
            [target_ML, target_data_groups] = generate_target_data(target_data_class, tile);

            if tile.PARA.read_from_store
                tile.ML_STORE = tile.RUN_INFO.TILE.ML_STORE; %this is the data in real space
                feature_data_in_ML = [tile.ML_STORE.in; feature_data_in_ML];
                tile.ML_STORE.in = feature_data_in_ML;
                target_ML = [tile.ML_STORE.out; target_ML];
                tile.ML_STORE.out = target_ML;
            else
                tile.ML_STORE.in = feature_data_in_ML;
                tile.ML_STORE.out = target_ML;
            end

            tile.TRANSFORM_ML_out = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.in2out_ML_target_class){tile.PARA.in2out_ML_target_class_index,1}); %real data to NN input and NN output to real data
            tile.TRANSFORM_ML_out.STATVAR.data_groups = target_data_groups;
            tile.TRANSFORM_ML_out.STATVAR.data = target_ML;
            tile.TRANSFORM_ML_out = finalize_init(tile.TRANSFORM_ML_out, tile); %calculate the STDs
            tile.TRANSFORM_ML_out = calculate_mean_std_derivatives(tile.TRANSFORM_ML_out, target_ML, tile);

            tile.TRANSFORM_ML_in = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.in2out_ML_features_class){tile.PARA.in2out_ML_features_class_index,1}); %real data to NN input and NN output to real data
            tile.TRANSFORM_ML_in.STATVAR.data_groups = feature_data_groups;
            tile.TRANSFORM_ML_in.STATVAR.data = feature_data_in_ML;
            tile.TRANSFORM_ML_in = finalize_init(tile.TRANSFORM_ML_in, tile);%calculate the STDs

            disp('training neural net')

            tile.STATVAR.in = feature_data_in_ML;
            tile.STATVAR.out = target_ML;
            
            for j=1:tile.PARA.number_of_MLs
                tile.ML = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.ml_class){tile.PARA.ml_class_index,1});
                tile.ML = finalize_init(tile.ML, tile);

                training_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.training_class){tile.PARA.training_class_index,1});
                training_class = finalize_init(training_class, tile);
                tile.ML = train_ML(training_class, tile);
                tile.ML_STORE.ML{1,j} = tile.ML;
            end
            tile.STATVAR.out = target_ML;

            tile = assign_timestamp_prediction(feature_data_class, tile);

            tile.OUT = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.out_class){tile.PARA.out_class_index,1});
            tile.OUT = finalize_init(tile.OUT, tile);

        end

        function tile = run_model(tile)
            disp('running neural net forward')
            for i = 1:size(tile.STATVAR.timestamp_prediction,1)
                tile.t = tile.STATVAR.timestamp_prediction(i,1);
                input_ML = generate_feature_data_prediction_single_timestamp(tile.FEATURE_CLASS, tile);
                input_ML = realWorld2NNinput(tile.TRANSFORM_ML_in, input_ML, tile);

                tile.STATVAR.ML_out = [];
                % tile.STATVAR.ML_out_std = []; store std of NN ensemble at the DA level? -> not needed now yet
                
                for j=1:tile.PARA.number_of_MLs
                    tile.ML = tile.ML_STORE.ML{1,j};
                    [out_j, ~] = progapagate_ML(tile.ML, input_ML);
                    ML_ensemble_realWorld = [];
                    for ii=1:size(out_j,2)
                        pe_ML =  NNoutput2realWorld(tile.TRANSFORM_ML_out, out_j(:,ii), tile);  %transform to real space, then to  target space, like calculating derivatives
                        ML_ensemble_realWorld = cat(3, ML_ensemble_realWorld, pe_ML); %concatenate members of ensemble for individual NN outoputs in 3rd direction 
                    end
                    tile.STATVAR.ML_out = cat(3, tile.STATVAR.ML_out,  mean(ML_ensemble_realWorld,3));
                 %   tile.STATVAR.ML_out_std = cat(3, tile.STATVAR.ML_out_std, std(ML_ensemble_realWorld, [], 3));
                end
                tile = store_OUT_tile(tile);
            end

            tile = write_OUT_tile(tile); %write to SPATIAL and excahnge info between tiles if parallel
            
            if tile.PARA.strip4store
                tile.OUT = [];
                tile.ML = [];
                tile.FEATURE_CLASS = [];
                tile.TRANSFORM_ML_out = [];
                tile.TRANSFORM_ML_in = [];
                tile.PARA = [];
                tile.RUN_INFO = [];
                tile.FORCING = [];
                tile.CONST = [];
                tile.STORE = [];
                tile.TEMP = [];
                tile.STATVAR = [];
                tile.t = [];
            end
        end
        

        function tile = store_OUT_tile(tile)
            tile.OUT = store_OUT(tile.OUT, tile);
        end        

        function tile = write_OUT_tile(tile)
            tile.OUT = write_OUT(tile.OUT, tile);
        end  
        
  

    end
end



