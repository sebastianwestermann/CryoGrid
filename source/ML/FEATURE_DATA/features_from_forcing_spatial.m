classdef features_from_forcing_spatial < matlab.mixin.Copyable
    
    properties
        PARA
        STATVAR
        CONST
        TEMP
        FORCING
        DATA
    end

    methods
        function features = provide_PARA(features)
            features.PARA.forcing_class = []; 
            features.PARA.forcing_class_index = [];

            features.PARA.averaging_interval = [];

            features.PARA.start_date_training=[];
            features.PARA.end_date_training=[];

            features.PARA.start_date_prediction=[];
            features.PARA.end_date_prediction=[];
            
            features.PARA.forcing2features_class = [];
            features.PARA.forcing2features_class_index = [];

            features.PARA.spatial2features_class = [];
            features.PARA.spatial2features_class_index = [];

            features.PARA.mask_class_training = [];
            features.PARA.mask_class_training_index = [];

            features.PARA.mask_class_prediction = [];
            features.PARA.mask_class_prediction_index = [];
            %features.PARA.range = [];
        end


        function features = provide_CONST(features)

        end

        function features = provide_STATVAR(features)

        end

        function features = finalize_init(features, tile)
            %establish time raster
            if strcmp(features.PARA.averaging_interval, 'annual')
                features.DATA.timestamp = datenum(features.PARA.start_date_training(1,1), features.PARA.start_date_training(2,1), features.PARA.start_date_training(3,1));
                current_date = datenum(year(features.DATA.timestamp)+1, features.PARA.start_date_training(2,1), features.PARA.start_date_training(3,1));
                while current_date <= datenum(features.PARA.end_date_training(1,1), features.PARA.start_date_training(2,1), features.PARA.start_date_training(3,1))
                    features.DATA.timestamp = [features.DATA.timestamp; current_date];
                    current_date = datenum(year(current_date)+1, features.PARA.start_date_training(2,1), features.PARA.start_date_training(3,1));
                end
            elseif strcmp(features.PARA.averaging_interval, 'monthly')
                %to DO
            elseif strcmp(features.PARA.averaging_interval, 'daily')
                features.DATA.timestamp = [datenum(features.PARA.start_date_training(1,1), features.PARA.start_date_training(2,1), features.PARA.start_date_training(3,1)): ...
                    datenum(features.PARA.end_date_training(1,1), features.PARA.end_date_training(2,1), features.PARA.end_date_training(3,1))]';
            end

            features.DATA.features_forcing = [];
            %average
            for i=1:size(features.PARA.forcing_class,1)
                forcing = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(features.PARA.forcing_class{i,1}){features.PARA.forcing_class_index(i,1),1}); %read out
                forcing = finalize_init(forcing, tile);
                for j=1:size(features.PARA.forcing2features_class, 1)
                    forcing2features = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(features.PARA.forcing2features_class{j,1}){features.PARA.forcing2features_class_index(j,1),1});
                    forcing2features = finalize_init(forcing2features, tile);
                    features.DATA.features_forcing = [features.DATA.features_forcing features_from_forcing(forcing2features, features, forcing)];
                end
            end
        end

        function out = generate_feature_data_training(features, tile)
            %spatial features 
            valid = ones(size(tile.RUN_INFO.SPATIAL.STATVAR.latitude,1),1); 
            for i=1:size(features.PARA.mask_class_training,1)
                mask = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(features.PARA.mask_class_training{i,1}){features.PARA.mask_class_training_index(i,1),1});
                mask = finalize_init(mask, tile);
                valid = apply_mask2(mask, tile, valid); %both additive and exclusive masks are possible
            end
            features_spatial = [];
            for i=1:size(features.PARA.spatial2features_class, 1)
                spatial2features = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(features.PARA.spatial2features_class{i,1}){features.PARA.spatial2features_class_index(i,1),1});
                spatial2features = finalize_init(spatial2features, tile);
                fs = features_from_spatial(spatial2features, tile, valid);
                features_spatial = [features_spatial fs];
            end
            features.DATA.features_spatial = features_spatial;

            %combine with forcing
            out = repmat(features.DATA.features_forcing, size(features_spatial,1),1);
            out_spatial = [];
            for j=1:size(features_spatial,1)
                out_spatial = [out_spatial; repmat(features_spatial(j,:), size(features.DATA.features_forcing,1),1)];
            end
            out = [out out_spatial];
        end

        function tile = assign_timestamp_prediction(features, tile)
            valid_forcing = features.DATA.timestamp >= datenum(features.PARA.start_date_prediction(1,1), features.PARA.start_date_prediction(2,1), features.PARA.start_date_prediction(3,1)) ...
                & features.DATA.timestamp <= datenum(features.PARA.end_date_prediction(1,1), features.PARA.end_date_prediction(2,1), features.PARA.end_date_prediction(3,1));
            tile.STATVAR.timestamp_prediction = features.DATA.timestamp(valid_forcing,1);
        end

        function out = generate_feature_data_prediction_single_timestamp(features, tile)

            %spatial features
            valid = zeros(size(tile.RUN_INFO.SPATIAL.STATVAR.latitude,1),1); 
            valid(tile.PARA.range,1) = 1;
            for i=1:size(features.PARA.mask_class_prediction,1)
                mask = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(features.PARA.mask_class_prediction{i,1}){features.PARA.mask_class_prediction_index(i,1),1});
                mask = finalize_init(mask, tile);
                valid = apply_mask2(mask, tile, valid); %both additive and exclusive masks are possible
            end
            valid = logical(valid);
            features_spatial = [];
            for i=1:size(features.PARA.spatial2features_class, 1)
                spatial2features = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(features.PARA.spatial2features_class{i,1}){features.PARA.spatial2features_class_index(i,1),1});
                spatial2features = finalize_init(spatial2features, tile);
                fs = features_from_spatial(spatial2features, tile, valid);
                features_spatial = [features_spatial fs];
            end

            %combine with forcing
            % valid_forcing = features.DATA.timestamp >= datenum(features.PARA.start_date_prediction(1,1), features.PARA.start_date_prediction(2,1), features.PARA.start_date_prediction(3,1)) ...
            %     & features.DATA.timestamp <= datenum(features.PARA.end_date_prediction(1,1), features.PARA.end_date_prediction(2,1), features.PARA.end_date_prediction(3,1));
            valid_forcing = find(abs(features.DATA.timestamp-tile.t) <1e-5); 
            out = repmat(features.DATA.features_forcing(valid_forcing,:), size(features_spatial,1),1);
            out = [out features_spatial];
          
        end

    end
end

