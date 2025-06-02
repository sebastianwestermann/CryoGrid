classdef features_from_forcing_spatial2 < matlab.mixin.Copyable
    
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

            features.DATA.features_forcing = {}; %cell array, so that the std and mean can be independently computed for each variable separately
            features.DATA.features_forcing_mean = [];
            features.DATA.features_forcing_std = [];            
            %average
            for i=1:size(features.PARA.forcing_class,1)
                forcing = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(features.PARA.forcing_class{i,1}){features.PARA.forcing_class_index(i,1),1}); %read out
                forcing = finalize_init(forcing, tile);
                for j=1:size(features.PARA.forcing2features_class, 1)
                    forcing2features = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(features.PARA.forcing2features_class{j,1}){features.PARA.forcing2features_class_index(j,1),1});
                    forcing2features = finalize_init(forcing2features, tile);
                    feature_data = features_from_forcing(forcing2features, features, forcing);
                    features.DATA.features_forcing = [features.DATA.features_forcing; feature_data];
                end
            end
            % for i=1:size(features.DATA.features_forcing,1)
            %     features.DATA.features_forcing_mean = [features.DATA.features_forcing_mean; mean(features.DATA.features_forcing{i,1}(:))];
            %     features.DATA.features_forcing_std = [features.DATA.features_forcing_std; std(features.DATA.features_forcing{i,1}(:))];
            % end
        end

        function [out, data_groups] = generate_feature_data_training(features, tile)

            %rearrange forcing
            data_groups = [];
            features_forcing = [];
            for i=1:size(features.DATA.features_forcing,1)
                features_forcing = [features_forcing features.DATA.features_forcing{i,1}];
                data_groups = [data_groups ones(1,size(features.DATA.features_forcing{i,1},2)).*i];
            end
            maxID_data_groups = i;

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
            data_groups = [data_groups maxID_data_groups+[1:size(features_spatial,2)]];
            % features.DATA.features_spatial_mean = mean(features_spatial, 1);
            % features.DATA.features_spatial_std = std(features_spatial, [], 1);
            % features_spatial = (features_spatial - features.DATA.features_spatial_mean) ./ features.DATA.features_spatial_std;


            %combine with forcing
            out = repmat(features_forcing, size(features_spatial,1),1);
            out_spatial = [];
            for j=1:size(features_spatial,1)
                out_spatial = [out_spatial; repmat(features_spatial(j,:), size(features_forcing,1),1)];
            end
            out = [out out_spatial];
        end

        function out = generate_feature_data_prediction(features, tile)

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
            features_spatial = (features_spatial - features.DATA.features_spatial_mean) ./ features.DATA.features_spatial_std;

            %combine with forcing
            features_forcing = [];
            for i=1:size(features.DATA.features_forcing,1)
                features_forcing = [features_forcing (features.DATA.features_forcing{i,1} - features.DATA.features_forcing_mean(i,1)) ./ features.DATA.features_forcing_std(i,1)];
            end
            out = repmat(features_forcing, size(features_spatial,1),1);
            out = [out features_spatial];
          
        end

    end
end

