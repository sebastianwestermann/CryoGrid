%========================================================================
% CryoGrid CLUSTERING class K_MEANS_add_clusters
% CLUSTERING class applying K-means clustering to spatial data
%
% S. Westermann, Aug 2025
%========================================================================

classdef K_MEANS_add_clusters < matlab.mixin.Copyable

    properties
        RUN_INFO
        PARA
        CONST
        STATVAR
    end
    
    methods
        function cluster = provide_PARA(cluster)
            cluster.PARA.number_of_clusters = [];
            cluster.PARA.cluster_variable_class = [];
            cluster.PARA.cluster_variable_class_index = [];
            cluster.PARA.weighting_variables = [];
            cluster.PARA.inflation_factor_weighting = [];
            cluster.PARA.add2SPATIAL = 0;
        end
        
        function cluster = provide_STATVAR(cluster)

        end
        
        function cluster = provide_CONST(cluster)
            
        end
        
        function cluster = finalize_init(cluster)
            
        end
        
        function cluster = compute_clusters(cluster)
            disp('computing clusters')
            data_cube = [];
            for i=1:size(cluster.PARA.cluster_variable_class,1)
                cluster_variable_class = copy(cluster.RUN_INFO.PPROVIDER.CLASSES.(cluster.PARA.cluster_variable_class{i,1}){cluster.PARA.cluster_variable_class_index(i,1),1});
                cluster_variable_class.SPATIAL = cluster.RUN_INFO.SPATIAL;
                data_cube = [data_cube provide_data(cluster_variable_class)];
            end

            Zs = zscore(data_cube,1,1); % Does this make sense for bounded variables, or should you transform first?
            
            sample_centroid_index = find(cluster.RUN_INFO.SPATIAL.STATVAR.is_cluster_centroid==1);
            cluster.STATVAR.sample_centroid_index = [];
            test=[];

            wn = spmdIndex;
            save(['test_' num2str(wn) '_' datestr(now, 'mm_dd_HH_MM') '.mat'])
            
            for it=1:cluster.PARA.number_of_clusters
                average_distance_to_cluster_center = zeros(size(data_cube,1),1)+Inf;
                max_distance_to_cluster_center = zeros(size(data_cube,1),1)+Inf;
                for i=1:size(sample_centroid_index,1)
                    average_distance_to_cluster_center = min(average_distance_to_cluster_center, sqrt(sum( (Zs-repmat(Zs(sample_centroid_index(i,1),:),size(data_cube,1),1)).^2,2)));
                    max_distance_to_cluster_center = min(max_distance_to_cluster_center, sqrt( max( (Zs - repmat(Zs(sample_centroid_index(i,1),:),size(data_cube,1),1)).^2, [], 2)));
                end
                score = (average_distance_to_cluster_center ./mean(average_distance_to_cluster_center,1)).^2 + (max_distance_to_cluster_center./mean(max_distance_to_cluster_center,1)).^2;
                score = sqrt(score);
                score = score ./ mean(score,1);
                score_save = score;
                for j=1:size(cluster.PARA.weighting_variables,1)
                    weighting_variable = cluster.RUN_INFO.SPATIAL.STATVAR.(cluster.PARA.weighting_variables{j,1})./mean(cluster.RUN_INFO.SPATIAL.STATVAR.(cluster.PARA.weighting_variables{j,1}),1);
                    score = score + cluster.PARA.inflation_factor_weighting .* sqrt(sum(weighting_variable.^2,2)) .*  double(score_save>1e-3); %last term eliminates the problem of getting existing cluster centers again
                end
                [~,id_new_centroid] = max(score);
                %test=[test;[score_save(id_new_centroid) weighting_variable(id_new_centroid)]];
                cluster.STATVAR.sample_centroid_index = [cluster.STATVAR.sample_centroid_index; id_new_centroid];
                sample_centroid_index = [sample_centroid_index; id_new_centroid];

            end
            
            cluster.STATVAR.key_centroid_index = cluster.RUN_INFO.SPATIAL.STATVAR.key(cluster.STATVAR.sample_centroid_index ,1);
            cluster.RUN_INFO.SPATIAL.STATVAR.is_cluster_centroid = cluster.RUN_INFO.SPATIAL.STATVAR.key .*0;
            cluster.RUN_INFO.SPATIAL.STATVAR.is_cluster_centroid(sample_centroid_index,1) = 1; %overwrite and add to existing

         end
        
        
         
         %-------------param file generation-----
         function cluster = param_file_info(cluster)
             cluster = provide_PARA(cluster);
             
             cluster.PARA.STATVAR = [];
             cluster.PARA.class_category = 'CLUSTERING';
             
             cluster.PARA.comment.number_of_clusters = {'number of clusters'};
             
             cluster.PARA.comment.max_iterations = {'maximum number of interations, interrupts k-means algorithm if no convergence is reached'};
             cluster.PARA.default_value.max_iterations = {1000};
             
             cluster.PARA.comment.cluster_variable_class = {'list of classes providing the data to which the clustering is applied'};
             cluster.PARA.options.cluster_variable_class.name = 'H_LIST';
             cluster.PARA.options.cluster_variable_class.entries_x = {'CLUSTER_RAW_VARIABLES' 'CLUSTER_SLOPE_ASPECT'};
             
             cluster.PARA.options.cluster_variable_class_index.name = 'H_LIST';
             cluster.PARA.options.cluster_variable_class_index.entries_x = {'1' '1'};
             
         end
         
    end
end

