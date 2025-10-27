classdef normalize_multicolumn_derivatives_time_depth < matlab.mixin.Copyable
    
    properties
        PARA
        STATVAR
        CONST
    end

    methods
        function in_out = provide_PARA(in_out)
            in_out.PARA.number_of_time_slices = [];
            in_out.PARA.variance_inflation_time_derivatives = []; %variance (not std!) multiplied with this factor
            in_out.PARA.variance_inflation_depth_derivatives = [];
        end

        function in_out = provide_CONST(in_out)

        end
        function in_out = provide_STATVAR(in_out)

        end

        function in_out = finalize_init(in_out, tile) %->CORRECTED

            in_out.STATVAR.data_size = size(in_out.STATVAR.data);       

            in_out.STATVAR.data_mean = [1:max(in_out.STATVAR.data_groups)].*NaN;
            in_out.STATVAR.data_std = [1:max(in_out.STATVAR.data_groups)].*NaN;
            for i=1:max(in_out.STATVAR.data_groups)
                range = in_out.STATVAR.data_groups==i;
                data = in_out.STATVAR.data(:,range);
                in_out.STATVAR.data_mean(1,i) = mean(data(:)); %mean and std can be used as normal in the equation, no need to reshape altough it is 3D
                in_out.STATVAR.data_std(1,i) = std(data(:));
            end
        end


        function out = realWorld2NNinput(in_out, in, tile) %->CORRECTED
            out = in.*NaN;
            for i=1:max(in_out.STATVAR.data_groups)
                range = find(in_out.STATVAR.data_groups==i);
                out(:,range) = (in(:,range) - in_out.STATVAR.data_mean(1, i)) ./ in_out.STATVAR.data_std(1,i);
            end
        end

        function out = NNoutput2realWorld(in_out, in, tile) %->CORRECTED
            %reshape
            in = reshape(in, size(in,1)./in_out.STATVAR.data_size(2), in_out.STATVAR.data_size(2));
            for i=1:max(in_out.STATVAR.data_groups)
                range = find(in_out.STATVAR.data_groups==i);
                out(:, range) = in(:, range) .* in_out.STATVAR.data_std(1,i) + in_out.STATVAR.data_mean(1, i);
            end
        end

        function in_out = calculate_mean_std_derivatives(in_out, in, tile)
            derivatives_time = [];
            for i=1:in_out.PARA.number_of_time_slices:size(in,1)
                in_ts = in(i:i+in_out.PARA.number_of_time_slices-1,:);
                derivatives_time = [derivatives_time; in_ts(2:end,:) - in_ts(1:end-1,:)];
            end
            in_out.STATVAR.derivatives_time_mean = mean(derivatives_time,1);
            in_out.STATVAR.derivatives_time_std = std(derivatives_time,1);
            in_out.STATVAR.derivatives_time_size = size(derivatives_time);

            derivatives_depth = in(:,2:end) - in(:,1:end-1);
            in_out.STATVAR.derivatives_depth_mean = mean(derivatives_depth,1);
            in_out.STATVAR.derivatives_depth_std = std(derivatives_depth,1);
            in_out.STATVAR.derivatives_depth_size = size(derivatives_depth);
        end

        function [out, in_out] = realWorld2training(in_out, in, tile) %This must calculate the derivatives, probably good in several stages
            derivatives_time = [];
            for i=1:in_out.PARA.number_of_time_slices:size(in,1)
                in_ts = in(i:i+in_out.PARA.number_of_time_slices-1,:);
                derivatives_time = [derivatives_time; in_ts(2:end,:) - in_ts(1:end-1,:)];
            end
            derivatives_time = (derivatives_time - in_out.STATVAR.derivatives_time_mean) ./ in_out.STATVAR.derivatives_time_std;
            derivatives_depth = in(:,2:end) - in(:,1:end-1);
            derivatives_depth = (derivatives_depth - in_out.STATVAR.derivatives_depth_mean) ./ in_out.STATVAR.derivatives_depth_std;

            out = in.*NaN;
            for i=1:max(in_out.STATVAR.data_groups)
                range = find(in_out.STATVAR.data_groups==i);
                out(:,range) = (in(:,range) - in_out.STATVAR.data_mean(1, i)) ./ in_out.STATVAR.data_std(1,i);
            end
            out = [out(:); derivatives_time(:); derivatives_depth(:)];
        end

        function variance = realWorld2training_variance(in_out, error_term, tile)

            data_uncertainty = in_out.STATVAR.data_groups .*0;
            for i=1:max(in_out.STATVAR.data_groups)
                range = find(in_out.STATVAR.data_groups==i);
                data_uncertainty(1,range) = error_term ./ in_out.STATVAR.data_std(1,i);
            end
            data_uncertainty = repmat(data_uncertainty, in_out.STATVAR.data_size(1), 1);
             
            derivatives_time_uncertainty = error_term ./ in_out.STATVAR.derivatives_time_std;
            derivatives_time_uncertainty = repmat(derivatives_time_uncertainty, in_out.STATVAR.derivatives_time_size(1),1);

            derivatives_depth_uncertainty = error_term ./ in_out.STATVAR.derivatives_depth_std;
            derivatives_depth_uncertainty = repmat(derivatives_depth_uncertainty, in_out.STATVAR.derivatives_depth_size(1),1);

            variance = [data_uncertainty(:).^2; in_out.PARA.variance_inflation_time_derivatives .* derivatives_time_uncertainty(:).^2; ...
                in_out.PARA.variance_inflation_depth_derivatives.*derivatives_depth_uncertainty(:).^2];

        end


    end
end

