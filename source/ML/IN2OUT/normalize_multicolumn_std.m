classdef normalize_multicolumn_std < matlab.mixin.Copyable
    
    properties
        PARA
        STATVAR
        CONST
    end

    methods
        function in_out = provide_PARA(in_out)

        end

        function in_out = provide_CONST(in_out)

        end
        function in_out = provide_STATVAR(in_out)

        end

        function in_out = finalize_init(in_out, tile)

            in_out.STATVAR.data_size = size(in_out.STATVAR.data);       

            in_out.STATVAR.data_mean{1,1} = in_out.PARA.data_groups.*0;
            in_out.STATVAR.data_std{1,1} = in_out.PARA.data_groups.*0;
            for i=1:max(in_out.PARA.data_groups(:))
                range = in_out.PARA.data_groups==i;
                data = in_out.STATVAR.data(:,range);
                in_out.STATVAR.data_mean{1,1}(range) = mean(data(:)); %mean and std can be used as normal in the equation, no need to reshape altough it is 3D
                in_out.STATVAR.data_std{1,1}(range) = std(data(:));
            end
        end


        function out = realWorld2NN(in_out, in, tile)
            out = (in - in_out.STATVAR.data_mean{1,1}) ./ in_out.STATVAR.data_std{1,1};
        end


        function [in_out, target] = realWorld2training(in_out, in, tile) %This must calculate the derivatives, probably good in several stages
            target.values{1,1} = in;
            % target.values{2,1} = in(:,1:end-1,:)-in(:,2:end,:);
            % target.values{3,1} = in(:,:,1:end-1)-in(:,:,2:end);
            in_out.STATVAR.sizes{1,1} = size(target.values{1,1},1,2,3);
            % in_out.STATVAR.sizes{2,1} = size(target.values{2,1},1,2,3);
            % in_out.STATVAR.sizes{3,1} = size(target.values{3,1},1,2,3);
            target.values{1,1} = (target.values{1,1} - in_out.STATVAR.data_mean{1,1}) ./ in_out.STATVAR.data_std{1,1};
            % in_out.STATVAR.data_mean{2,1} = mean(target.values{2,1}(:));
            % in_out.STATVAR.data_std{2,1} = std(target.values{2,1}(:));
            % target.values{2,1} = (target.values{2,1} - in_out.STATVAR.data_mean{2,1}) ./ in_out.STATVAR.data_std{2,1};
            % in_out.STATVAR.data_mean{3,1} = mean(target.values{3,1}(:));
            % in_out.STATVAR.data_std{3,1} = std(target.values{3,1}(:));
            % target.values{3,1} = (target.values{3,1} - in_out.STATVAR.data_mean{3,1}) ./ in_out.STATVAR.data_std{3,1};            

            %target = [target.values{1,1}(:); target.values{2,1}(:); target.values{3,1}(:)];
            target = target.values{1,1}(:);
        end

        function uncertainty = realWorld2training_uncertainty(in_out, relative_error_term, tile)
            uncertainty.values{1,1} = relative_error_term ./ repmat(in_out.STATVAR.data_std{1,1}, in_out.STATVAR.sizes{1,1}(1)./size(in_out.STATVAR.data_std{1,1},1), ...
                in_out.STATVAR.sizes{1,1}(2)./size(in_out.STATVAR.data_std{1,1},2), in_out.STATVAR.sizes{1,1}(3)./size(in_out.STATVAR.data_std{1,1},3));
            % uncertainty.values{2,1} = relative_error_term .*100./ repmat(in_out.STATVAR.data_std{2,1}, in_out.STATVAR.sizes{2,1}(1)./size(in_out.STATVAR.data_std{2,1},1), ...
            %     in_out.STATVAR.sizes{2,1}(2)./size(in_out.STATVAR.data_std{2,1},2), in_out.STATVAR.sizes{2,1}(3)./size(in_out.STATVAR.data_std{2,1},3));
            % uncertainty.values{3,1} = relative_error_term .*100./ repmat(in_out.STATVAR.data_std{3,1}, in_out.STATVAR.sizes{3,1}(1)./size(in_out.STATVAR.data_std{3,1},1), ...
            %     in_out.STATVAR.sizes{3,1}(2)./size(in_out.STATVAR.data_std{3,1},2), in_out.STATVAR.sizes{3,1}(3)./size(in_out.STATVAR.data_std{3,1},3));

            %uncertainty = [uncertainty.values{1,1}(:); uncertainty.values{2,1}(:); uncertainty.values{3,1}(:)];
            uncertainty = uncertainty.values{1,1}(:); 
            uncertainty = uncertainty .^2;
        end

        function target = NNout2training(in_out, in, tile) 
            target.values{1,1} = reshape(in, [in_out.STATVAR.sizes{1,1} size(in,2)]);
            target.values{1,1} = target.values{1,1} .*  in_out.STATVAR.data_std{1,1} + in_out.STATVAR.data_mean{1,1};
            % target.values{2,1} = target.values{1,1}(:,1:end-1,:,:) - target.values{1,1}(:,2:end,:,:);
            % target.values{3,1} = target.values{1,1}(:,:,1:end-1,:) -  target.values{1,1}(:,:,2:end,:);

            target.values{1,1} = (target.values{1,1} - in_out.STATVAR.data_mean{1,1}) ./ in_out.STATVAR.data_std{1,1};
            % target.values{2,1} = (target.values{2,1} - in_out.STATVAR.data_mean{2,1}) ./ in_out.STATVAR.data_std{2,1};
            % target.values{3,1} = (target.values{3,1} - in_out.STATVAR.data_mean{3,1}) ./ in_out.STATVAR.data_std{3,1};            

            % target = [reshape(target.values{1,1}, size(target.values{1,1}(:),1)./size(in,2), size(in,2)); ...
            %     reshape(target.values{2,1}, size(target.values{2,1}(:),1)./size(in,2), size(in,2)); ...
            %     reshape(target.values{3,1}, size(target.values{3,1}(:),1)./size(in,2), size(in,2))];
            target = reshape(target.values{1,1}, size(target.values{1,1}(:),1)./size(in,2), size(in,2));
        end
        
        
        % function out = real_world2NN_normalize_reshape(in_out, in, tile)
        %     out = (in - in_out.STATVAR.data_mean) ./ in_out.STATVAR.data_std;
        %     out = out(:);
        % end        
        % 
        % function out = NN2real_world_normalize(in_out, in, tile)
        %     in = reshape(in, in_out.STATVAR.var_size);
        %     out = in .* in_out.STATVAR.data_std + in_out.STATVAR.data_mean;
        % end
        % 
        % function out = real_world2NN_no_normalize(in_out, in, tile)
        %     out = in(:);
        % end

        function out = NN2real_world_no_normalize(in_out, in, tile)
            in = reshape(in, in_out.STATVAR.var_size);
        end
    end
end

