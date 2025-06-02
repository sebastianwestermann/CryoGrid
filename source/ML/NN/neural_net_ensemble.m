classdef neural_net_ensemble < matlab.mixin.Copyable
    
    properties
        PARA
        STATVAR
        CONST
        TEMP
        ACTIVATION
    end

    methods
        function nn = provide_PARA(nn)
            nn.PARA.activation_functions = []; 
            nn.PARA.number_of_neurons = [];
            nn.PARA.ensemble_size = [];
        end

        function nn = provide_CONST(nn)

        end
        function nn = provide_STATVAR(nn)

        end

        function nn = finalize_init(nn, tile)
           
            for i=1:size(nn.PARA.activation_functions, 1)
                af = str2func(nn.PARA.activation_functions{i,1});
                nn.ACTIVATION{i,1} = af();
            end
            nn.PARA.number_of_neurons(end, 1) = size(tile.TEMP.in2out_NN_target.STATVAR.data(:), 1) ./ size(tile.TEMP.in2out_NN_target.STATVAR.data, 1);
            neurons_per_layer = [size(tile.TEMP.in2out_NN_features.STATVAR.data(:), 1) ./ size(tile.TEMP.in2out_NN_features.STATVAR.data, 1); nn.PARA.number_of_neurons];
            size_weights = 0;
            for i=2:size(neurons_per_layer,1)
                size_weights = size_weights + neurons_per_layer(i-1).*neurons_per_layer(i);
            end
            nn.TEMP.neurons_per_layer = neurons_per_layer;
            nn.STATVAR.weights = randn(size_weights, nn.PARA.ensemble_size);
            nn.STATVAR.bias_terms = randn(sum(nn.PARA.number_of_neurons), nn.PARA.ensemble_size);
        end


        function [out_all, nn_parameters] = progapagate_ML(nn, in)

            out_all=zeros(size(in,1).*nn.PARA.number_of_neurons(end,1), nn.PARA.ensemble_size);
            for count = 1:nn.PARA.ensemble_size
                out = in;
                start_id_w=1;
                start_id_b=1;
                for i=1:size(nn.PARA.activation_functions,1)
                    end_id_w = start_id_w - 1 + nn.TEMP.neurons_per_layer(i,1) .* nn.TEMP.neurons_per_layer(i+1,1);
                    weights = reshape(nn.STATVAR.weights(start_id_w:end_id_w, count), nn.TEMP.neurons_per_layer(i,1), nn.TEMP.neurons_per_layer(i+1,1));
                    start_id_w = end_id_w + 1;
                    end_id_b = start_id_b - 1 + nn.TEMP.neurons_per_layer(i+1,1); %first
                    bias_term = nn.STATVAR.bias_terms(start_id_b:end_id_b,count);
                    start_id_b = end_id_b + 1;
                    out = out*weights+bias_term';
                    out = propagate(nn.ACTIVATION{i,1}, out);
                end
                %out = out'; %NEW: columns different values, like time dependence; rows: different samples
                out_all(:,count) = out(:);
            end
            nn_parameters = [nn.STATVAR.weights; nn.STATVAR.bias_terms];
        end


        function nn = reset_parameters(nn, new_parameters)
            nn.STATVAR.weights = new_parameters(1:size(nn.STATVAR.weights,1),:);
            nn.STATVAR.bias_terms = new_parameters(size(nn.STATVAR.weights,1)+1:end,:);
        end
    end
end


% ypred=forward_MLP(X,npl,afuncs,w,b)
%% Forward MLP 
% Inputs:
% X = Feature matrix (Ns x Nf) 
% npl = Vector containing number of neurons per layer ((Nl+1) x 1) 
% afuncs = Cell with activation functions for each layer ((Nl +1) x 1)
% w = Weight vector (concatenation of weight matrices for each layer) 
% b = Bias vector (concatenation of bias vector for each layer)
% Short hands: Ns = Number of sampled feature vectors, Nf = Number of
% features (inputs), Nl = Number of hidden layers.

% tanh seems promising: https://doi.org/10.1016/j.neunet.2021.08.015


%Ns=size(X,1);
% Nf=size(X,2);
% X=X'; % Transpose, so features as rows.
% 
% % Forward pass through layers.
% % Coded as explicitly as possible so it's possible to follow what's
% % happening.
% for j=1:numel(npl)
%     Nj=npl(j);
%     afj=afuncs{j};
%     if j==1
%         Nprev=Nf;
%         Zold=X;
%         wind=0;
%         bind=0;
%     else
%         Nprev=npl(j-1);
%         Zold=Zj;
%         wind=max(thesew);
%         bind=max(theseb); %bind
% 
%     end
%     thesew=1:(Nprev*Nj);
%     thesew=thesew+wind;
%     Wj=w(thesew);
%     if Nprev>1
%         Wj=reshape(Wj,Nj,Nprev);
%     end
%     theseb=1:Nj;
%     theseb=theseb+bind;
%     bj=b(theseb);
% 
% 
%     Aj=Wj*Zold+bj; % Activity of this layer
%     % Note, bj addition still works for each sample through broadcasting.
% 
%     % Activation 
%     % There's a whole zoo of these, check: https://en.wikipedia.org/wiki/Activation_function
%     if strcmp('ReLU',afj)
%         Zj=max(Aj,0);
%     elseif strcmp('LeakyReLU',afj)
%         here=Aj>0;
%         Zj=zeros(size(Aj));
%         Zj=Aj;
%         neg=Aj<0;
%         Zj(neg)=0.01.*Aj(neg);
%         %Zj(here)=Aj(here);
%         %Zj(~here)=0.01.*Aj(~here);
%     elseif strcmp('swish',afj)
%         Zj=Aj./(1+exp(-Aj));
%     elseif strcmp('sigmoid',afj)
%         Zj=1./(1+exp(-Aj));
%     elseif strcmp('tanh',afj)
%         Zj=tanh(Aj); 
%         %Zj=LUTtanh(Aj); % < Seems to scale best!
%         %Zj=max(Aj,0);
%         %Zj=Aj.*(945+105.*Aj.^2+Aj.^4)./(945+420.*Aj.^2+15.*Aj.^4);
%         %https://www.math.utah.edu/~beebe/software/ieee/tanh.pdf
%         %rap=@(x) x.*(945+105.*x.^2+x.^4)./(945+420.*x.^2+15.*x.^4);
%     elseif strcmp('linear',afj) % Technically afine, but ok
%         Zj=Aj; % Do nothing
%     end
% 
% end
% 
% ypred=Zj; % Output


