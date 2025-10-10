classdef train_EnKA3 < matlab.mixin.Copyable

    properties
        PARA
        STATVAR
        CONST
        TEMP
    end

    methods
        function train = provide_PARA(train)
            train.PARA.minibatch_size = [];
            train.PARA.number_of_iterations = [];
            train.PARA.relative_error_term = 1e-5;
        end

        function train = provide_CONST(train)

        end
        function train = provide_STATVAR(train)

        end

         function train = finalize_init(train, tile)

         end


        function ml = train_ML(train, tile)
            ml = tile.ML;

            input_ML = realWorld2NNinput(tile.TRANSFORM_ML_in, tile.STATVAR.in, tile); %only normalizes features
            
            target = realWorld2training(tile.TRANSFORM_ML_out, tile.STATVAR.out, tile);  %compute the physics constraints for target and normalize
            variance = realWorld2training_variance(tile.TRANSFORM_ML_out, train.PARA.relative_error_term, tile);

            for i=1:train.PARA.number_of_iterations
                
                [predicted_ensemble_ML, ml_parameters] = progapagate_ML(tile.ML, input_ML); %get all of the output, that should not be supermuch time
                predicted_ensemble_processed = [];
                for ii=1:size(predicted_ensemble_ML,2)
                    pe_ML =  NNoutput2realWorld(tile.TRANSFORM_ML_out, predicted_ensemble_ML(:,ii), tile);  %transform to real space, then to  target space, like calculating derivatives 
                    predicted_ensemble_processed = [predicted_ensemble_processed realWorld2training(tile.TRANSFORM_ML_out, pe_ML, tile)];
                end
                %reduce to a subspace for training
                training_IDs = randsample([1:size(target,1)]', train.PARA.minibatch_size); %this could also be all of the observations
                new_parameters = EnKA(train, ml_parameters, target(training_IDs,1), predicted_ensemble_processed(training_IDs,:), train.PARA.number_of_iterations, variance(training_IDs,:), 0);

                tile.ML = reset_parameters(tile.ML, new_parameters);
            end
        end


        function ml = train_ML2(train, ml, tile)

            for i=1:train.PARA.number_of_iterations
             %   training_IDs = randsample([1:size(tile.STATVAR.out,1)]', round(train.PARA.training_fraction .* size(tile.STATVAR.out,1)));
                %[predicted_ensemble_ML, ml_parameters] = progapagate_ML(ml, tile.STATVAR.in(training_IDs,:));
                
                [predicted_ensemble_ML, ml_parameters] = progapagate_ML(ml, tile.STATVAR.in);
                new_parameters = EnKA(train, ml_parameters, tile.STATVAR.out, predicted_ensemble_ML, train.PARA.number_of_iterations, train.TEMP.R, 0);

                ml = reset_parameters(ml, new_parameters);

            end

        end

        function Xu=EnKA(train, X,y,Yp,alpha,R,dostoch)
            %% Xu=EnKA(X,y,Yp,alpha,R,dostoch)
            % Implementation of the ensemble Kalman analysis step.
            % Based on Evensen et al. (2022) and Emerick (2016).

            No=size(y,1);
            Ns=size(X,1);
            Ne=size(X,2); %corresponds to number of ensemble members

            % Convenient matrix and vector operations
            INe=eye(Ne);
            ONe=ones(Ne,1);
            Pi=(INe-ONe*ONe'./Ne);


            if numel(R)==No
                R=diag(R);
            elseif numel(R)==1
                R=R.*eye(No);
            elseif size(R,1)~=No||size(R,2)~=No
                error('Check the dimensions of R');
            end

            Xa=X*Pi; %this works, both just compute the deviation of each entry from the ensemble mean
            Ypa=Yp*Pi; %this works only if Yp has same horizontal dimension as X, so not if there are several outputs


            Ypat=Ypa';
            C_XY=Xa*Ypat;
            C_YY=Ypa*Ypat;
            %aR=(Ne*alpha)*R;

            if No>Ne
                % Subspace inversion using truncated SVD.
                % Decomposition R=S*Rhat*S
                S=diag(diag(R).^(1/2)); % Assuming diagonal R
                Sinv=diag(diag(S).^(-1));
                Rhat=eye(No);
                [U,W,~]=svd(Sinv*Ypa,'econ'); % input is No x Ne
                svr=cumsum(diag(W))./sum(diag(W));
                Nn=1:Ne; % Or Ne?
                svthresh=0.99;
                try
                    Nr=min(Nn(svr>svthresh));
                catch
                    disp('hi')
                end
                these=1:Nr;
                Ur=U(:,these);
                Wr=W(these,these);
                Wrinv=diag(diag(Wr.^(-1)));
                Q=(alpha*(Ne-1))*(Wrinv*Ur')*Rhat*(Ur*Wrinv);
                [Z,H,~]=svd(Q,'econ'); % Q is "R" in Emerick's paper.
                L=(Sinv*Ur*Wrinv*Z);
                IHinv=diag((diag(H)+1).^(-1)); % H is diagonal
                Cinv=L*IHinv*L';
                K=C_XY*Cinv;
            else
                K=C_XY/(C_YY+((Ne-1)*alpha).*R); %could be done for each of the dates, either as 3rd dim, or through a loop -> error increaes with each round, on the order of 10^4.*10^-3*std  at the end 
            end

            if dostoch
                perts=sqrt(alpha).*sqrt(R)*randn(No,Ne);
                Y=y+perts;
                Xu=X+K*(Y-Yp);
            else
                Innm=y-mean(Yp,2);
                Xm=mean(X,2);
                Xum=Xm+K*Innm;
                Xua=Xa-0.5.*(K*Ypa); %I don't understand that, observations do not seem to come in - > comes out as 0.5 Xa when error term is small, pulls X closer to the mean each time, becomes (1/(big number * error)* Xa *Ya*Ya' 
                Xu=Xum+Xua; % Strange bias?
            end
       end

    end

end