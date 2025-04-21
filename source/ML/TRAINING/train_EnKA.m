classdef train_EnKA < matlab.mixin.Copyable

    properties
        PARA
        STATVAR
        CONST
        TEMP
    end

    methods
        function train = provide_PARA(train)
            train.PARA.training_fraction = [];
            train.PARA.number_of_iterations = [];
            train.PARA.relative_error_term = 1e-3;
        end

        function train = provide_CONST(train)

        end
        function train = provide_STATVAR(train)

        end

        function train = finalize_init(train, tile)
            sigy = train.PARA.relative_error_term ./ tile.STATVAR.out_std(1, tile.TEMP.var_ID);
            sigy=repmat(sigy, size(tile.STATVAR.out,1), 1);
            train.TEMP.R = sigy.^2;
        end

        function ml = train_ML(train, tile)
            ml = tile.ML;
            for i=1:train.PARA.number_of_iterations
                training_IDs = randsample([1:size(tile.STATVAR.out,1)]', round(train.PARA.training_fraction .* size(tile.STATVAR.out,1)));
                [predicted_ensemble_ML, ml_parameters] = progapagate_ML(tile.ML, tile.STATVAR.in(training_IDs,:));

                new_parameters = EnKA(train, ml_parameters, tile.STATVAR.out(training_IDs,:), predicted_ensemble_ML, train.PARA.number_of_iterations, train.TEMP.R(training_IDs,:), 0);

                tile.ML = reset_parameters(tile.ML, new_parameters);

            end

        end



        function Xu=EnKA(train, X,y,Yp,alpha,R,dostoch)
            %% Xu=EnKA(X,y,Yp,alpha,R,dostoch)
            % Implementation of the ensemble Kalman analysis step.
            % Based on Evensen et al. (2022) and Emerick (2016).

            No=size(y,1);
            Ns=size(X,1);
            Ne=size(X,2);

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

            Xa=X*Pi;
            Ypa=Yp*Pi;


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
                K=C_XY/(C_YY+((Ne-1)*alpha).*R);
            end

            if dostoch
                perts=sqrt(alpha).*sqrt(R)*randn(No,Ne);
                Y=y+perts;
                Xu=X+K*(Y-Yp);
            else
                Innm=y-mean(Yp,2);
                Xm=mean(X,2);
                Xum=Xm+K*Innm;
                Xua=Xa-0.5.*(K*Ypa);
                Xu=Xum+Xua; % Strange bias?
            end
       end

    end

end