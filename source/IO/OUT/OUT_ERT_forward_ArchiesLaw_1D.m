classdef OUT_ERT_forward_ArchiesLaw_1D < matlab.mixin.Copyable

    properties
        PARA
        CONST
        STATVAR
        TEMP
        TIMESTAMP
        OUTPUT_TIME

    end

    methods

        function obs = provide_PARA(obs)
            obs.PARA.electrode_folder = [];
            obs.PARA.electrode_filename = [];
            obs.PARA.horizontal_grid_spacing = []; %assumed the same in x and z-direction
            %Archies law parameters
            obs.PARA.resistance_water = [];
            obs.PARA.cementation_exponent = []; %m = 1.3 for soils, make dependent
            obs.PARA.tortuosity = []; %a = 0.6-1
            obs.PARA.saturation_exponent = 2; %n

            obs.PARA.output_timestep = [];
            obs.PARA.save_date = [];
            obs.PARA.save_interval = [];
            obs.PARA.tag = [];
        end

        function obs = provide_CONST(obs)

        end

        function obs = provide_STATVAR(obs)

        end

        function obs = finalize_init(obs, tile)
            a = load([obs.PARA.electrode_folder obs.PARA.electrode_filename]);
            obs.STATVAR.electrode_positions = a.electrode_positions;
            %load observations

            forcing = tile.FORCING;

            % Set the next (first) output time. This is the next (first) time output
            % is collected (in memory) for later storage to disk.
            obs.OUTPUT_TIME = forcing.PARA.start_time + obs.PARA.output_timestep;

            % Set the next (first) save time. This is the next (first) time all the
            % collected output is saved to disk.
            if isempty(obs.PARA.save_interval) || isnan(obs.PARA.save_interval)
                obs.SAVE_TIME = forcing.PARA.end_time;
            else
                obs.SAVE_TIME = min(forcing.PARA.end_time,  datenum([obs.PARA.save_date num2str(str2num(datestr(forcing.PARA.start_time,'yyyy')) + obs.PARA.save_interval)], 'dd.mm.yyyy'));
            end
            obs.TEMP.potentials = [];

        end


        function obs = store_OUT(obs, tile)

            t = tile.t;
            forcing = tile.FORCING;
            run_name = tile.PARA.run_name; %tile.RUN_NUMBER;
            result_path = tile.PARA.result_path;
            timestep = tile.timestep;
            out_tag = obs.PARA.tag;

            if t>=obs.OUTPUT_TIME

                disp([datestr(t)])
                obs.TIMESTAMP=[obs.TIMESTAMP t];

                CURRENT = tile.TOP.NEXT;
                while ~is_ground_surface(CURRENT) && ~(strcmp(class(CURRENT), 'Bottom'))
                    CURRENT = CURRENT.NEXT;
                end

                layerThick = [];
                water = [];
                ice = [];
                mineral = [];
                organic = [];

                while ~(strcmp(class(CURRENT), 'Bottom'))
                    layerThick = [layerThick; CURRENT.STATVAR.layerThick];
                    water = [water; CURRENT.STATVAR.water ./CURRENT.STATVAR.layerThick./CURRENT.STATVAR.area];
                    ice = [ice; CURRENT.STATVAR.ice./CURRENT.STATVAR.layerThick./CURRENT.STATVAR.area];
                    mineral = [mineral; CURRENT.STATVAR.mineral./CURRENT.STATVAR.layerThick./CURRENT.STATVAR.area];
                    organic = [organic; CURRENT.STATVAR.organic./CURRENT.STATVAR.layerThick./CURRENT.STATVAR.area];
                    CURRENT = CURRENT.NEXT;
                end

                porosity = 1-ice-mineral-organic;
                saturation = water ./ porosity;
                depths = cumsum([0; layerThick]);
                depths = (depths(1:end-1,1)+depths(2:end,1))./2;
                resistances =  obs.PARA.tortuosity .* obs.PARA.resistance_water .* porosity.^(-obs.PARA.cementation_exponent) .* saturation.^(-obs.PARA.saturation_exponent);

                %calculation of potentials

                %assign electrode-positions
                electrode_positions = obs.STATVAR.electrode_positions;
                el_pos=[];
                for i=1:4
                    el_pos=[el_pos electrode_positions(:,i) electrode_positions(:,i).*0];
                end
                electrode_positions = el_pos;
                electrode_positions_x = electrode_positions(:,1:2:7);

                %define grid
                delta_x = obs.PARA.horizontal_grid_spacing; %in meter
                length_of_survey = max(electrode_positions_x(:)) - min(electrode_positions_x(:));
                midpoint_of_survey = (max(electrode_positions_x(:)) - min(electrode_positions_x(:)))./2 + min(electrode_positions_x(:));
                electrode_positions_shifted = electrode_positions;
                electrode_positions_shifted(:,1:2:7) = electrode_positions_shifted(:,1:2:7)-midpoint_of_survey;

                dx = repmat(delta_x, round(length_of_survey./delta_x)+2, 1); %append some cells on the side
                i=1;
                while 2^i<=64
                    dx = [delta_x .* 2.^i; dx; delta_x .* 2.^i];
                    i=i+1;
                end

                dz = layerThick; %take grid from model, does not seem
                resistances = repmat(resistances, 1, size(dx,1));
                resistances = resistances';

                %current electrodes
                tic;
                [fwd_model_para] = get_2_5Dpara(obs, electrode_positions_shifted(:,1:4),dx,dz,1./resistances,4,electrode_positions_shifted(:,5:8),[1:size(electrode_positions,1)]');
                toc;
                %Run the forward model
                [result, potential_field] = dcfw2_5D(obs, 1./resistances, fwd_model_para);

                %current electrodes first, then potential electrodes
                %potential_field could be saved for diagnosis
                obs.TEMP.potentials = [obs.TEMP.potentials; result'];

                obs.OUTPUT_TIME = min(obs.SAVE_TIME, obs.OUTPUT_TIME + obs.PARA.output_timestep);
                if t>=obs.SAVE_TIME
                    % It is time to save all the collected model output to disk
                     
                    CG_out = obs.TEMP.potentials;

                    if ~(exist([result_path run_name])==7)
                        mkdir([result_path run_name])
                    end
                    if isempty(out_tag) || all(isnan(out_tag))
                        save([result_path run_name '/' run_name '_' datestr(t,'yyyymmdd') '.mat'], 'CG_out')
                    else
                        save([result_path run_name '/' run_name '_' out_tag '_' datestr(t,'yyyymmdd') '.mat'], 'CG_out')
                    end
                    
                    % Clear the out structure
                    obs.TIMESTAMP=[];
                    obs.TEMP.potentials = [];
                    if ~isnan(obs.PARA.save_interval)
                        % If save_interval is defined, uptate SAVE_TIME for next save opertion 
                        obs.SAVE_TIME = min(forcing.PARA.end_time,  datenum([obs.PARA.save_date num2str(str2num(datestr(obs.SAVE_TIME,'yyyy')) + obs.PARA.save_interval)], 'dd.mm.yyyy'));
                        % If save_interval is not defined, we will save at the very end of the model run
                        % and thus do not need to update SAVE_TIME (update would fail because save_interval is nan)
					end
                end
            end
        end


        function [b] = boundary_correction_2_5(obs, dx,dz,s,srcterm,k,g)
            %%function [b] = boundary_correction_2_5(dx,dz,srcterm,k,g);
            %% dx is the x direction cell discretization, a vector length [nx,1];
            %% dz is the z direction cell discretization, a vector length [nz,1];
            %% srcterm is matrix size [Number of RHS, 4] the rows contain the coords
            %of the positive and negative sources [xp zp xn zn] Sources do not need to be
            %located at cell centers. This code assumes that x = 0 is the center of the
            %model space and the z = 0 is the surface, and z is positve down.
            %k is the fourier coeffs we are using in the 2.5D approximation
            %g is the weighting coeffs we are using in the 2.5D approximation

            %%Removes singualrity and decreases boundary effects by modifying source
            %%term based on the anaylitical solution for a homogeneous halfspace
            %%written by Adam Pidlisecky June 2005; Last modified Feb 2006.


            %Initialize a few quantities for later
            nx = length(dx);
            nz = length(dz);
            FOR = zeros(nx*nz,nx*nz);

            %%First Create a grid in realspace

            %build the 2d grid - numbered fromn 0 to maximum extent
            z(1) = 0; for i=1:length(dz); z(i+1) = z(i)+dz(i); end;
            x(1) = 0; for i=1:length(dx); x(i+1) = x(i)+dx(i); end;

            %Center the grid about zero
            x = shiftdim(x) - max(x)/2;

            %Set surface to Z = 0
            z= shiftdim(z);

            %find the cell centers
            xc = x(1:end-1) +dx/2;
            zc = z(1:end-1) +dz/2;

            %Make a couple matrices so we don't have to loop through each location
            %below
            [X,Z] = ndgrid(xc,zc);

            U = zeros(nx*nz,size(srcterm,1));
            %solve for u on this grid using average mref;

            %Now we need to average the conductivity structure
            area = dx(:)*dz(:)';

            savg = area.*s;
            savg = sum(savg(:))./sum(area(:));

            %turn the warning off b/c we know there is a divide by zero, we will fix it
            %later.
            warning('off');

            %loop over all sources
            for i = 1:size(srcterm,1);

                %norm of positive current electrode and 1st potential electrode
                pve1 = ((X - srcterm(i,1)).^2 + (Z - srcterm(i,2)).^2).^0.5;
                %norm of negative current electrode and 1st potential electrode
                nve1 = ((X - srcterm(i,3)).^2 + (Z - srcterm(i,4)).^2).^0.5;
                %norm of imaginary positive current electrode and 1st potential electrode
                pveimag1 = ((X - srcterm(i,1)).^2 + (Z + srcterm(i,2)).^2).^0.5;
                %norm of imaginary negative current electrode and 1st potential electrode
                nveimag1 = ((X - srcterm(i,3)).^2 + (Z + srcterm(i,4)).^2).^0.5;
                U(:,i) = reshape(1/(savg*4*pi)*(1./pve1-1./nve1+1./pveimag1-1./nveimag1),nx*nz,1);

            end
            warning('on');

            %%now check for singularites due to the source being on a node
            for i = 1:size(srcterm,1);
                I  = find(isinf(U(:,i)));

                if max(size(I)) > 0;
                    for j = 1:length(I);
                        [a,c] = ind2sub([nx,nz], I(j));
                        %%Check to see if this a surface electrode
                        if c ==1
                            %if it is average over 3 cells
                            U(I(j),i) = mean(U(sub2ind([nx,nz],a+1,c),i) + U(sub2ind([nx,nz],a,c+1),i)...
                                +U(sub2ind([nx,nz],a-1,c),i));
                        else
                            %otherwise average over 4 cells
                            U(I(j),i) = mean(U(sub2ind([nx,nz],a+1,c),i) + U(sub2ind([nx,nz],a,c+1),i)...
                                +  U(sub2ind([nx,nz],a-1,c),i) +  U(sub2ind([nx,nz],a,c-1),i));
                        end;
                    end;
                end;
            end;

            %now that we have the "true" potentials, we need to crank that through our
            %forward operator so we can create a corrected source term
            D = div2d(obs, dx,dz);
            G = grad2d(obs, dx,dz);

            %% Assemble a homogeneous forward operator
            %put it all together to create operator matrix

            savg = savg*ones(size(s));

            R = massf2d(obs, 1./savg,dx,dz);
            S = spdiags(1./diag(R),0,size(R,1),size(R,2));
            %put it all together to create operator matrix

            Ahomo = -D*S*G;

            %Now we form the operator that will yield our new RHS
            for j = 1:length(k);
                %
                %modify the operator based on the fourier transform
                L = Ahomo +(k(j)^2)*spdiags(mkvc(obs, savg),0,size(Ahomo,1),size(Ahomo,2));

                %Assemble the forwad operator, including the Fourier
                %integration
                FOR = FOR + 0.5*g(j)*inv(L);

            end;


            %now we get bnew by solving bnew = FOR*U;
            b = FOR\U;
            disp('Finished Source correction ');
        end

        function[q] = calcRHS_2_5(obs, dx,dz,src);
            % [q] = calcRHS(dx,dz,src)
            %% dx is the x direction cell discretization, a vector length [nx,1];
            %% dz is the z direction cell discretization, a vector length [nz,1];
            %% src is matrix size [Number of RHS, 4] the rows contain the coords
            %of the positive and negative sources [xp zp xn zn] Sources do not need to be
            %located at cell centers. This code assumes that x = 0 is the center of the
            %model space and the z = 0 is the surface, and z is positve down.

            % this algoritm interpolates the source locations onto the grid defined by
            % dx and dz, then it assembles things into an RHS vector for computation of
            % the potential field in dcfw2dmod
            % Adam Pidlisecky, Aug 2005

            nx = length(dx);
            nz = length(dz);
            %find the area of all the cells by taking the kronecker product
            da = dx(:)*dz(:)';
            da =da(:);

            %%%%%%%% Loop over sources - where sources are a dipole %%%%%%%%%%%
            % allocate space
            nh = nx*nz;
            q = spalloc(nh,size(src,1),size(src,1)*16 );

            fprintf('   Calc RHS (for source)   ');

            %if there is more than one source
            for k=1:size(src,1);
                %interpolate the location of the sources to the nearest cell nodes
                Q= interpmat_N2_5(obs, dx,dz,src(k,1),src(k,2)); %%%%%
                Q = Q - interpmat_N2_5(obs, dx,dz,src(k,3),src(k,4));

                q(:,k) = mkvc(obs, Q)./da;

            end;
            disp('Done ');
        end

        function [xc,zc] = cell_centre2d(obs, dx,dz);
            %%[xc,zc] = cell_centre(dx,dz);
            %%Finds the LH coordsystem cartisien coords of the cell centers;
            %%%Adam Pidlisecky, Modified July 2004

            %calculate the cartesian cell centered grid for inversion
            dx = mkvc(obs, dx);
            dz = mkvc(obs, dz);
            %build the 3d grid - numbered fromn 0 to maximum extent
            z = [0; cumsum(dz)];
            x = [0; cumsum(dx)];

            %Center the grid about zero
            x = shiftdim(x) - max(x)/2;

            %z = shiftdim(z) - max(z)/2;
            %Set surface to Z = 0
            z= shiftdim(z);

            %find the cell centers
            xc = x(1:end-1) +dx/2;
            zc = z(1:end-1) +dz/2;
        end

        function [dobs,U] = dcfw2_5D(obs, s,Para);
            %[U] = dcfw2d(s,Para);
            %Solves the equation (-Div*sigma*grad)*u =q;
            %2.5D forward model using fourier cosine transform in y direction to
            %produce an approximation to 3D current flow, while assuming cthe
            %conductivity structure in the y direction is invarient.

            %% s is the conductivity structure in 2d with dimensions size[nx, nz];
            %% This means imagesc(s') would correspond to a crossection view.

            %% Para is a matlab structure that contains all the information required for
            %% forward modelling. this sturcture is generated using get_2_5Dpara.


            %Code written by Adam Pidlisecky, July 2005; last update Aug 2005;

            dx =Para.dx;
            dz =Para.dz;
            %%make the sigma matrix
            R = massf2d(obs, 1./s,dx,dz);
            S = spdiags(1./diag(R),0,size(R,1),size(R,2));
            %put it all together to create operator matrix
            A = -Para.D*S*Para.G;

            %initialize the solution vector
            U = zeros(size(Para.b));

            %%Enter a loop to solve the forward problem for all fourier coeffs
            for i = 1:length(Para.k);

                %%Now we solve the forward problem

                %modify the operator based on the fourier transform
                L = A +(Para.k(i)^2)*spdiags(mkvc(obs, s),0,size(A,1),size(A,2));

                %now integrate for U;
                U = U+Para.g(i)*(L\(0.5*Para.b));

            end;
            disp('Finished forward calc  ');

            %%see if the Q is around and pick data otherwise return an empty vector
            try
                dobs= Qu(obs, Para.Q,U,Para.srcnum);
            catch
                dobs = [];
            end
        end

        function[D,Dx,Dz] = div2d(obs, dx,dz)
            %[D,Dx,Dz] = div2d(dx,dz)
            %% dx is the x direction cell discretization, a vector length [nx,1];
            %% dz is the z direction cell discretization, a vector length [nz,1];
            %Adam Pidlisecky 2005; based on an implementation of E.Haber 1999 for the
            %3D problem

            dx = shiftdim(dx);
            dz = shiftdim(dz);

            Nx = length(dx)-2;
            Nz = length(dz)-2;

            % Number the phi grid
            np = (Nx+2)*(Nz+2);
            GRDp = reshape(1:1:np,Nx+2,Nz+2);


            % Number the Ax grid
            nax = (Nx+1)*(Nz+2);
            GRDax = reshape(1:1:nax, (Nx+1),(Nz+2));

            % Number the Az grid
            naz = (Nx+2)*(Nz+1);
            GRDaz = reshape(1:1:naz, (Nx+2),(Nz+1));

            % Generates the grid
            ex = ones(Nx+2,1);
            ez = ones(Nz+2,1);

            Dx =kron(dx',ez)';
            Dz = kron(ex',dz)';



            %%%%   Generate d/dx  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            lx = []; jx = []; kx = [];

            % Entries (l,j,k)

            lx = mkvc(obs, GRDp(2:end-1,:));
            jx = mkvc(obs, GRDax(1:end-1,:));
            kx = mkvc(obs, -1./Dx(2:end-1,:));

            % Entries (l+1,j,k)

            lx = [lx;lx];
            jx = [jx;mkvc(obs, GRDax(2:end,:))];
            kx = [kx;-kx];

            % BC at x = 0

            lx = [lx;mkvc(obs, GRDp(1,:))];
            jx = [jx;mkvc(obs, GRDax(1,:))];
            kx = [kx;mkvc(obs, 1./Dx(1,:))];

            % BC at x = end

            lx = [lx;mkvc(obs, GRDp(end,:))];
            jx = [jx;mkvc(obs, GRDax(end,:))];
            kx = [kx;mkvc(obs, -1./Dx(end,:))];


            %%%%   Generate d/dz  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            lz = []; jz = []; kz = [];

            % Entries (l,j,k)

            lz = mkvc(obs, GRDp(:,2:end-1));
            jz = mkvc(obs, GRDaz(:,1:end-1));
            kz = mkvc(obs, -1./Dz(:,2:end-1));

            % Entries (l+1,j,k)

            lz = [lz;lz];
            jz = [jz;mkvc(obs, GRDaz(:,2:end))];
            kz = [kz;-kz];

            % BC on z = 0
            lz = [lz; mkvc(obs, GRDp(:,1))];
            jz = [jz; mkvc(obs, GRDaz(:,1))];
            kz = [kz; mkvc(obs, 1./Dz(:,1))];

            % BC on z = end
            lz = [lz; mkvc(obs, GRDp(:,end))];
            jz = [jz; mkvc(obs, GRDaz(:,end))];
            kz = [kz; mkvc(obs, -1./Dz(:,end))];

            %%%%%%%% Generate the div %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            Dx = sparse(lx,jx,kx,np,nax);
            Dz = sparse(lz,jz,kz,np,naz);

            D = [Dx,Dz];
        end

        function [Para] = get_2_5Dpara(obs, srcloc,dx,dz,BC_cor,num,recloc,srcnum);
            %%function [Para] = get_2_5Dpara(srcloc,dx,dz,BC_cor,num,recloc,srcnum);
            %%Function generates the structure required for calculating the potential field
            %%using dcfw_2_5d.m%


            %%srcloc is matrix size [Number of RHS, 4] the rows contain the coords
            %of the positive and negative sources [xp zp xn zn]
            %%Sources and recievers do not need to be located at cell centers.
            %%This code assumes that x = 0 is the center of the
            %%model space and the z = 0 is the surface, and z is positve down.

            %dx is the x direction cell discretization, a vector length [nx,1];

            %dz is the z direction cell discretization, a vector length [nz,1];

            %Bc_cor enter '[]' for no correction, pass in the s matrix to apply correction.

            %num - a scalar containing the number of fourier coefficents you want to
            %use. Enter '0' to use the default parameters (from Xu et al, 2000)

            %recloc - Enter [] if you do not want to calculate data at receiver
            %locations
            %%for a pole-dipole survey this is an [ndata, 2] matrix
            %%containing the x,z coords of the +ve reciever
            %%for a dipole-dipole survey this is an [ndata, 4] matrix
            %%containing the x,z coords of the +ve and -ve reciever locations

            %srcnum - a vector length(recloc) that lists the source term number that
            %a given reciever corresponds to. e.g. if we have 4 receiver locations and
            %the first 2 receivers correspond to the first source term, and the second
            %receivers correspond to the second source term, then srcnum = [1 1 2 2].


            %%Adam Pidlisecky Created Nov, 2005.

            nx = length(dx);
            nz = length(dz);

            %%Assign grid information
            Para.dx = dx;
            Para.dz = dz;
            Para.nx = nx;
            Para.nz = nz;

            %%Create Diveragnce and Gradient operators once so we don't need to calculate them
            %%again
            Para.D = div2d(obs, dx,dz);
            Para.G = grad2d(obs, dx,dz);

            %%optimize k and g for the given survey geometry.
            if num == 0;
                disp('Using default Fourier Coeffs');
                Para.k = [0.0217102 .2161121 1.0608400 5.0765870];
                Para.g = [0.0463660 0.2365931 1.0382080 5.3648010];
            else
                disp('Optimizing for Fourier Coeffs');
                [k,g,obj,err] = get_k_g_opt(obs, dx,dz,srcloc,num);
                %%Assign the k and g values to Para
                Para.k = k;
                Para.g = g;
            end;


            %%Create the right hand side of the forward modeling equation
            %%See if we are applying the BC correction

            if isempty(BC_cor);
                %no correctoion, so we interpolate the src locations onto the grid
                disp('Interpolating source locations');
                Para.b = calcRHS_2_5(obs, dx,dz,srcloc);

            else
                %Calculate the RHS with a BC coorection applied
                disp('Applying BC/Singularity correction');
                Para.b = boundary_correction_2_5(obs, dx,dz,BC_cor,srcloc,Para.k,Para.g);
            end;



            %%See if we are creating a reciever term
            try
                %%Get the Q matrix for the observation points
                Para.Q = interpmat_N2_5(obs, dx,dz,recloc(:,1),recloc(:,2)); %%%%%

                %% See if it is a dipole survey - if not the other electrode is assumed to
                %% be at infinitey
                try

                    Para.Q = Para.Q - interpmat_N2_5(obs, dx,dz,recloc(:,3),recloc(:,4));

                end;

            end;
            %%Assign srcnumbers (empty vector if no receivers were supplied)
            Para.srcnum=srcnum;
        end

        function [k,g,obj,err] = get_k_g_opt(obs, dx,dz,srcterm,num);
            %function [k,g,obj,err]= get_k_g_opt(dx,dz,srcterm,num);
            %%Function for generating a selection optimized wavenumbers
            %%and weighting coeff's for use with dcfw_2_5d.m.

            %dx is the x direction cell discretization, a vector length [nx,1];
            %dz is the z direction cell discretization, a vector length [nz,1];

            %%srcterm is matrix size [Number of RHS, 4] the rows contain the coords
            %of the positive and negative sources [xp zp xn zn]
            %%Sources and recievers do not need to be located at cell centers.
            %%This code assumes that x = 0 is the center of the
            %%model space and the z = 0 is the surface, and z is positve down.

            %%num - scalar value indicating the number of Fourier coeffs we solve for.

            %%Adam Pidlisecky, Created Nov 2005, Last Modified Feb 2006.


            %set the maximum number of iterations for the optimization routine
            itsmax = 25;

            %Max number of radii to search over
            Max_num =2000;

            %Number of linesearch steps
            lsnum = 10;
            %Line Search parameters
            %lower bound
            ls_low_lim = 0.01;
            %upper bound
            ls_up_lim =1;


            %Define observation distances
            rpos =[]; rneg =[]; rpos_im = []; rneg_im = [];


            %hard wired search radius for determining k and g.
            Xradius = [0.1 0.5 1 5 10 20 30];
            Xradius = [zeros(size(Xradius)) Xradius]';
            Zradius = flipud(Xradius);



            for i = 1:size(srcterm,1);
                Xr = Xradius+srcterm(i,1);
                Zr = Zradius+srcterm(i,2);

                %norm of positive current electrode and 1st potential electrode
                rpost = ((Xr-srcterm(i,1)).^2 + (Zr-srcterm(i,2)).^2).^.5;

                %norm of negative current electrode and 1st potential electrode
                rnegt = ((Xr-srcterm(i,3)).^2 + (Zr-srcterm(i,4)).^2).^.5;

                %norm of imaginary positive current electrode and 1st potential electrode
                rpos_imt = ((Xr-srcterm(i,1)).^2 + (Zr+srcterm(i,2)).^2).^.5;

                %norm of imaginary negative current electrode and 1st potential electrode
                rneg_imt = ((Xr-srcterm(i,3)).^2 + (Zr+srcterm(i,4)).^2).^.5;

                rpos = [rpos; rpost(:)]; rneg = [rneg; rnegt(:)];
                rpos_im = [rpos_im; rpos_imt(:)]; rneg_im = [rneg_im; rneg_imt(:)];
            end;


            %%Now we remove all non-unique radii
            rtot = [rpos rneg rpos_im rneg_im];
            rtot = unique(rtot,'rows');


            %Trim the number of radii down to the size Max_num
            tnum = length(rpos);
            if tnum > Max_num
                rpos=rpos(1:ceil(tnum./Max_num):end);
                rneg=rneg(1:ceil(tnum./Max_num):end);
                rpos_im=rpos_im(1:ceil(tnum./Max_num):end);
                rneg_im=rneg_im(1:ceil(tnum./Max_num):end);
            end

            %initialize a starting guess for k0
            k0 = logspace(-2,0.5,num);

            %Calculate the A matrix
            %Set up a matrix of radii
            rinv = (1./rtot(:,1)-1./rtot(:,2)+1./rtot(:,3)-1./rtot(:,4)).^-1;
            %check for any divide by zeros and remove them
            i = find([ ~isinf(sum(1./rtot,2)+rinv(:))]);
            rtot= rtot(i,:); rinv=rinv(i);


            %Form matrices for computation
            rinv1 = rinv*ones(1,num);
            rpos1 = rtot(:,1)*ones(1,num);
            rneg1 = rtot(:,2)*ones(1,num);
            rpos_im1 = rtot(:,3)*ones(1,num);
            rneg_im1 = rtot(:,4)*ones(1,num);

            %Identity vector
            I = ones(size(rpos1,1),1);
            %K values matrix
            Km = ones(size(rpos1,1),1) *(k0(:))';

            %Calculate the A matrix
            A = rinv1.*real(besselk(0,rpos1.*Km)-besselk(0,rneg1.*Km)+besselk(0,rpos_im1.*Km)-besselk(0,rneg_im1.*Km));

            %%Estimate g for the given K values
            v = A*((A'*A)\(A'*I));
            %Evaluate the objective function for the initial guess
            obj(1) = (1-v)'*(1-v);

            %Start counter and initialize the optimization
            its = 1; %iteration counter
            knew=k0; %updated k vector
            stop =0; %Stopping toggle incase A becomes illconditioned
            reduction = 1; %Variable for ensure sufficent decrease between iterations
            % Optimization terminates if objective function is not
            % reduced by at least 5% at each iteration

            while obj >1e-5 & its <itsmax & stop ==0 & reduction > 0.05;

                %%Create the derivative matrix
                dvdk = zeros(length(v),num);
                for i = 1:num;
                    Ktemp = Km;
                    Ktemp(:,i) = 1.05*(Ktemp(:,i));
                    %form an new A matrix
                    A = rinv1.*real(besselk(0,rpos1.*Ktemp)-besselk(0,rneg1.*Ktemp)...
                        +besselk(0,rpos_im1.*Ktemp)-besselk(0,rneg_im1.*Ktemp));

                    L = A'*A;

                    %%Estimate g for the given K values
                    vT = A*((L)\(A'*I));

                    %Calculate the derivative for the appropriate column
                    dvdk(:,i) = (vT-v)./(Ktemp(:,i)-Km(:,i));
                end;

                %Apply some smallness regularization
                h = dvdk'*(I-v)+1e-8*eye(length(knew))*knew(:);
                dk = (dvdk'*dvdk+1e-8*eye(length(knew)))\h;

                %Perform a line-search to maximize the descent
                for j =[1:lsnum];
                    warning off;
                    ls =linspace(ls_low_lim,ls_up_lim,lsnum);
                    ktemp =knew(:) +ls(j)*dk(:);

                    Km = ones(size(rpos1,1),1) *(ktemp(:))';
                    %Matrix of ones
                    %Calculate the A matrix
                    A = rinv1.*real(besselk(0,rpos1.*Km)-besselk(0,rneg1.*Km)...
                        +besselk(0,rpos_im1.*Km)-besselk(0,rneg_im1.*Km));
                    L = A'*A;
                    %%Estimate g for the given K values

                    v = A*(L\(A'*I));
                    objt = (1-v)'*(1-v);
                    ls_res(j,:) = [objt ls(j)];

                    warning on;
                end;

                %Find the smallest objective function from the line-search
                [b,c] = (min(ls_res(:,1)));

                %Create a new guess for k
                knew =knew(:) +ls(c)*dk(:);
                %eval obj funct
                Km = ones(size(rpos1,1),1) *(knew(:))';

                %Calculate the A matrix
                A = rinv1.*real(besselk(0,rpos1.*Km)-besselk(0,rneg1.*Km)...
                    +besselk(0,rpos_im1.*Km)-besselk(0,rneg_im1.*Km));
                %%Estimate g for the given K values

                v = A*((A'*A)\(A'*I));
                obj(its+1)= (1-v)'*(1-v);
                reduction = obj(its)./obj(its+1)-1;
                its = its+1;
                %Check the conditioning of the matrix
                if rcond(A'*A) < 1e-20;
                    knew =knew(:) -ls(c)*dk(:);
                    stop = 1;
                end;
            end;

            %Get the RMS fit
            err = sqrt(obj./length(rpos));
            %The final k values
            k =abs(knew);
            Km = ones(size(rpos1,1),1) *(k(:))';
            %Reform A to obtian the final g values
            A = rinv1.*real(besselk(0,rpos1.*Km)-besselk(0,rneg1.*Km)...
                +besselk(0,rpos_im1.*Km)-besselk(0,rneg_im1.*Km));
            %Calculate g values
            g = ((A'*A)\(A'*I));
        end

        function[G,Gx,Gz] = grad2d(obs, dx,dz)
            % [G,Gx,Gz] = grad2d(dx,dz)
            %% dx is the x direction cell discretization, a vector length [nx,1];
            %% dz is the z direction cell discretization, a vector length [nz,1];
            %Adam Pidlisecky 2005; based on an implementation of E.Haber 1999 for the
            %3D problem

            dx = shiftdim(dx);
            dz = shiftdim(dz);

            Nx = length(dx)-2;
            Nz = length(dz)-2;

            % Number the phi grid
            np = (Nx+2)*(Nz+2);
            GRDp = reshape(1:1:np,Nx+2,Nz+2);

            % Number the Ax grid
            nax = (Nx+1)*(Nz+2);
            GRDax = reshape(1:1:nax, (Nx+1),(Nz+2));

            % Number the Az grid
            naz = (Nx+2)*(Nz+1);
            GRDaz = reshape(1:1:naz, (Nx+2),(Nz+1));

            %%%%   Generate d/dx  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            lx = []; jx = []; kx = [];

            % Generate grid
            ex = ones(Nx+2,1);
            ez = ones(Nz+2,1);
            Dx =kron(dx',ez)';
            % Entries (l,j,k)

            lx = mkvc(obs, GRDax);
            jx = mkvc(obs, GRDp(1:end-1,:));
            kx = mkvc(obs, -2./(Dx(1:end-1,:) + Dx(2:end,:)));

            % Entries (l+1,j,k)

            lx = [lx; lx];
            jx = [jx;mkvc(obs, GRDp(2:end,:))];
            kx = [kx;-kx];


            %%%%   Generate d/dz  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            lz = []; jz = []; kz = [];
            Dz = kron(ex',dz)';
            % Entries (l,j,k)

            lz = mkvc(obs, GRDaz);
            jz = mkvc(obs, GRDp(:,1:end-1));
            kz = mkvc(obs, -2./(Dz(:,1:end-1) + Dz(:,2:end)));;

            % Entries (l,j,k+1)

            lz = [lz; lz];
            jz = [jz;mkvc(obs, GRDp(:,2:end))];
            kz = [kz; -kz];

            Gx = sparse(lx,jx,kx,nax,np);
            Gz = sparse(lz,jz,kz,naz,np);

            G = [Gx;Gz];
        end

        function [Q] = interpmat_N2_5(obs, dx,dz,xr,zr)
            % [Q] = interpmat(dx,dz,xr,yr,zr)
            %
            % Interpolation program for creating interpoaltion operator
            % dx is the x direction cell discretization, a vector length [nx,1];
            % dz is the z direction cell discretization, a vector length [nz,1];
            % xr x- xlocation of the source, assuming dx is centered about zero
            % zr z-location of the source assuming surface is zero, and z is positive
            % down.
            %Will blow up if reciever is located on the outside cell
            % modified july 2005 - Adam Pidlisecky


            dx = shiftdim(dx);
            dz = shiftdim(dz);

            %build the 3d grid - numbered fromn 0 to maximum extent
            z(1) = 0; for i=1:length(dz); z(i+1) = z(i)+dz(i); end;
            x(1) = 0; for i=1:length(dx); x(i+1) = x(i)+dx(i); end;

            %Center the grid about zero
            x = shiftdim(x) - max(x)/2;

            %Set surface to Z = 0
            z= shiftdim(z);

            %find the cell centers
            xc = x(1:end-1) +dx/2;
            zc = z(1:end-1) +dz/2;

            %take care of surface sources by shifting them to first cell centre
            I = find(zr < min(zc));
            zr(I) = min(zc);

            %call linear interp scheme below

            [Q] = linint(obs, xc,zc,xr,zr);

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function[Q] = linint(obs, x,z,xr,zr)
            %
            % This function does local linear interpolation
            % computed for each receiver point in turn
            %
            % calls mkvc
            %
            % [Q] = linint(x,z,xr,zr)
            % Interpolation matrix
            %


            nx = length(x) ;
            nz = length(z) ;

            np = length(xr);
            Q = sparse(np,nx*nz);

            for i = 1:np,

                % fprintf('Point %d\n',i);

                [dd,im] = min(abs(xr(i)-x));

                if  xr(i) - x(im) >= 0,  % Point on the left
                    ind_x(1) = im;
                    ind_x(2) = im+1;
                elseif  xr(i) - x(im) < 0,  % Point on the right
                    ind_x(1) = im-1;
                    ind_x(2) = im;
                end;
                dx(1) = xr(i) - x(ind_x(1));
                dx(2) = x(ind_x(2)) - xr(i);


                [dd,im] = min(abs(zr(i) - z));
                if  zr(i) -z(im) >= 0,  % Point on the left
                    ind_z(1) = im;
                    ind_z(2) = im+1;
                elseif  zr(i) -z(im) < 0,  % Point on the right
                    ind_z(1) = im-1;
                    ind_z(2) = im;
                end;
                dz(1) = zr(i) - z(ind_z(1));
                dz(2) = z(ind_z(2)) - zr(i);

                Dx =  (x(ind_x(2)) - x(ind_x(1)));
                Dz =  (z(ind_z(2)) - z(ind_z(1)));


                % Build the row for the Q matrix
                v = zeros(nx, nz);

                v( ind_x(1), ind_z(1)) = (1-dx(1)/Dx)*(1-dz(1)/Dz);
                v( ind_x(1), ind_z(1)) = (1-dx(1)/Dx)*(1-dz(1)/Dz);
                v( ind_x(2), ind_z(1)) = (1-dx(2)/Dx)*(1-dz(1)/Dz);
                v( ind_x(2), ind_z(1)) = (1-dx(2)/Dx)*(1-dz(1)/Dz);
                v( ind_x(1), ind_z(2)) = (1-dx(1)/Dx)*(1-dz(2)/Dz);
                v( ind_x(1), ind_z(2)) = (1-dx(1)/Dx)*(1-dz(2)/Dz);
                v( ind_x(2), ind_z(2)) = (1-dx(2)/Dx)*(1-dz(2)/Dz);
                v( ind_x(2), ind_z(2)) = (1-dx(2)/Dx)*(1-dz(2)/Dz);

                % Insert Row into Q matrix
                Q(i,:) = mkvc(obs, v)';

            end
        end

        function[S] = massf2d(obs, s,dx,dz)
            % [S] = massf(s,dx,dz)
            %% s is the conductivity structure in 2d with dimensions size[nx nz];
            %% dx is the x direction cell discretization, a vector length [nx,1];
            %% dz is the z direction cell discretization, a vector length [nz,1];

            %Adam Pidlisecky 2005; based on an implementation of E.Haber 1999 for the
            %3D problem

            dx = shiftdim(dx);
            dz = shiftdim(dz);

            Nx = length(dx)-2;
            Nz = length(dz)-2;

            % Number the Ax grid
            nax = (Nx+1)*(Nz+2);
            GRDax = reshape(1:1:nax, (Nx+1),(Nz+2));


            % Number the Az grid
            naz = (Nx+2)*(Nz+1);
            GRDaz = reshape(1:1:naz, (Nx+2),(Nz+1));

            % Generates the 2D grid
            ex = ones(Nx+2,1);
            ez = ones(Nz+2,1);

            Dx =kron(dx',ez)';
            Dz = kron(ex',dz)';

            dA = Dx.*Dz;

            %%%% Generate x Coeficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            l = 2:Nx+2; k = 1:Nz+2;

            % Avarage rho on x face
            rhof = (dA(l,k).*s(l,k) + dA(l-1,k).*s(l-1,k))/2;

            dVf = (dA(l,k) + dA(l-1,k))/2;

            rhof = rhof./dVf;

            lx = []; jx = []; kx = []; rx = [];

            %% Coef (i,j,k)
            lx = mkvc(obs, GRDax);
            jx = mkvc(obs, GRDax);
            kx = mkvc(obs, rhof);

            Sx = sparse(lx,jx,kx);

            %%%% Generate z Coeficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            l = 1:Nx+2; k = 2:Nz+2;

            % Avarage rho on z face
            rhof = (dA(l,k-1).*s(l,k-1) + dA(l,k).*s(l,k))/2;

            dVf = (dA(l,k-1) + dA(l,k))/2;

            rhof = rhof./dVf;

            lz = []; jz = []; kz = [];

            %% Coef (i,j,k)
            lz = mkvc(obs, GRDaz);
            jz = mkvc(obs, GRDaz);
            kz = mkvc(obs, rhof);

            Sz = sparse(lz,jz,kz);
            %%%% Assemble Matrix  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            Oxz = sparse(nax,naz);

            S = [Sx,  Oxz; ...
                Oxz',  Sz];
        end

        function[v] = mkvc(obs, A)
            % v = mkvc(A)
            %Rearrange a matrix into a vector.
            %Can be substituted for v = v(:)
            v = reshape(A,size(A,1)*size(A,2)*size(A,3),1);
        end

        function[v] = Qu(obs, Q,u,srcnum)
            % [v] = Qu(Q,u,srcnum)
            %Selects a subset of data from the entire potential field.
            %Adam Pidlisecky, modified November 2005

            v = [];

            for i = 1:size(u,2)
                %find q cells related to the source config
                j = find(srcnum == i)  ;
                vv = Q(j,:)*u(:,i);
                v = [v;vv];
            end

        end

    end
end

