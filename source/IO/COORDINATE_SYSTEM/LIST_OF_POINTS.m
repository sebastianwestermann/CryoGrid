%========================================================================
% CryoGrid SPATIAL_REFERENCE class LIST_OF_POINTS
% S. Westermann, Dec 2022
%========================================================================


classdef LIST_OF_POINTS < matlab.mixin.Copyable

    properties
        RUN_INFO
        PARA
        CONST
        STATVAR
        TEMP
        ACTION
    end
    
    methods
        function proj = provide_PARA(proj)
            proj.PARA.info_file_folder = []; %first row contains variable names
            proj.PARA.info_file_name = []; 
            
            proj.PARA.assign_tile_properties_class = [];
            proj.PARA.assign_tile_properties_class_index = [];

            proj.PARA.new_reference = 1;
        end
        
        function proj = provide_STATVAR(proj)

        end
        
        function proj = provide_CONST(proj)
            
        end
        
        function proj = finalize_init(proj)
            
            [~,~,raw] = xlsread([proj.PARA.info_file_folder proj.PARA.info_file_name]);
            
            for i=1:size(raw,2)
                if strcmp(raw{2,i}, 'num')
                    proj.STATVAR.(raw{1,i}) = cell2mat(raw(3:end,i));
                else
                    proj.STATVAR.(raw{1,i}) = raw(3:end,i);
                end
            end
            
            variables = fieldnames(proj.STATVAR);
            lat_yes=0;
            lon_yes=0;
            alt_yes=0;
            for i=1:size(variables,1)
                if strcmp(variables{i,1}, 'latitude')
                    lat_yes = 1;
                end
                if strcmp(variables{i,1}, 'longitude')
                    lon_yes = 1;
                end
                if strcmp(variables{i,1}, 'altitude')
                    alt_yes = 1;
                end
            end
            if ~lat_yes
                proj.STATVAR.latitude = repmat(70, size(proj.STATVAR.(variables{1,1}),1),1);
            end
            if ~lon_yes
                proj.STATVAR.longitude = repmat(70, size(proj.STATVAR.(variables{1,1}),1),1);
            end

            if ~alt_yes
                proj.STATVAR.altitude = repmat(70, size(proj.STATVAR.(variables{1,1}),1),1);
            end
            proj.STATVAR.key = [1:size(proj.STATVAR.latitude,1)]';


            % Calculate skyview_factor, originally code from POINT_SLOPE.
            % NOTE: proj.STATVAR are cells, not floats
            proj.STATVAR.horizon_bins = repmat([0],size(proj.STATVAR.latitude,1),1); % [0, 360]; %
            proj.STATVAR.horizon_angles = repmat([0],size(proj.STATVAR.latitude,1),1); % [0, 0]; %
 
            if isempty(proj.STATVAR.skyview_factor) || sum(isnan(proj.STATVAR.skyview_factor))>0
                
                for i = 1:length(proj.STATVAR.latitude)
                    bin_increment = 360/proj.STATVAR.number_of_horizon_bins(i);
                    horizon_bins = 0:bin_increment:360;
                    horizon_angles = zeros(size(horizon_bins));
                % horizon_bins = proj.STATVAR.horizon_bins(i,:); %[0:360]';
                % horizon_angles = interp1(proj.STATVAR.horizon_bins(i,:), proj.STATVAR.horizon_angles(i,:), horizon_bins);
                
                azmRadian = (pi/180).*horizon_bins;
     
                % convert output from horizon program to radians and translate to angle
                % from zenith
                H = (pi/180).*(90 - horizon_angles);
                
                aspectRadian = (pi/180).*(proj.STATVAR.aspect);
                % modify limits of integration for slopes facing away from horizons
              
                t = cosd(proj.STATVAR.aspect(i)-proj.STATVAR.horizon_bins(i))<0;
                %Simplified trig, the original was H(t) = min(H(t),...
                %  acos(-cos(azmRadian(t)-aspectRadian)*sind(slopeDegrees)./...
                %  sqrt(cosd(slopeDegrees)^2+sind(slopeDegrees)^2*cos(azmRadian(t)-aspectRadian).^2)));
                % but same as
                H(t) = min(H(t), acos(sqrt(1-1./(1+tand(proj.STATVAR.slope_angle(i)).^2.*cos(azmRadian(t)-aspectRadian(i)).^2))));
                    qIntegrand = (cosd(proj.STATVAR.slope_angle(i)).*sin(H).^2 + sind(proj.STATVAR.slope_angle(i)).*cos(aspectRadian(i)-azmRadian).*(H-cos(H).*sin(H)))/2;
              
                % shouldn't be any negative, except perhaps rounding error, so just in case
                    qIntegrand(qIntegrand<0) = 0;
                
                % integrate
                    proj.STATVAR.skyview_factor(i) = trapz(azmRadian,qIntegrand)./pi;
                end
            end
            
%             %apply masks before data sets
%             proj.STATVAR.mask = logical(proj.STATVAR.longitude.*1);
%             for i=1:size(proj.PARA.mask_class_index,1)
%                 mask_class = copy(proj.RUN_INFO.PPROVIDER.CLASSES.(proj.PARA.mask_class{i,1}){proj.PARA.mask_class_index(i,1),1});
%                 mask_class.PARENT = proj;
%                 mask_class = finalize_init(mask_class);
%                 mask_class = apply_mask(mask_class); %can be additive or subtractive
%             end
%             
%             %reduce the list to the ones inside the masks
%             mask = proj.STATVAR.mask;
%             fn = fieldnames(proj.STATVAR);
%             for i=1:size(fn,1)  
%                     proj.STATVAR.(fn{i,1})(~mask) = [];
%             end
%             
% 
%             %load data sets
%             for i=1:size(proj.PARA.data_class,1)
%                 data_class = copy(proj.RUN_INFO.PPROVIDER.CLASSES.(proj.PARA.data_class{i,1}){proj.PARA.data_class_index(i,1),1});
%                 data_class.PARENT = proj;
%                 data_class = finalize_init(data_class);
%                 data_class = load_data(data_class); %can be additive or subtractive
%             end
%             
%             for i=1:size(proj.PARA.data_mask_class_index,1)
%                 mask_class = copy(proj.RUN_INFO.PPROVIDER.CLASSES.(proj.PARA.data_mask_class{i,1}){proj.PARA.data_mask_class_index(i,1),1});
%                 mask_class.PARENT = proj;
%                 mask_class = finalize_init(mask_class);
%                 mask_class = apply_mask(mask_class); %can be additive or subtractive
%             end
%             
%             mask = proj.STATVAR.mask;
%             fn = fieldnames(proj.STATVAR);
%             for i=1:size(fn,1)  
%                     proj.STATVAR.(fn{i,1})(~mask) = [];
%             end
            
            for i=1:size(proj.PARA.assign_tile_properties_class,1)
                proj.ACTION{i,1} = copy(proj.RUN_INFO.PPROVIDER.CLASSES.(proj.PARA.assign_tile_properties_class{i,1}){proj.PARA.assign_tile_properties_class_index(i,1),1});
                proj.ACTION{i,1} = finalize_init(proj.ACTION{i,1});
                proj.ACTION{i,1}.PROJ = proj;
            end

        end
        
 
%         %-------------param file generation-----
%         function proj = param_file_info(proj)
%             proj = provide_PARA(proj);
%             
%             proj.PARA.STATVAR = [];
%             proj.PARA.class_category = 'SPATIAL_REFERENCE';
%             proj.PARA.default_value = [];
%             
%             proj.PARA.comment.max_lat = {'maximum latitude of model domain'}; 
%             
%             proj.PARA.comment.min_lat = {'minimum latitude of model domain'};
%             
%             proj.PARA.comment.lat_grid_cell_size = {'latitudinal spacing of grid'};
%             
%             proj.PARA.comment.max_lon = {'maximum longitude of model domain'};
%             
%             proj.PARA.comment.min_lon = {'minimum longitude of model domain'};
%             
%             proj.PARA.comment.lon_grid_cell_size = {'longitudinal spacing of grid'};
%             
%             proj.PARA.comment.mask_class = {'list of mask classes, constains the region of interest based on the coordinates'};
%             proj.PARA.options.mask_class.name = 'H_LIST';
%             proj.PARA.options.mask_class_index.name = 'H_LIST';
%             
%             proj.PARA.comment.data_class = {'list of data provider classes, provide data for each target location'};
%             proj.PARA.options.data_class.name = 'H_LIST';
%             proj.PARA.options.data_class_index.name = 'H_LIST';
%             
%             proj.PARA.comment.data_mask_class = {'list of data mask classes, constrains the region of interest based on the data provided by data provider classes'};
%             proj.PARA.options.data_mask_class.name = 'H_LIST';
%             proj.PARA.options.data_mask_class_index.name = 'H_LIST';
%             
%             proj.PARA.comment.assign_tile_properties_class = {'translates the data to changes to the classes used in the simulations for each point, basically cutomizing the "parameter file" for each taret location'};
%             proj.PARA.options.assign_tile_properties_class.name = 'H_LIST';
%             proj.PARA.options.assign_tile_properties_class.entries_x = {'update_one2one' 'tag_out_w_run_number'};            
%         end
        
    end
end

