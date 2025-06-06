%========================================================================
% CryoGrid SPATIAL_REFERENCE class UTM
% Defines the locations of the target points as a regular grid in 
% geographic coordinates, with fixed resolution. 
% The region of interest can be selected by optional
% mask classes. For each point, data (e.g. calculated from DEM, or
% other sources) are read by data provder classes. Based on these data, 
% data mask classes can be used to further constrain the modeled domain. 
%
% S. Westermann, Dec 2022
%========================================================================


classdef SPATIAL_FROM_STORE <  matlab.mixin.Copyable

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
            
            proj.PARA.mask_class = []; %acts on the entire 2d matirx
            proj.PARA.mask_class_index = [];
            
            proj.PARA.data_class = [];
            proj.PARA.data_class_index = [];
            
            proj.PARA.data_mask_class = []; %
            proj.PARA.data_mask_class_index = [];
            
            proj.PARA.assign_tile_properties_class = [];
            proj.PARA.assign_tile_properties_class_index = [];

            proj.PARA.new_reference = 0;
        end
        
        function proj = provide_STATVAR(proj)

        end
        
        function proj = provide_CONST(proj)
            
        end
        
        function proj = finalize_init(proj)
           
            
            %apply masks before data sets
            proj.STATVAR.mask = logical(proj.STATVAR.longitude.*1);
            for i=1:size(proj.PARA.mask_class_index,1)
                mask_class = copy(proj.RUN_INFO.PPROVIDER.CLASSES.(proj.PARA.mask_class{i,1}){proj.PARA.mask_class_index(i,1),1});
                mask_class.PARENT = proj;
                mask_class = finalize_init(mask_class);
                mask_class = apply_mask(mask_class); %can be additive or subtractive
            end
            
            %reduce the list to the ones inside the masks
            mask = proj.STATVAR.mask;
            fn = fieldnames(proj.STATVAR);
            for i=1:size(fn,1)  
                    proj.STATVAR.(fn{i,1})(~mask) = [];
            end
            

            %load data sets
            for i=1:size(proj.PARA.data_class,1)
                data_class = copy(proj.RUN_INFO.PPROVIDER.CLASSES.(proj.PARA.data_class{i,1}){proj.PARA.data_class_index(i,1),1});
                data_class.PARENT = proj;
                data_class = finalize_init(data_class);
                data_class = load_data(data_class); %can be additive or subtractive
            end
            
            for i=1:size(proj.PARA.data_mask_class_index,1)
                mask_class = copy(proj.RUN_INFO.PPROVIDER.CLASSES.(proj.PARA.data_mask_class{i,1}){proj.PARA.data_mask_class_index(i,1),1});
                mask_class.PARENT = proj;
                mask_class = finalize_init(mask_class);
                mask_class = apply_mask(mask_class); %can be additive or subtractive
            end
            
            mask = proj.STATVAR.mask;
            fn = fieldnames(proj.STATVAR);
            for i=1:size(fn,1)  
                    proj.STATVAR.(fn{i,1})(~mask) = [];
            end
            
            for i=1:size(proj.PARA.assign_tile_properties_class,1)
                proj.ACTION{i,1} = copy(proj.RUN_INFO.PPROVIDER.CLASSES.(proj.PARA.assign_tile_properties_class{i,1}){proj.PARA.assign_tile_properties_class_index(i,1),1});
                proj.ACTION{i,1} = finalize_init(proj.ACTION{i,1});
                proj.ACTION{i,1}.PROJ = proj;
            end

        end
        
 
        %-------------param file generation-----
        function proj = param_file_info(proj)
            proj = provide_PARA(proj);
            
            proj.PARA.STATVAR = [];
            proj.PARA.class_category = 'SPATIAL_REFERENCE';
            proj.PARA.default_value = [];
            
            proj.PARA.comment.max_lat = {'maximum latitude of model domain'}; 
            
            proj.PARA.comment.min_lat = {'minimum latitude of model domain'};
            
            proj.PARA.comment.lat_grid_cell_size = {'latitudinal spacing of grid'};
            
            proj.PARA.comment.max_lon = {'maximum longitude of model domain'};
            
            proj.PARA.comment.min_lon = {'minimum longitude of model domain'};
            
            proj.PARA.comment.lon_grid_cell_size = {'longitudinal spacing of grid'};
            
            proj.PARA.comment.mask_class = {'list of mask classes, constains the region of interest based on the coordinates'};
            proj.PARA.options.mask_class.name = 'H_LIST';
            proj.PARA.options.mask_class_index.name = 'H_LIST';
            
            proj.PARA.comment.data_class = {'list of data provider classes, provide data for each target location'};
            proj.PARA.options.data_class.name = 'H_LIST';
            proj.PARA.options.data_class_index.name = 'H_LIST';
            
            proj.PARA.comment.data_mask_class = {'list of data mask classes, constrains the region of interest based on the data provided by data provider classes'};
            proj.PARA.options.data_mask_class.name = 'H_LIST';
            proj.PARA.options.data_mask_class_index.name = 'H_LIST';
            
            proj.PARA.comment.assign_tile_properties_class = {'translates the data to changes to the classes used in the simulations for each point, basically cutomizing the "parameter file" for each taret location'};
            proj.PARA.options.assign_tile_properties_class.name = 'H_LIST';
            proj.PARA.options.assign_tile_properties_class.entries_x = {'update_one2one' 'tag_out_w_run_number'};            
        end
        
    end
end

