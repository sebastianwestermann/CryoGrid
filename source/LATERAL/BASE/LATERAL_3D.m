%========================================================================
% CryoGrid LATERAL_3D class which manages all LATERAL_IA classes for multi tile (3D) runs 
% LATERAL_IA classes that can be used for multi tile runs are stored in the folder LAT_3D.  
% The field IA_CLASSES stores all active LATERAL_IA classes in a cell array
% and LATERAL_3D evaluates this list top-down for each interaction timestep.
% S. Westermann, T. Ingeman-Nielsen, J. Scheer, Oct 2020
%========================================================================

classdef LATERAL_3D < matlab.mixin.Copyable
 
    properties
        IA_TIME_INCREMENT
        IA_TIME
        ACTIVE
        IA_CLASSES
        TOP
        BOTTOM
        PARA
        CONST
        STATVAR % the state variables of the realization itself
        STATVAR2ALL %the state variables of the realization itself that should be sent to all other workers, not only the connected dones
        STATVAR_PRIVATE
        ENSEMBLE %the state variables of the other ensemble members (cell aray)
    end
    
    methods
        
        %---initialization---------
        
        
        
        function lateral = provide_PARA(lateral)
            lateral.PARA.hill_slope = []; %0 or 1 
            lateral.PARA.ia_time_increment = [];
        end
        
        
        function lateral = provide_CONST(lateral)
            lateral.CONST.day_sec = [];
            lateral.CONST.c_w = [];
            lateral.CONST.c_i = [];        
        end
        
        
        function lateral = provide_STATVAR(lateral)
%             lateral.STATVAR.depths = [];
%             lateral.STATVAR.water_status = [];
%             lateral.STATVAR.hydraulicConductivity = [];
%             lateral.STATVAR.water_table_elevation = [];
%             lateral.STATVAR.water_available = [];
%             lateral.STATVAR.T_water = [];
        end
       
        
        function lateral = initialize_excel(lateral)
            
        end
        

        
        function lateral = finalize_init(lateral, tile)
            
            lateral.IA_TIME_INCREMENT = lateral.PARA.ia_time_increment;
            
            lateral.IA_TIME = tile.FORCING.PARA.start_time + lateral.IA_TIME_INCREMENT;
            lateral.TOP = tile.TOP;
            lateral.BOTTOM = tile.BOTTOM;
            
            %user-defined in the parameter file
            for i=1:size(tile.PARA.lateral_IA_classes, 1)
                lat_ia_class = tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.lateral_IA_classes{i,1});
                lat_ia_class = lat_ia_class{tile.PARA.lateral_IA_classes_index, 1};
                lateral.IA_CLASSES{i,1} = copy(lat_ia_class);
                lateral.IA_CLASSES{i} = finalize_init(lateral.IA_CLASSES{i}, tile);
                lateral.IA_CLASSES{i}.PARENT = lateral;
            end
            
%             lateral.ENSEMBLE={};
%             lateral.STATVAR.index = 0; %set index to zero for 1D runs
%             lateral.PARA.num_realizations = 1;

            for i=1:size(lateral.IA_CLASSES,1)
                lateral.IA_CLASSES{i} = set_ia_time(lateral.IA_CLASSES{i}, tile.t);
                lateral.IA_CLASSES{i} = set_ACTIVE(lateral.IA_CLASSES{i}, i, tile.t - lateral.IA_TIME_INCREMENT);
            end

            lateral.ENSEMBLE={};
            lateral.ACTIVE = zeros(size(lateral.IA_CLASSES,1),1);
            
%             lateral.PARA.connected = tile.RUN_INFO.PARA.connected;
%             lateral.PARA.contact_length = tile.RUN_INFO.PARA.contact_length;
%             lateral.PARA.distance = tile.RUN_INFO.PARA.distance;
%             lateral.PARA.worker_number = tile.RUN_INFO.PARA.worker_number;
%             lateral.PARA.num_realizations = tile.RUN_INFO.PARA.number_of_tiles;
            
            lateral.PARA.connected = tile.RUN_INFO.SPATIAL.PARA.connected;
            lateral.PARA.contact_length = tile.RUN_INFO.SPATIAL.PARA.contact_length;
            lateral.PARA.distance = tile.RUN_INFO.SPATIAL.PARA.distance;
            
            lateral.PARA.num_realizations = tile.RUN_INFO.SPATIAL.PARA.number_of_tiles;
            
            lateral.PARA.worker_number = tile.RUN_INFO.PARA.worker_number;
            lateral.STATVAR.index = lateral.PARA.worker_number; % this is double, but it is necessary to store the index also in STATVAR
        
            %lateral.PARA.is_active = 1; %can be used to switch off lateral interactions temporarily
        end

 
        function lateral = update_lateral(lateral, tile)
            for i=1:size(lateral.IA_CLASSES,1)
                lateral.IA_CLASSES{i} = set_ia_time(lateral.IA_CLASSES{i}, tile.t);
                lateral.IA_CLASSES{i} = set_ACTIVE(lateral.IA_CLASSES{i}, i, tile.t - lateral.IA_TIME_INCREMENT);
            end
        end
        
    
        
        function lateral = assign_number_of_realizations(lateral, num_realizations)
             lateral.PARA.num_realizations = num_realizations;           
        end
        
%         function lateral = get3d_PARA(lateral) %initializes 3D parameters, move this to Excel, etc. file later
%             lateral.PARA.run_number = [1; 2; 3];
%             lateral.PARA.connected = [0 1 0; 1 0 1; 0 1 0];
%             
%             lateral.PARA.contact_length = [0 24.1240351380764 0; 24.1240351380764 0 41.7840545427251;0 41.7840545427251 0]; %[0 1; 1  0];
%             lateral.PARA.distance = [0 3.549647869859770 0;3.549647869859770 0 2.366431913239846;0 2.366431913239846 0]; %[0 10; 10 0];
%             %lateral.PARA.class_list ={{'LAT3D_WATER_UNCONFINED_AQUIFER'; 'LAT3D_HEAT'; 'LAT3D_SNOW_CROCUS'}; {'LAT3D_WATER_UNCONFINED_AQUIFER'; 'LAT3D_HEAT'; 'LAT3D_SNOW_CROCUS'}; ...
%            %                             {'LAT3D_WATER_UNCONFINED_AQUIFER'; 'LAT3D_WATER_SEEPAGE_FACE'; 'LAT3D_HEAT'; 'LAT3D_SNOW_CROCUS'}};
% % 
% %             lateral.PARA.run_number = [1; 2];
% %             lateral.PARA.connected = [0 1; 1 0];
% %             lateral.PARA.contact_length = [0 1; 1  0];
% %             lateral.PARA.distance = [0 2; 2 0];
% %             lateral.PARA.class_list ={{'LAT3D_WATER'}; {'LAT3D_WATER'}};      
% 
% %                 lateral.PARA.run_number = [2; 1; 1; 3];
% %                 lateral.PARA.connected = [0 1 0 0; 1 0 1 0; 0 1 0 1; 0 0 1 0]; %[0 1 1; 1 0 1; 1 1 0];
% %                 lateral.PARA.contact_length = [0 1 0 0; 1 0 1 0; 0 1 0 1; 0 0 1 0]; % [0 1 1; 1 0 1; 1 1 0];
% %                 lateral.PARA.distance = [0 2 0 0; 2 0 2 0; 0 2 0 2; 0 0 2 0]; %[0 2 5; 2 0 2; 5 2 0];
% %                 lateral.PARA.class_list ={{'LAT3D_WATER'; 'LAT3D_WATER_RESERVOIR2'}; {'LAT3D_WATER'};  {'LAT3D_WATER'}; {'LAT3D_WATER'; 'LAT3D_WATER_SEEPAGE_FACE2'}};
%             %lateral.PARA.class_list ={{'LAT3D_WATER'}; {'LAT3D_WATER'};  {'LAT3D_WATER'}; {'LAT3D_WATER'}};
% %             lateral.PARA.class_list ={{'LAT3D_WATER_UNCONFINED_AQUIFER'; 'LAT3D_HEAT'; 'LAT3D_SNOW_CROCUS'}; ...
% %                                         {'LAT3D_WATER_UNCONFINED_AQUIFER'; 'LAT3D_HEAT'; 'LAT3D_SNOW_CROCUS'; 'LAT3D_WATER_SEEPAGE_FACE'}};             
%             lateral.PARA.central_worker = 2; %index of worker performing global computations (needed for LATERAL_3D_water) - should be set so that it has connections to as many as possible other workers
%         end
        
        function lateral = get_index(lateral)
            if lateral.PARA.num_realizations > 1
                lateral.STATVAR.index = labindex;
            else
                lateral.STATVAR.index = 1;
            end
            %lateral.PARA.num_realizations = numlabs;
        end
        
                
        function run_number = get_run_number(lateral, run_number)
            if lateral.PARA.num_realizations > 1
                run_number = [run_number '_' num2str(lateral.PARA.run_number(lateral.STATVAR.index,1))]; 
            end
                
        end

        function output_number = get_output_number(lateral, run_number)
            if lateral.PARA.num_realizations > 1
                output_number = [run_number '_' num2str(lateral.STATVAR.index,1)]; 
            end 
        end
        
      
        %---time integration----------------
        %five steps: 1. pull for stratigraphy, 2. pack, send, receive and unpack, 3.
        %calculate fluxes/derivatives, 4. push to stratigraphy (prgnostic step), 5.
        %recompute diagnostic step
        function lateral = interact(lateral, tile)
            t=tile.t;
            if t>=lateral.IA_TIME %lateral.PARA.is_active && 
                if sum(lateral.ACTIVE) > 0
                    %disp(t-floor(t))
                    
                    %re-compute elevations of the different classes
                    CURRENT = lateral.BOTTOM.PREVIOUS;
                    elevation = CURRENT.STATVAR.lowerPos;
                    while ~(strcmp(class(CURRENT), 'Top'))
                        CURRENT.STATVAR.lowerPos = elevation;
                        elevation = elevation + sum(CURRENT.STATVAR.layerThick,1);
                        CURRENT.STATVAR.upperPos = elevation;
                        CURRENT = CURRENT.PREVIOUS;
                    end
                    
                    %compute ground surface elevation 
                    CURRENT = lateral.TOP.NEXT;
                    lateral.STATVAR.ground_surface_elevation = [];
                    while ~(strcmp(class(CURRENT), 'Bottom')) && isempty(lateral.STATVAR.ground_surface_elevation)
                        lateral.STATVAR.ground_surface_elevation = get_groundSurfaceElevation(CURRENT);
                        CURRENT = CURRENT.NEXT;
                    end
                    
                    
                    %PULL information from individual stratigraphy classes
                    for i=1:size(lateral.IA_CLASSES,1)
                        if lateral.ACTIVE(i,1)
                            lateral.IA_CLASSES{i} = pull(lateral.IA_CLASSES{i}, tile); %After that, each lateral class has all the info it needs in STATVAR-> assign all the variables to class.PARENT in pull
                        end
                    end
                    
                    
                    labBarrier;
                    %PACK AND SEND -  only send a single data package to all connected workers
                    if lateral.PARA.num_realizations > 1
                        data_package_out = pack(lateral, 'STATVAR');
                        if ~isempty(lateral.STATVAR2ALL)
                            lateral.STATVAR2ALL.index=lateral.STATVAR.index;
                            data_package_out_all = pack(lateral, 'STATVAR2ALL');
                        else
                            data_package_out_all = [];
                        end
                        %add the STATVAR2ALL info
                        for i = 1:lateral.PARA.num_realizations
                            if lateral.PARA.connected(lateral.STATVAR.index, i)
                                labSend([data_package_out_all; data_package_out], i, 1);
                            elseif i~=lateral.STATVAR.index
                                labSend(data_package_out_all, i, 1);
                            end
                        end
                        for i = 1:lateral.PARA.num_realizations
                            if lateral.PARA.connected(lateral.STATVAR.index, i) || i~=lateral.STATVAR.index
                                data_package_in = labReceive(i, 1);
                                if ~isempty(data_package_in)
                                    lateral = unpack(lateral, data_package_in); %read received column vector and transform into STATVAR
                                end
                            end
                        end
                    end
                    labBarrier;
                    
                    
                    %calculate all derivatives/fluxes
                    for i=1:size(lateral.IA_CLASSES,1)
                        if lateral.ACTIVE(i,1)
                            lateral.IA_CLASSES{i} = get_derivatives(lateral.IA_CLASSES{i}, tile);
                        end
                    end
                    
                    %PUSH information (fluxes) back to individual stratigraphy classes
                    for i=1:size(lateral.IA_CLASSES,1)
                        if lateral.ACTIVE(i,1)
                            lateral.IA_CLASSES{i} = push(lateral.IA_CLASSES{i}, tile);
                        end
                    end
                    
                    CURRENT = lateral.BOTTOM.PREVIOUS;
                    while ~(strcmp(class(CURRENT), 'Top'))
                        CURRENT = compute_diagnostic(CURRENT, tile);
                        CURRENT = check_trigger(CURRENT, tile);
                        CURRENT = CURRENT.PREVIOUS;
                    end

                end
              
               %set ACTIVE for next timestep
               for i=1:size(lateral.IA_CLASSES,1)
                   lateral.IA_CLASSES{i} = set_ACTIVE(lateral.IA_CLASSES{i}, i, t);
               end
               
               lateral.ENSEMBLE ={};
               lateral.STATVAR2ALL = [];
               lateral.IA_TIME = lateral.IA_TIME + lateral.IA_TIME_INCREMENT;

            end
        end
        
        
        %service functions
        function data_package = pack(lateral, VAR) %transform lateral.STATVAR into column vector ready to send
            variables = fieldnames(lateral.(VAR));
            data_package = [];
            for i=1:size(variables,1)
                %variables{i,1}
                data_package=[data_package; size(variables{i,1},2); double(variables{i,1})']; % # of characters followed by characters as doubles
                data_package=[data_package; size(lateral.(VAR).(variables{i,1}),1); lateral.(VAR).(variables{i,1})]; % # of entries followed by values
            end
        end
        
        function lateral = unpack(lateral, data_package) %read received column vector and transform into STATVAR
            i=1;
            while i<=size(data_package,1)
               variable_name = char(data_package(i+1:i+data_package(i,1),1)');
               i = i + data_package(i,1) + 1;
               STATVAR.(variable_name) = data_package(i+1:i+data_package(i,1),1);
               i = i + data_package(i,1) + 1;
            end
            lateral.ENSEMBLE{size(lateral.ENSEMBLE,1)+1,1} = STATVAR;
        end
        
        function lateral = get_overlap_cells(lateral, variable, variable_out) %no need to loop through stratigraphy, al the information should be in lateral
            for i=1:size(lateral.ENSEMBLE,1)
                if lateral.PARA.connected(lateral.STATVAR.index, lateral.ENSEMBLE{i,1}.index)
                    cell_1 = -lateral.STATVAR.(variable);
                    cell_2 = -lateral.ENSEMBLE{i,1}.(variable);
                    
                    lateral.ENSEMBLE{i,1}.(variable_out) = [];%zeros(size(cell_1,1)-1,size(cell_2,1)-1);
                    
                    if size(cell_1,1) > 1 && size(cell_2,1) > 1
                        for i1=1:size(cell_1,1)-1
                            i2=1;
                            a = max(0, - max(cell_1(i1,1), cell_2(i2,1)) + min(cell_1(i1+1,1), cell_2(i2+1,1)));
                            
                            while a <= 0 && i2 < size(cell_2,1)-1
                                %overlap2(i1,i2) = a;
                                i2 = i2+1;
                                a = max(0, - max(cell_1(i1,1), cell_2(i2,1)) + min(cell_1(i1+1,1), cell_2(i2+1,1)));
                            end
                            if a>0
                                lateral.ENSEMBLE{i,1}.(variable_out) = [lateral.ENSEMBLE{i,1}.(variable_out);  [i1  i2 a]];
                            end
                            
                            i2_start = i2;
                            while a > 0 && i2 < size(cell_2,1)-1
                                
                                i2 = i2+1;
                                a = max(0, - max(cell_1(i1,1), cell_2(i2,1)) + min(cell_1(i1+1,1), cell_2(i2+1,1)));
                                if a>0
                                    lateral.ENSEMBLE{i,1}.(variable_out) = [lateral.ENSEMBLE{i,1}.(variable_out);  [i1  i2 a]];
                                end
                            end
                            i2 = i2_start;
                        end
                    end
                end
            end
        end

        function lateral = get_overlap_cells2(lateral, variable, variable_out) %no need to loop through stratigraphy, all the information should be in lateral
            for i=1:size(lateral.ENSEMBLE,1)
                if lateral.PARA.connected(lateral.STATVAR.index, lateral.ENSEMBLE{i,1}.index)
                    cell_1 = -(lateral.STATVAR.(variable) - double(lateral.PARA.hill_slope) .* lateral.STATVAR.ground_surface_elevation);
                    cell_2 = -(lateral.ENSEMBLE{i,1}.(variable) - double(lateral.PARA.hill_slope) .* lateral.ENSEMBLE{i,1}.ground_surface_elevation);
                    
                    lateral.ENSEMBLE{i,1}.(variable_out) = [];%zeros(size(cell_1,1)-1,size(cell_2,1)-1);
                    
                    if size(cell_1,1) > 1 && size(cell_2,1) > 1
                        for i1=1:size(cell_1,1)-1
                            i2=1;
                            a = max(0, - max(cell_1(i1,1), cell_2(i2,1)) + min(cell_1(i1+1,1), cell_2(i2+1,1)));
                            
                            while a <= 0 && i2 < size(cell_2,1)-1
                                %overlap2(i1,i2) = a;
                                i2 = i2+1;
                                a = max(0, - max(cell_1(i1,1), cell_2(i2,1)) + min(cell_1(i1+1,1), cell_2(i2+1,1)));
                            end
                            if a>0
                                lateral.ENSEMBLE{i,1}.(variable_out) = [lateral.ENSEMBLE{i,1}.(variable_out);  [i1  i2 a]];
                            end
                            
                            i2_start = i2;
                            while a > 0 && i2 < size(cell_2,1)-1
                                
                                i2 = i2+1;
                                a = max(0, - max(cell_1(i1,1), cell_2(i2,1)) + min(cell_1(i1+1,1), cell_2(i2+1,1)));
                                if a>0
                                    lateral.ENSEMBLE{i,1}.(variable_out) = [lateral.ENSEMBLE{i,1}.(variable_out);  [i1  i2 a]];
                                end
                            end
                            i2 = i2_start;
                        end
                    end
                end
            end
        end
        
        function overlap = get_overlap_aquifers(lateral, cell_1, cell_2) %same function as get_overlap_cells, but with flexible input
            cell_1 = -cell_1;
            cell_2 = -cell_2;
            
            overlap = [];
            
            if size(cell_1,1) > 1 && size(cell_2,1) > 1
                for i1=1:size(cell_1,1)-1
                    i2=1;
                    a = max(0, - max(cell_1(i1,1), cell_2(i2,1)) + min(cell_1(i1+1,1), cell_2(i2+1,1)));
                    b = (- max(cell_1(i1,1), cell_2(i2,1)) - min(cell_1(i1+1,1), cell_2(i2+1,1))) ./ 2;
                    
                    while a <= 0 && i2 < size(cell_2,1)-1
                        %overlap2(i1,i2) = a;
                        i2 = i2+1;
                        a = max(0, - max(cell_1(i1,1), cell_2(i2,1)) + min(cell_1(i1+1,1), cell_2(i2+1,1)));
                        b = (- max(cell_1(i1,1), cell_2(i2,1)) - min(cell_1(i1+1,1), cell_2(i2+1,1)))./2;
                    end
                    if a>0
                        overlap = [overlap;  [i1  i2 a b]];
                    end
                    
                    i2_start = i2;
                    while a > 0 && i2 < size(cell_2,1)-1
                        
                        i2 = i2+1;
                        a = max(0, - max(cell_1(i1,1), cell_2(i2,1)) + min(cell_1(i1+1,1), cell_2(i2+1,1)));
                        b = (- max(cell_1(i1,1), cell_2(i2,1)) - min(cell_1(i1+1,1), cell_2(i2+1,1))) ./ 2;
                        if a>0
                            overlap = [overlap;  [i1  i2 a b]];
                        end
                    end
                    i2 = i2_start;
                end
            end
        end
        
        function lateral = get_overlap_cells_mass(lateral)
            for i=1:size(lateral.ENSEMBLE,1)
                if lateral.PARA.connected(lateral.STATVAR.index, lateral.ENSEMBLE{i,1}.index)
                    variable_out = 'overlap_mass';
                    
                    thickness1 = sum(lateral.STATVAR.layerThick);
                    thickness2 = sum(lateral.ENSEMBLE{i,1}.layerThick);
                    
                    %make a unified layerThick data set that is then used to get
                    %cell overlap: mean glacier thickness = (g1+g2)/2. s1 =
                    %(g1+g2)/(2 g1); sf2 = (g1+g2)/(2 g2)

                    cell_1 = (thickness1 + thickness2) ./ 2 ./ thickness1 .* cumsum([0; lateral.STATVAR.layerThick]);
                    cell_2 = (thickness1 + thickness2) ./ 2 ./ thickness2 .* cumsum([0; lateral.ENSEMBLE{i,1}.layerThick]);
                    
                    lateral.ENSEMBLE{i,1}.(variable_out) = [];
                    
                    if size(cell_1,1) > 1 && size(cell_2,1) > 1
                        for i1=1:size(cell_1,1)-1
                            i2=1;
                            a = max(0, - max(cell_1(i1,1), cell_2(i2,1)) + min(cell_1(i1+1,1), cell_2(i2+1,1)));
                            
                            while a <= 0 && i2 < size(cell_2,1)-1
                                i2 = i2+1;
                                a = max(0, - max(cell_1(i1,1), cell_2(i2,1)) + min(cell_1(i1+1,1), cell_2(i2+1,1)));
                            end
                            if a>0
                                lateral.ENSEMBLE{i,1}.(variable_out) = [lateral.ENSEMBLE{i,1}.(variable_out);  [i1  i2 a]];
                            end
                            
                            i2_start = i2;
                            while a > 0 && i2 < size(cell_2,1)-1
                                
                                i2 = i2+1;
                                a = max(0, - max(cell_1(i1,1), cell_2(i2,1)) + min(cell_1(i1+1,1), cell_2(i2+1,1)));
                                if a>0
                                    lateral.ENSEMBLE{i,1}.(variable_out) = [lateral.ENSEMBLE{i,1}.(variable_out);  [i1  i2 a]];
                                end
                            end
                            i2 = i2_start;
                        end
                    end
                end
            end
        end
        
        
        %-------------param file generation-----
        function lateral = param_file_info(lateral)
             lateral = provide_PARA(lateral);
             
             lateral.PARA.class_category = 'LATERAL';
             lateral.PARA.STATVAR = [];
             lateral.PARA.default_value.hill_slope = {1};
             lateral.PARA.comment.hill_slope = {'1: hillslope flow, cells connected with ground surface as reference; 0: cells connected with absolute elevation as reference'};
             lateral.PARA.default_value.ia_time_increment = {0.25};
             lateral.PARA.comment.ia_time_increment = {'minimum of constant timestep for lateral interaction classes, LATERAL_IA classes must have multiples of this one [day]'};
             lateral.PARA.options = [];
        end
        
    end
end

