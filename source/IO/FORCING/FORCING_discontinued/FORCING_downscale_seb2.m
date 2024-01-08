%========================================================================
% CryoGrid FORCING class FORCING_downscale_slope_seb2
% same as FORCING_downscale_slope_seb, but assumes that there is no trend 

%
% Authors:
% S. Westermann, T. Ingeman-Nielsen, J. Scheer, October 2020
% T. Ingeman-Nielsen, October 2022
%
%========================================================================

classdef FORCING_downscale_seb2 < FORCING_base & READ_FORCING_mat
    
    methods
        
        function forcing = provide_PARA(forcing)         
            % INITIALIZE_PARA  Initializes PARA structure, setting the variables in PARA.  

            forcing.PARA.carrier_forcing_class = [];  
            forcing.PARA.carrier_forcing_class_index = [];
            forcing.PARA.offset_from_GMT_carrier = []; %in hours
            forcing.PARA.reference_forcing_class = [];   
            forcing.PARA.reference_forcing_class_index = [];
            forcing.PARA.offset_from_GMT_reference = []; %in hours
            
            forcing.PARA.start_overlap = [];
            forcing.PARA.end_overlap = [];
            
            forcing.PARA.variables = [];
            forcing.PARA.overlap_interval = []; %interval for which the list of overalp-pairs is created
           % forcing.PARA.downscaling_function = [];
            
            %all the stuff needed for the "normal" part of the class
            forcing.PARA.start_time = []; % start time of the simulations (must be within the range of data in forcing file)
            forcing.PARA.end_time = [];   % end time of the simulations (must be within the range of data in forcing file)
            forcing.PARA.rain_fraction = [];  %rainfall fraction assumed in sumulations (rainfall from the forcing data file is multiplied by this parameter)
            forcing.PARA.snow_fraction = [];  %snowfall fraction assumed in sumulations (snowfall from the forcing data file is multiplied by this parameter)
            forcing.PARA.albedo_surrounding_terrain = [];
            forcing.PARA.heatFlux_lb = [];  % heat flux at the lower boundary [W/m2] - positive values correspond to energy gain
            forcing.PARA.airT_height = [];  % height above ground at which air temperature (and wind speed!) from the forcing data are applied
            forcing.PARA.all_rain_T = [];
            forcing.PARA.all_snow_T = [];
            
            forcing.PARA.post_proc_class = [];  %optional post-processing classes
            forcing.PARA.post_proc_class_index = [];
        end
        
        
        function forcing = provide_CONST(forcing)
            forcing.CONST.Tmfw = [];
            forcing.CONST.sigma = [];            
        end
        

        function forcing = provide_STATVAR(forcing)
            
        end
                
        
        function forcing = finalize_init(forcing, tile)
            
            forcing = set_start_and_end_time(forcing); % assign start/end time
           
            %load carrier and reference
            carrier_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(forcing.PARA.carrier_forcing_class){forcing.PARA.carrier_forcing_class_index,1});
            carrier_class = finalize_init(carrier_class, tile);
            carrier_class.DATA.timeForcing = carrier_class.DATA.timeForcing - forcing.PARA.offset_from_GMT_carrier ./ 24;
            reference_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(forcing.PARA.reference_forcing_class){forcing.PARA.reference_forcing_class_index,1});
            reference_class.PARA.start_time = [forcing.PARA.start_overlap(1,1);1; 1];
            reference_class.PARA.end_time = [forcing.PARA.end_overlap(1,1)+1; 1; 1];
            reference_class = finalize_init(reference_class, tile);
            reference_class.DATA.timeForcing = reference_class.DATA.timeForcing - forcing.PARA.offset_from_GMT_reference ./ 24;
            
            if ~isfield(carrier_class.DATA, 'precip')
                carrier_class.DATA.precip = carrier_class.DATA.snowfall + carrier_class.DATA.rainfall;
            end
            if ~isfield(reference_class.DATA, 'precip')
                reference_class.DATA.precip = reference_class.DATA.snowfall + reference_class.DATA.rainfall;
            end
            
            
            %find indices for overlap period
            if ~isempty(forcing.PARA.start_overlap) && sum(isnan(forcing.PARA.start_overlap))==0
                forcing.PARA.start_overlap = datenum(forcing.PARA.start_overlap(1,1), forcing.PARA.start_overlap(2,1), forcing.PARA.start_overlap(3,1));
                start_index_carrier = find(carrier_class.DATA.timeForcing(:,1) >= forcing.PARA.start_overlap - (carrier_class.DATA.timeForcing(2,1)-carrier_class.DATA.timeForcing(1,1))./2,1); %subtract half a timestep
                start_index_reference = find(reference_class.DATA.timeForcing(:,1) >= forcing.PARA.start_overlap - (reference_class.DATA.timeForcing(2,1)-reference_class.DATA.timeForcing(1,1))./2,1); %subtract half a timestep
            else %determine which class starts first
                if carrier_class.DATA.timeForcing(1,1) - (carrier_class.DATA.timeForcing(2,1)-carrier_class.DATA.timeForcing(1,1))./2 > reference_class.DATA.timeForcing(1,1) - (reference_class.DATA.timeForcing(2,1)-reference_class.DATA.timeForcing(1,1))./2
                    start_index_carrier = 1;
                    start_index_reference = find(reference_class.DATA.timeForcing(:,1) >= carrier_class.DATA.timeForcing(1,1) - (carrier_class.DATA.timeForcing(2,1)-carrier_class.DATA.timeForcing(1,1))./2 - ...
                        (reference_class.DATA.timeForcing(2,1)-reference_class.DATA.timeForcing(1,1))./2,1); %subtract half a timestep

                else
                    start_index_carrier = find(carrier_class.DATA.timeForcing(:,1) >= reference_class.DATA.timeForcing(1,1) - (carrier_class.DATA.timeForcing(2,1)-carrier_class.DATA.timeForcing(1,1))./2 - ...
                        (reference_class.DATA.timeForcing(2,1)-reference_class.DATA.timeForcing(1,1))./2,1); %subtract half a timestep
                    start_index_reference = 1;
                end
            end
            if ~isempty(forcing.PARA.end_overlap) && sum(isnan(forcing.PARA.end_overlap))==0
                forcing.PARA.end_overlap = datenum(forcing.PARA.end_overlap(1,1), forcing.PARA.end_overlap(2,1), forcing.PARA.end_overlap(3,1));
                end_index_carrier = find(carrier_class.DATA.timeForcing(:,1) < forcing.PARA.end_overlap - (carrier_class.DATA.timeForcing(2,1)-carrier_class.DATA.timeForcing(1,1))./2,1, 'last'); %subtract half a timestep
                end_index_reference = find(reference_class.DATA.timeForcing(:,1) < forcing.PARA.end_overlap - (reference_class.DATA.timeForcing(2,1)-reference_class.DATA.timeForcing(1,1))./2,1, 'last'); %subtract half a timestep
            else
                if carrier_class.DATA.timeForcing(end,1) + (carrier_class.DATA.timeForcing(2,1)-carrier_class.DATA.timeForcing(1,1))./2 > reference_class.DATA.timeForcing(end,1) + (reference_class.DATA.timeForcing(2,1)-reference_class.DATA.timeForcing(1,1))./2
                    end_index_carrier = find(carrier_class.DATA.timeForcing(:,1) < reference_class.DATA.timeForcing(end,1) + (carrier_class.DATA.timeForcing(2,1)-carrier_class.DATA.timeForcing(1,1))./2 + ...
                        (reference_class.DATA.timeForcing(2,1)-reference_class.DATA.timeForcing(1,1))./2,1, 'last'); %subtract half a timestep
                    end_index_reference = size(reference_class.DATA.timeForcing,1);
                else
                    end_index_carrier = size(carrier_class.DATA.timeForcing,1);
                    end_index_reference = find(reference_class.DATA.timeForcing(:,1) < carrier_class.DATA.timeForcing(end,1) + (carrier_class.DATA.timeForcing(2,1)-carrier_class.DATA.timeForcing(1,1))./2 + ...
                        (reference_class.DATA.timeForcing(2,1)-reference_class.DATA.timeForcing(1,1))./2, 1, 'last'); %subtract half a timestep
                end
            end
            
             
            random_sequence_obtained = 0;
            random_sequence = [];
            pos_q=find(strcmp(forcing.PARA.variables, 'q'));
            if ~isempty(pos_q)
                forcing.PARA.variables(pos_q,:) =[];
                forcing.DATA.q = [];
                get_q=1;
            else
                get_q=0;
            end
            
            mean_value_reference = [];
            
            for i=1:size(forcing.PARA.variables,1)
                if strcmp(forcing.PARA.overlap_interval, 'year')
                
                elseif strcmp(forcing.PARA.overlap_interval, 'month')
                    %compile overlap pairs
                    overlap_pairs_time = []; 
                    overlap_pairs = [];

                    y = year(max(carrier_class.DATA.timeForcing(start_index_carrier,1), reference_class.DATA.timeForcing(start_index_reference,1)));
                    m = month(max(carrier_class.DATA.timeForcing(start_index_carrier,1), reference_class.DATA.timeForcing(start_index_reference,1)));
                    end_y=year(min(carrier_class.DATA.timeForcing(end_index_carrier), reference_class.DATA.timeForcing(end_index_reference)));
                    end_m=month(min(carrier_class.DATA.timeForcing(end_index_carrier), reference_class.DATA.timeForcing(end_index_reference)));
                    while datenum(y,m,1) <= datenum(end_y, end_m, 1)
                        overlap_pairs_time = [overlap_pairs_time; datenum(y,m,15)];
                        pos_carrier = find(carrier_class.DATA.timeForcing(:,1)>=datenum(y,m,1) & carrier_class.DATA.timeForcing(:,1)<datenum(y,m+1,1));
                        pos_reference = find(reference_class.DATA.timeForcing(:,1)>=datenum(y,m,1) & reference_class.DATA.timeForcing(:,1)<datenum(y,m+1,1));
                        overlap_pairs = [overlap_pairs; mean(carrier_class.DATA.(forcing.PARA.variables{i,1})(pos_carrier,1)) ...
                            mean(reference_class.DATA.(forcing.PARA.variables{i,1})(pos_reference,1))];
                        m=m+1;
                        if m==13
                            m=1;
                            y=y+1;
                        end
                    end
                    reference_before_time = [];
                    reference_before = [];
                    y = year(reference_class.DATA.timeForcing(1,1));
                    m = month(reference_class.DATA.timeForcing(1,1));
                    end_y = year(max(carrier_class.DATA.timeForcing(start_index_carrier,1), reference_class.DATA.timeForcing(start_index_reference,1)));
                    end_m = month(max(carrier_class.DATA.timeForcing(start_index_carrier,1), reference_class.DATA.timeForcing(start_index_reference,1)));
                    while datenum(y,m,1) < datenum(end_y, end_m, 1)
                        reference_before_time = [reference_before_time; datenum(y,m,15)];
                        pos_reference = find(reference_class.DATA.timeForcing(:,1)>=datenum(y,m,1) & reference_class.DATA.timeForcing(:,1)<datenum(y,m+1,1));
                        reference_before = [;  reference_before; mean(reference_class.DATA.(forcing.PARA.variables{i,1})(pos_reference,1))];
                        m=m+1;
                        if m==13
                            m=1;
                            y=y+1;
                        end
                    end
                    reference_after_time = []; 
                    reference_after = [];                    
                    y=year(min(carrier_class.DATA.timeForcing(end_index_carrier), reference_class.DATA.timeForcing(end_index_reference)));
                    m=month(min(carrier_class.DATA.timeForcing(end_index_carrier), reference_class.DATA.timeForcing(end_index_reference)));
                    while datenum(y,m,1) < datenum(year(reference_class.DATA.timeForcing(end,1)), month(reference_class.DATA.timeForcing(end,1))+1, 1)
                        reference_after_time = [reference_after_time; datenum(y,m,15)];
                        pos_reference = find(reference_class.DATA.timeForcing(:,1)>=datenum(y,m,1) & reference_class.DATA.timeForcing(:,1)<datenum(y,m+1,1));
                        reference_after = [;  reference_after; mean(reference_class.DATA.(forcing.PARA.variables{i,1})(pos_reference,1))];
                        m=m+1;
                        if m==13
                            m=1;
                            y=y+1;
                        end
                    end
                    

                    mvr = [];
                    
                    for m=1:12
                        pos = find(month(overlap_pairs_time(:,1))==m);
                        mvr = [mvr; mean(overlap_pairs(pos, 2))];
%                         absolute_offset = [absolute_offset; mean(overlap_pairs(pos, 2) - overlap_pairs(pos,1))];
%                         relative_offset = [relative_offset; mean(overlap_pairs(pos, 2)./overlap_pairs(pos,1))];
                    end
                    
                    mean_value_reference = [mean_value_reference mvr];
                    
                elseif strcmp(forcing.PARA.overlap_interval, 'day')
                    
                
                end
            end
            
            %load new reference
            reference_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(forcing.PARA.reference_forcing_class){forcing.PARA.reference_forcing_class_index,1});
            reference_class.PARA.start_time = [year(forcing.PARA.start_time);1; 1];
            reference_class.PARA.end_time = [year(forcing.PARA.end_time)+1; 1; 1];
            reference_class = finalize_init(reference_class, tile);
            reference_class.DATA.timeForcing = reference_class.DATA.timeForcing - forcing.PARA.offset_from_GMT_reference ./ 24;
            
            
            for i=1:size(forcing.PARA.variables,1)
                
                if strcmp(forcing.PARA.overlap_interval, 'month')
                                        
                    forcing.DATA.(forcing.PARA.variables{i,1}) = [];
                    y = year(reference_class.DATA.timeForcing(1,1));
                    m = month(reference_class.DATA.timeForcing(1,1));
                    count = 1;
                    while datenum(y,m,1) < datenum(year(reference_class.DATA.timeForcing(end,1)), month(reference_class.DATA.timeForcing(end,1))+1, 1)
                        
                        %select random month from carrier
                        if random_sequence_obtained == 0
                            pos = find(month(overlap_pairs_time) == m);
                            this_pos = randperm(size(pos,1),1);
                            pos = pos(this_pos,1);
                            random_sequence = [random_sequence; pos];
                        else
                            pos = random_sequence(count,1);
                            count = count+1;
                        end
                        
                        pos_start_carrier = find(carrier_class.DATA.timeForcing(:,1) >= datenum(year(overlap_pairs_time(pos,1)), m, 1), 1);
                        pos_end_carrier = find(carrier_class.DATA.timeForcing(:,1) < datenum(year(overlap_pairs_time(pos,1)), m+1, 1), 1, 'last');

                        new_month = carrier_class.DATA.(forcing.PARA.variables{i,1})(pos_start_carrier:pos_end_carrier,1);
                        
                        range = find(reference_class.DATA.timeForcing>=datenum(y,m,1) & reference_class.DATA.timeForcing<datenum(y,m+1,1));
                        if strcmp(forcing.PARA.variables{i,1},'precip') || strcmp(forcing.PARA.variables{i,1},'Sin') || strcmp(forcing.PARA.variables{i,1},'wind')
                            relative_change = mean(reference_class.DATA.(forcing.PARA.variables{i,1})(range,1)) ./ mean_value_reference(m,i);
                            if isnan(relative_change)
                                relative_change = 1;
                            end
                            new_month = max(0, relative_change .* new_month);
                        else
                            absolute_offset_reference = mean(reference_class.DATA.(forcing.PARA.variables{i,1})(range,1)) - mean_value_reference(m,i);
                            new_month = new_month + absolute_offset_reference;
                         end
                         
                         forcing.DATA.(forcing.PARA.variables{i,1}) = [forcing.DATA.(forcing.PARA.variables{i,1}); new_month];

                        if get_q && strcmp(forcing.PARA.variables{i,1}, 'Tair')
                            new_month_q = carrier_class.DATA.q(pos_start_carrier:pos_end_carrier,1);
                            old_month = carrier_class.DATA.Tair(pos_start_carrier:pos_end_carrier,1);
                            forcing.DATA.q = [forcing.DATA.q; new_month_q .* exp(17.62.*new_month./(243.12+new_month)) ./ exp(17.62.*old_month./(243.12+old_month))];
                        end
                        
                        m=m+1;
                        if m==13
                            m=1;
                            y=y+1;
                        end
                    end
                    random_sequence_obtained = 1;
                    starttime = datenum(year(reference_class.DATA.timeForcing(1,1)), month(reference_class.DATA.timeForcing(1,1)),1);
                    timestamp = carrier_class.DATA.timeForcing(2,1) - carrier_class.DATA.timeForcing(1,1);
                    forcing.DATA.timeForcing = [starttime:timestamp:starttime+timestamp.*(size(forcing.DATA.(forcing.PARA.variables{i,1}),1)-1)]';
                    

                end
            end
            
            
            %forcing = adjust_precip_quantile_mapping(forcing, carrier_class, tile);
            forcing = split_precip_Tair(forcing);

            forcing.DATA.rainfall = forcing.DATA.rainfall.*forcing.PARA.rain_fraction;
            forcing.DATA.snowfall = forcing.DATA.snowfall.*forcing.PARA.snow_fraction;
            
            forcing = check_and_correct(forcing); % Remove known errors
             
            forcing = initialize_TEMP(forcing);
%             forcing = initialize_TEMP_slope(forcing);
%             
%             forcing = reduce_precip_slope(forcing, tile);
%             
%             forcing = convert_accumulated2instantaneous_Sin(forcing, tile);
%             
%             forcing = SolarAzEl(forcing, tile);
%             
%             %make own function?
%             if ~isfield(forcing.DATA, 'S_TOA') %this will not really work in this case, would require paleo-S_TOA from 
%                 mu0=max(sind(forcing.DATA.sunElevation),0); % Trunacte negative values.
%                 sunset=mu0<cosd(89);%(mu0==0); % Sunset switch.
%                 forcing.DATA.S_TOA = 1370.*mu0;
%             end
%             
%             forcing = split_Sin(forcing); % split Sin in dir and dif
%             forcing = terrain_corr_Sin_dif(forcing, tile);
%             forcing = reproject_Sin_dir(forcing, tile);
%             forcing = terrain_shade(forcing, tile);
%             forcing.DATA.Sin = forcing.DATA.Sin_dir + forcing.DATA.Sin_dif;
%             
%             set pressure to mean pressure at corresponding altitude (international
%             altitude formula) if not provided by the forcing time series
            if ~isfield(forcing.DATA, 'p')
                altitude = tile.PARA.altitude;
                forcing.DATA.p=forcing.DATA.Tair.*0 + 1013.25.*100.*(1-0.0065./288.15.*altitude).^5.255;
            end
            
            %optional post-processing with dedicated classes
            if ~isempty(forcing.PARA.post_proc_class) && sum(isnan(forcing.PARA.post_proc_class_index)) == 0
                for i=1:size(forcing.PARA.post_proc_class,1)
                    post_proc_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(forcing.PARA.post_proc_class{i,1}){forcing.PARA.post_proc_class_index(i,1),1});
                    post_proc_class = finalize_init(post_proc_class, tile);
                    forcing = post_process(post_proc_class, forcing, tile);
                end
            end
            
        end
        
        
        function forcing = interpolate_forcing(forcing, tile)
            forcing = interpolate_forcing@FORCING_base(forcing, tile);
                        
%             forcing.TEMP.rainfall = forcing.TEMP.rainfall + double(forcing.TEMP.Tair > 2) .* forcing.TEMP.snowfall;  %reassign unphysical snowfall
%             forcing.TEMP.snowfall = double(forcing.TEMP.Tair <= 2) .* forcing.TEMP.snowfall;
        end

        
        
        function forcing = adjust_precip_quantile_mapping(forcing, carrier_class, tile)
            
            carrier_class.DATA.precip(carrier_class.DATA.precip<0.001)=0;
            forcing.DATA.precip(forcing.DATA.precip<0.001)=0;
            forcing_index = [1:size(forcing.DATA.timeForcing,1)]';
            
            for month=1:12
                
                carrier_res=[];
                result_res=[];
                result_index = [];
                for i=year(forcing.PARA.start_overlap):year(forcing.PARA.end_overlap)
                    range=carrier_class.DATA.timeForcing>=datenum(i,month,1) & carrier_class.DATA.timeForcing<datenum(i,month+1,1);
                    carrier_res=[carrier_res; carrier_class.DATA.precip(range)];
                    
                    range=forcing.DATA.timeForcing>=datenum(i,month,1) & forcing.DATA.timeForcing<datenum(i,month+1,1);
                    result_res=[result_res; forcing.DATA.precip(range)];
                    result_index=[result_index; forcing_index(range)];
                end
                
                carrier_res(carrier_res==0)=[];
                result_index(result_res==0)=[];
                result_res(result_res==0)=[];
                
                carrier_res=carrier_res./mean(carrier_res);
                
                carrier_sorted=sort(carrier_res);
                [result_sorted, index_result_sorted] = sort(result_res);
                
                carrier_sorted=carrier_sorted./mean(carrier_sorted);
                result_sorted=result_sorted ./mean(result_sorted);
                
                
                number_of_quantiles = round(min(size(result_sorted,1), size(carrier_sorted,1))./4);
                
                incremement_carrier_res = size(carrier_res,1)./number_of_quantiles;
                incremement_result_res = size(result_res,1)./number_of_quantiles;
                
                res=[];
                
                for i=1:number_of_quantiles
                    range_carrier = [round((i-1).* incremement_carrier_res+1):round(i.*incremement_carrier_res)]';
                    range_result = [round((i-1).* incremement_result_res+1):round(i.*incremement_result_res)]';
                    res=[res; [i mean(result_sorted(range_result)) mean(carrier_sorted(range_carrier,1)) mean(result_sorted(range_result))./mean(carrier_sorted(range_carrier,1))]];
                end
                
                %apply
                
                start_year=year(forcing.DATA.timeForcing(1,1));
                end_year=year(forcing.DATA.timeForcing(end,1));
                number_of_chunks = max(1, round((end_year-start_year)./(year(forcing.PARA.end_overlap)-year(forcing.PARA.start_overlap))));
                year_slices=round(start_year+[0:number_of_chunks].*((end_year-start_year)./number_of_chunks));
                year_slices(end)=year_slices(end)+1;
                
                for k=1:size(year_slices,2)-1
                    
                    result_res=[];
                    result_index = [];
                    for i=year_slices(k):year_slices(k+1)-1
                        range=forcing.DATA.timeForcing>=datenum(i,month,1) & forcing.DATA.timeForcing<datenum(i,month+1,1);
                        result_res=[result_res; forcing.DATA.precip(range)];
                        result_index=[result_index; forcing_index(range)];
                    end
                    
                    result_index(result_res==0)=[];
                    result_res(result_res==0)=[];
                    
                    [result_sorted, index_result_sorted] = sort(result_res);
                    incremement_result_res = size(result_res,1)./number_of_quantiles;
                    
                    for i=1:number_of_quantiles
                        range_result = [round((i-1).* incremement_result_res+1):round(i.*incremement_result_res)]';
                        forcing.DATA.precip(result_index(index_result_sorted(range_result))) = forcing.DATA.precip(result_index(index_result_sorted(range_result))) ./ res(i,4);
                    end
                end
            end
            
        end
        
        


        %-------------param file generation-----
        function forcing = param_file_info(forcing)
            forcing = provide_PARA(forcing);

            forcing.PARA.STATVAR = [];
            forcing.PARA.class_category = 'FORCING';
            
            forcing.PARA.comment.filename = {'filename of Matlab file containing forcing data'};
            
            forcing.PARA.default_value.forcing_path = {'../CryoGridCommunity_forcing/'};
            forcing.PARA.comment.forcing_path = {'path where forcing data file is located'};
            
            forcing.PARA.comment.start_time = {'start time of the simulations (must be within the range of data in forcing file) - year month day'};
            forcing.PARA.options.start_time.name =  'H_LIST';
            forcing.PARA.options.start_time.entries_x = {'year' 'month' 'day'};
            
            forcing.PARA.comment.end_time = {'end_time time of the simulations (must be within the range of data in forcing file) - year month day'};
            forcing.PARA.options.end_time.name =  'H_LIST'; % 
            forcing.PARA.options.end_time.entries_x = {'year' 'month' 'day'};
            
            forcing.PARA.default_value.rain_fraction = {1};  
            forcing.PARA.comment.rain_fraction = {'rainfall fraction assumed in sumulations (rainfall from the forcing data file is multiplied by this parameter)'};
            
            forcing.PARA.default_value.snow_fraction = {1};  
            forcing.PARA.comment.snow_fraction = {'snowfall fraction assumed in sumulations (rainfall from the forcing data file is multiplied by this parameter)'};

            forcing.PARA.default_value.heatFlux_lb = {0.05};
            forcing.PARA.comment.heatFlux_lb = {'heat flux at the lower boundary [W/m2] - positive values correspond to energy gain'};
            
            forcing.PARA.default_value.airT_height = {2};  
            forcing.PARA.comment.airT_height = {'height above ground surface where air temperature from forcing data is applied'};
            
            forcing.PARA.comment.post_proc_class = {'list of postprocessing classes to modify forcing data in user-defined ways; no post-processing applied when empty'};
            forcing.PARA.options.post_proc_class.name = 'H_LIST';
            
            forcing.PARA.comment.post_proc_class_index = {''};
            forcing.PARA.options.post_proc_class_index.name = 'H_LIST';
        end
    end
end