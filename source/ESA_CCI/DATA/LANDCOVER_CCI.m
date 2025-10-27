classdef LANDCOVER_CCI < matlab.mixin.Copyable

    properties
        PARENT
        CONST
        STATVAR
        PARA
    end

    
    methods
        
        function lc = provide_PARA(lc)
            lc.PARA.landcover_file = [];
            lc.PARA.landcover_path = [];
            
            lc.PARA.yedoma_file = [];
            lc.PARA.yedoma_path = [];
            
            lc.PARA.accumulated{1,1} = [140 150 152 153]; %sparse vegetation
            lc.PARA.accumulated{2,1} = [10  11  12  20 130]; %grasslands and croplands
            lc.PARA.accumulated{3,1} = [30 40 100 110 120 121 122]; %shrubs
            lc.PARA.accumulated{4,1} = [50 60 61 62 80 81 82 90]; %deciduous forest
            lc.PARA.accumulated{5,1} = [70 71 72]; %evergreen forest
            lc.PARA.accumulated{6,1} = [160 170 180]; %wetlands
            lc.PARA.accumulated{7,1} = [190 200 201 202 220]; %bare areas and urban
            
            %Yedoma
            lc.PARA.accumulated{8,1} = [140 150 152 153]; %sparse vegetation
            lc.PARA.accumulated{9,1} = [50 60 61 62 80 81 82 90]; %deciduous forest
            lc.PARA.accumulated{10,1} = [160 170 180]; %wetlands
        end
        
        function lc = provide_STATVAR(lc)

        end
        
        function lc = provide_CONST(lc)
            
        end
        
        function lc = finalize_init(lc)
            
        end
        
        function lc = load_data(lc)
            
            tic
            lc.PARENT.STATVAR.landcover = repmat(lc.PARENT.STATVAR.latitude .* 0, 1,10);
            
            max_lat = max(lc.PARENT.STATVAR.latitude);
            min_lat = min(lc.PARENT.STATVAR.latitude);
            max_lon = max(lc.PARENT.STATVAR.longitude);
            min_lon = min(lc.PARENT.STATVAR.longitude);
            delta_1km_lat  = lc.PARENT.PARA.delta_lat;
            delta_1km_lon  = lc.PARENT.PARA.delta_lon;
            
            delta_300m_lat = 1/360;
            delta_300m_lon = 1/360;
           
           % roi_start_index_lat = max(1, round((90 - (max_lat + delta_1km_lat/2))./delta_300m_lat + 1)-1);
           % roi_start_index_lon = max(1, round((180 + (min_lon - delta_1km_lon/2))./delta_300m_lon + 1)-1);
            roi_start_index_lat = max(1, round((90 - (max_lat + delta_1km_lat/2))./delta_300m_lat + 1)-7);
            roi_start_index_lon = max(1, round((180 + (min_lon - delta_1km_lon/2))./delta_300m_lon + 1)-7);
            
          %  number_of_elements_lat = min(180/delta_300m_lat, round((max_lat - min_lat + delta_1km_lat) ./ delta_300m_lat)+2);
          %  number_of_elements_lon = min(360/delta_300m_lon, round((max_lon - min_lon + delta_1km_lon) ./ delta_300m_lon)+2);
            number_of_elements_lat = min(180/delta_300m_lat-roi_start_index_lat+1, round((max_lat - min_lat + delta_1km_lat) ./ delta_300m_lat)+14);
            number_of_elements_lon = min(360/delta_300m_lon-roi_start_index_lon+1, round((max_lon - min_lon + delta_1km_lon) ./ delta_300m_lon)+14);
            
            landcover = double(ncread([lc.PARA.landcover_path lc.PARA.landcover_file], 'lc', [roi_start_index_lat roi_start_index_lon], [number_of_elements_lat number_of_elements_lon], [1 1]));
            
            for i=1:size(lc.PARENT.STATVAR.latitude,1)
                target_latitude = lc.PARENT.STATVAR.latitude(i,1);
                target_longitude = lc.PARENT.STATVAR.longitude(i,1); %make sure longitude is indeed -180->180
                
                index_1km_lat = max(1, round((90 - target_latitude + delta_1km_lat/2) ./ delta_1km_lat));
                index_1km_lon = max(1, round((180 + target_longitude + delta_1km_lon/2) ./ delta_1km_lon));
                
                start_index_300m_lat = round((index_1km_lat - 1) .*  delta_1km_lat ./ delta_300m_lat.*1e6) ./1e6; %0 ; 3.6
                end_index_300m_lat = round(index_1km_lat .*  delta_1km_lat ./ delta_300m_lat .* 1e6) ./1e6; % 3.6  ; 7.2
	%	end_index_300m_lat(floor(end_index_300m_lat)-floor(start_index_300m_lat)+1 > number_of_elements_lat) = floor(start_index_300m_lat) + number_of_elements_lat - 2;
                start_index_300m_lon = round((index_1km_lon - 1) .*  delta_1km_lon ./ delta_300m_lon .*1e6) ./1e6;
                end_index_300m_lon = round(index_1km_lon .*  delta_1km_lon ./ delta_300m_lon .* 1e6) ./1e6;
	%	end_index_300m_lon(floor(end_index_300m_lon)-floor(start_index_300m_lon)+1 > number_of_elements_lon) = floor(start_index_300m_lon) + number_of_elements_lon - 2;

                
                %start_index_300m_lat2 = ceil(start_index_300m_lat); % floor(start_index_300m_lat + 1); 
                %end_index_300m_lat2 = ceil(end_index_300m_lat); %floor(end_index_300m_lat + 1); 
                %start_index_300m_lon2 = ceil(start_index_300m_lon); %floor(start_index_300m_lon + 1);
                %end_index_300m_lon2 = ceil(end_index_300m_lon); %floor(end_index_300m_lon + 1);
                
                start_index_300m_lat2 = floor(start_index_300m_lat + 1); %ceil(start_index_300m_lat); % ; 
                end_index_300m_lat2 = ceil(end_index_300m_lat); %floor(end_index_300m_lat + 1); 
                start_index_300m_lon2 = floor(start_index_300m_lon + 1); %ceil(start_index_300m_lon);
                end_index_300m_lon2 = ceil(end_index_300m_lon); %floor(end_index_300m_lon + 1);
                
                weight_300m = ones(end_index_300m_lat2 - start_index_300m_lat2 + 1, end_index_300m_lon2 - start_index_300m_lon2 + 1);
                weight_300m(1,:) = weight_300m(1,:) .*(start_index_300m_lat2 - start_index_300m_lat);
                weight_300m(end,:) = weight_300m(end,:) .*(-(end_index_300m_lat2 - 1) + end_index_300m_lat);
                weight_300m(:,1) = weight_300m(:,1) .*(start_index_300m_lon2 - start_index_300m_lon);
                weight_300m(:,end) = weight_300m(:,end) .*(-(end_index_300m_lon2 - 1) + end_index_300m_lon);
                
                weight_300m = weight_300m ./ (delta_1km_lat ./ delta_300m_lat) ./ (delta_1km_lon ./ delta_300m_lon); %normalize
                
                start_index_300m_lat3 = start_index_300m_lat2 - roi_start_index_lat + 1;
                end_index_300m_lat3 = end_index_300m_lat2 - roi_start_index_lat + 1;
                start_index_300m_lon3 = start_index_300m_lon2 - roi_start_index_lon + 1;
                end_index_300m_lon3 = end_index_300m_lon2 - roi_start_index_lon + 1;
                
                if start_index_300m_lat3 <1 || start_index_300m_lon3 <1
                    a=1;
                end
               
                for k = 1:7%size(lc.PARA.accumulated,1) %no yedoma yet
                    land_cover_fraction = weight_300m.* 0;
                    for j=1:size(lc.PARA.accumulated{k,1},2)
                        land_cover_fraction = land_cover_fraction + double(landcover(start_index_300m_lat3:end_index_300m_lat3,start_index_300m_lon3:end_index_300m_lon3) == lc.PARA.accumulated{k,1}(1,j));
                        
                    end
                    lc.PARENT.STATVAR.landcover(i,k) = sum(weight_300m(:) .* land_cover_fraction(:));
                end
                
            end
            %normalize to account for water bodies
            lc.PARENT.STATVAR.landcover = lc.PARENT.STATVAR.landcover ./ repmat(sum(lc.PARENT.STATVAR.landcover,2), 1, size(lc.PARENT.STATVAR.landcover,2));
            toc
            
            
            %add grassland class to tundra class for areas out of Central Asia
            %grasslandCentralAsia = lc.PARENT.STATVAR.landcover(:,2) .* double(lc.PARENT.STATVAR.longitude > 50 & lc.PARENT.STATVAR.latitude < 55);
            [gl_lat, gl_lon] = load_grassland_polygon(lc);
            
            % grasslandCentralAsia = lc.PARENT.STATVAR.landcover(:,2) .* double(lc.PARENT.STATVAR.longitude > -123 & lc.PARENT.STATVAR.latitude < 60); %chnaged after conversation with Brendan O'Neill, so that grasslands in Canada prairies are not considered tundra; 60 degrees seems a good compromise 
            grasslandCentralAsia = lc.PARENT.STATVAR.landcover(:,2) .* double(inpolygon(lc.PARENT.STATVAR.latitude, lc.PARENT.STATVAR.longitude, gl_lat, gl_lon)); %changed after conversation with Brendan O'Neill, so that grasslands in Canada prairies are not considered tundra; 60 degrees seems a good compromise

            grasslandOther = lc.PARENT.STATVAR.landcover(:,2) - grasslandCentralAsia;
            lc.PARENT.STATVAR.landcover(:,1) = lc.PARENT.STATVAR.landcover(:,1) + grasslandOther;
            lc.PARENT.STATVAR.landcover(:,2) = grasslandCentralAsia;
            
            
            %Yedoma map         
            delta_yedoma = double(ncread([lc.PARA.yedoma_path lc.PARA.yedoma_file], 'delta_lat_lon'));
            
            roi_start_index_lat = max(1, round((90 - (max_lat + delta_1km_lat/2))./delta_yedoma + 1)-1);
            roi_start_index_lon = max(1, round((180 + (min_lon - delta_1km_lon/2))./delta_yedoma + 1)-1);
            
            number_of_elements_lat = min(180/delta_yedoma-roi_start_index_lat+1, round((max_lat - min_lat + delta_1km_lat) ./ delta_yedoma)+2);
            number_of_elements_lon = min(360/delta_yedoma-roi_start_index_lon+1, round((max_lon - min_lon + delta_1km_lon) ./ delta_yedoma)+2);
            yedoma = double(ncread([lc.PARA.yedoma_path lc.PARA.yedoma_file], 'yedoma', [roi_start_index_lat roi_start_index_lon], [number_of_elements_lat number_of_elements_lon], [1 1]));
            
            lat_yedoma_start = 90-delta_yedoma/2 - (roi_start_index_lat-1) .* delta_yedoma; %latitude of 1st Yedoma pixel in ROI
            lon_yedoma_start = -180 + delta_yedoma/2 + (roi_start_index_lon-1) .* delta_yedoma; %longitude of 1st Yedoma pixel in ROI
            
            factor_lat = delta_1km_lat ./ delta_yedoma;
            factor_lon = delta_1km_lon ./ delta_yedoma;
            yedoma_lc = lc.PARENT.STATVAR.latitude.*0;
            delta_index_lat = (factor_lat-1)/2; 
            delta_index_lon = (factor_lon-1)/2;
            for i=1:size(lc.PARENT.STATVAR.latitude,1)
                center_pixel_lat = (lat_yedoma_start - lc.PARENT.STATVAR.latitude(i,1)) ./ delta_yedoma + 1;
                center_pixel_lon = (-lon_yedoma_start + lc.PARENT.STATVAR.longitude(i,1)) ./ delta_yedoma +1;
                offset_lat = delta_index_lat;
                offset_lon = delta_index_lon;
                yedoma_lc(i,1) = yedoma_lc(i,1) + sum(sum(yedoma(round(center_pixel_lat-offset_lat):round(center_pixel_lat+offset_lat), round(center_pixel_lon-offset_lon):round(center_pixel_lon+offset_lon))));
            end
            yedoma_lc = yedoma_lc ./ factor_lon ./ factor_lat; %check that normalized
            
            lc.PARENT.STATVAR.landcover(:,8) = lc.PARENT.STATVAR.landcover(:,1) .* yedoma_lc;
            lc.PARENT.STATVAR.landcover(:,1) = lc.PARENT.STATVAR.landcover(:,1) .* (1-yedoma_lc);
            lc.PARENT.STATVAR.landcover(:,9) = lc.PARENT.STATVAR.landcover(:,4) .* yedoma_lc;
            lc.PARENT.STATVAR.landcover(:,4) = lc.PARENT.STATVAR.landcover(:,4) .* (1-yedoma_lc);
            % wetland Yedoma not taken into account now
            %INFO.landcover(:,6) = INFO.landcover(:,6) .* (1-Yedoma_yes_no);
%             lc_list(:,10) = lc_list(:,10) .* 0; %no effect on wetland
            %INFO.landcover(:,6) = INFO.landcover(:,6) .* (1-Yedoma_yes_no);
            %INFO.landcover(:,10) = INFO.landcover(:,10) .* Yedoma_yes_no;

            lc.PARENT.STATVAR.landcover(lc.PARENT.STATVAR.landcover<0) = 0; %eliminate rounding errors
        end

        function [gl_lat, gl_lon] = load_grassland_polygon(lc)
            %grasslands2.kml in CryoGridCommunity_data
            gl_lat= [43.0361880356192
                48.4244008476347
                48.1977146181899
                48.8385034848756
                49.4736662322806
                49.4274439523764
                49.1714473072957
                49.4362297965566
                49.9943108216031
                50.1888471633907
                50.7261690358846
                51.1187349468152
                51.6972226260531
                52.0828005936034
                53.0283451510996
                54.0369388607992
                55.3680302903190
                55.8610745564074
                56.9110409616511
                58.8617357161552
                58.9408828544971
                58.2721540016672
                55.6739173054069
                54.8905989697458
                54.4019921521481
                54.2801730858122
                52.8195607033879
                50.2042237950491
                49.4474109249538
                48.8192577298068
                49.6299192524906
                49.1428550704734
                49.3783461610570
                47.1326803905548
                46.1092732511687
                45.7549925582800
                57.6114663305662
                61.1821669663195
                61.3894827951818
                60.0904783191077
                60.4852177851681
                60.2879315710086
                57.9971230360888
                59.0928914865590
                57.7691867990086
                57.6927445485404
                60.8152113100651
                59.1550929584152
                55.1023254088357
                54.4323695563531
                52.3160183361465
                46.9068649886461
                35.3642800947836
                30.0385371766006
                18.9976657363116
                7.85688758043571
                6.21853399514197
                11.3766129501533
                15.4285866751527
                11.3549782206481
                3.09607307859424
                18.6542293318535
                35.5685312722723];

            gl_lon = [-127.267545642716
                -124.849291394198
                -123.104034173139
                -123.141836390386
                -123.053442038687
                -121.593884934091
                -117.936708875195
                -114.431953444833
                -114.277066636303
                -114.447738724071
                -114.653061811080
                -115.024589124859
                -115.300494867217
                -115.555544157546
                -116.908695385551
                -118.668363652686
                -121.522074426332
                -122.360258989057
                -122.430194709620
                -119.728753436826
                -116.487482426734
                -114.078282362741
                -112.250800532316
                -107.226565575261
                -101.733806963581
                -97.8583889777430
                -91.4739727522816
                -86.9966502749212
                -83.4318557608985
                -77.3333288692785
                -73.1718806092381
                -67.1916814912072
                -64.3831436168145
                -59.3745144464321
                -57.4295285417875
                -52.1066328069659
                7.93638033716456
                17.4308683949312
                38.2485269411215
                49.4632670618069
                58.8808020766097
                72.6584637315669
                79.8278811644131
                84.7236956165751
                95.8059974846535
                102.365688564917
                115.754674987327
                128.603980357315
                136.028186461371
                141.367831508808
                141.442821167139
                140.499314375327
                131.159835298993
                126.268678879034
                119.609247565849
                107.097555440164
                90.4328020294331
                57.0712945196939
                12.6828728621334
                -57.9204543712336
                -81.2975208570115
                -111.048542570139
                -123.746879561166];
        end


    end
end

