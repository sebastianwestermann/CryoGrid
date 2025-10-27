
t_store=[];
Sin_store=[];
Sun_el_store = [];
tile.PARA.latitude = 47.989722000000000;


hour_angle = rad2deg(DATA.hour_angle);
hour_angle(hour_angle<0) = hour_angle(hour_angle<0)+360;
hour_angle(hour_angle>=360) = hour_angle(hour_angle>=360)-360;

for t=[0:1/100:365]+forcing.DATA.timeForcing(1,1) %[0:1/100:365]
    t_store=[t_store; t];
    posit=floor((t-forcing.DATA.timeForcing(1,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1)))+1;

    t_weight = (t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1)); % current distance from last timestep (0,1)

        
    if hour_angle(posit+1,1) <  hour_angle(posit,1)
        ha = hour_angle(posit,:) + (hour_angle(posit+1,:) + 360 - hour_angle(posit,:)).*t_weight;
    else
        ha = hour_angle(posit,:) + (hour_angle(posit+1,:) - hour_angle(posit,:)).*t_weight;
    end

    DecR = forcing.DATA.DecR(posit,:)+(forcing.DATA.DecR(posit+1,:)-forcing.DATA.DecR(posit,:)).*t_weight;
    sunElevation = asind(cosd(tile.PARA.latitude) .* cos(DecR) .* cosd(ha) + sind(tile.PARA.latitude) .* sin(DecR));
 
    S_TOA = 1370.*max(sind(sunElevation),0);

    S_TOA1 = 1370.*max(sind(forcing.DATA.sunElevation(posit,:)),0);
    S_TOA2 = 1370.*max(sind(forcing.DATA.sunElevation(posit+1,:)),0);

    attenuation1 = forcing.DATA.Sin(posit,:)./S_TOA1;
    attenuation2 = forcing.DATA.Sin(posit+1,:)./S_TOA2;
    if isnan(attenuation1) && isnan(attenuation2)
        Sin = 0;
    else
        if isnan(attenuation1)
            Sin = attenuation2 .* S_TOA;
        elseif isnan(attenuation2)
            Sin = attenuation1 .* S_TOA;
        else
            attenuation = attenuation1 + (attenuation2 - attenuation1).*t_weight;
            Sin = attenuation .* S_TOA;
        end
    end

    Sun_el_store = [Sun_el_store; sunElevation];
    Sin_store=[Sin_store; Sin];

end