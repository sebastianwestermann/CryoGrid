function  temperature = piece_wise_linear(fit_params, year_vector, fixed_params)


year_points = fixed_params;
unconstrained_years = find(isnan(year_points));
year_points(unconstrained_years) = fit_params(1:size(unconstrained_years,1));
temperature_points = fit_params(size(unconstrained_years,1)+1:end, 1);
temperature=[];

segment=1;
for i=1:size(year_vector,1)
    if year_vector(i,1)>= year_points(segment+1,1)
        segment = segment+1;
    end
    temperature = [temperature; temperature_points(segment,1) + ...
        (temperature_points(segment+1,1)-temperature_points(segment,1)) ./ (year_points(segment+1)-year_points(segment,1)) .* (year_vector(i,1) - year_points(segment,1))];
end


%fixed_params = vector length number_of_segments+1
