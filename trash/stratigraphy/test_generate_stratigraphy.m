T =-10;
layerThick = 0.1;

number_of_tiles = 1;
number_of_cells = 10;
layerThick = repmat(0.1, number_of_cells, number_of_tiles);


block_diameter = [1 0.2 0.04];
frost_weathering = -T .* double(T<0) .* 1e-3 .* [];
number_of_block_classes = size(block_diameter,2);


mineral = ones(number_of_cells, number_of_tiles, 1+number_of_block_classes+3); %bedrock, different block sizes, sand, silt, clay
porosity = [0 repmat(0.5, 1, number_of_block_classes) 0.5 0.6 0.6];

volume = (1-porosity(2:end-3)) .* layerThick;
volume_one_block = pi() ./6 .* block_diameter.^3;
number_of_blocks = volume ./ volume_one_block;
surface_area = number_of_blocks .* pi() .* block_diameter.^2;

