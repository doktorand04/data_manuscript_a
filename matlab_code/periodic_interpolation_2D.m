%% Creating a periodic interpolation of the  any size square grid to any size square grid
function interpolated_values = periodic_interpolation_2D(matrix_valued_function, grid_points, land_length)

    % If length is not specified, the default is 2*pi (period of
    % trigonometric function)
    if nargin < 3
        land_length = 2*pi;
    end
      
    grid_points = grid_points(:); % Convert the grid points into a vector
    grid_size = size(matrix_valued_function, 1); % Size of the grid
    grid_length = land_length/grid_size; % Length of each grid
    int_length = grid_length*(1:grid_size)'; % Converting the length to integer values
    
    period = (2*pi)/land_length; % Period of the function
    
    int_points = grid_points*ones( 1 , length(int_length) ) - ones( length(grid_points) , 1 ) * int_length'; % Interpolated grid points
    
    interpolant = sin( period*int_points*grid_size*0.5 ) ./ ( grid_size * ( tan( period*int_points*0.5 ) ) ); % Creating the interpolant matrix for interpolation (matrix used as a function to create the interpolation)
    interpolant( int_points == 0 ) = 1; % Need to make sure already existing values does not change.
    
    interpolated_values = interpolant*matrix_valued_function*(interpolant'); % The interpolation
end