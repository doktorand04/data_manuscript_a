% Convert an Nx3 matrix to 256x3 to create a smooth color spectrum

function spectrum = custom_spectrum(color_matrix,num_target_colors)
    % Normalize colors to [0,1] range
    x_colors = color_matrix / 255;
    
    if nargin == 1
        % Number of target colors
        num_target_colors = 256;
    end
    
    % Original indices
    x_indices = linspace(1, num_target_colors, size(x_colors, 1));
    
    % Target indices
    interp_indices = linspace(1, num_target_colors, num_target_colors);
    
    % Interpolate each color channel
    r_interp = interp1(x_indices, x_colors(:,1), interp_indices, 'linear');
    g_interp = interp1(x_indices, x_colors(:,2), interp_indices, 'linear');
    b_interp = interp1(x_indices, x_colors(:,3), interp_indices, 'linear');
    
    % Combine interpolated channels
    spectrum = [r_interp' g_interp' b_interp'];
end


