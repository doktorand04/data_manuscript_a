% Creating a divergent color spectrum based on the data
function spectrum = custom_spectrum_diverging(color_matrix, your_data)
    % Find the number of colors
    num_colors = size(color_matrix,1);

    max_data = max(your_data); % The maximum of your data
    min_data = min(your_data); % The minimum of your data
    
    %Create a proportion based on the maximum and minimum of your data (It
    %handles only negative or only positive data)
    pos_prop = round(256*max(max_data,0)/(max(max_data,0)-min(min_data,0)));
    neg_prop = 256-pos_prop;
   
    
    % Creating the positive and negative spectrum
    pos_color_matrix = color_matrix((num_colors+1)/2:num_colors,:);
    neg_color_matrix = color_matrix(1:(num_colors+1)/2,:);
    
    pos_spectrum = custom_spectrum(pos_color_matrix);
    neg_spectrum = custom_spectrum(neg_color_matrix);
    
    if pos_prop == 0
        spectrum = neg_spectrum;
    elseif neg_prop == 0
        spectrum = pos_spectrum;
    else
        spectrum = [neg_spectrum((257-neg_prop):256,:);pos_spectrum(1:pos_prop,:)]; % This creates a divergent spectrum 
    end
end