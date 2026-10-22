% Plotting the consumer abundance, resource abundance and resource supplies in space
clear; close all; clc
codeDir   = fileparts(mfilename("fullpath"));   % .../matlab_code
addpath(fullfile(codeDir,"crameri"))
parentDir = fileparts(codeDir);                 % shared parent folder
dataDir   = fullfile(parentDir,"affinity_tradeoff_substitutability");
a0 = -0.16526; s = [-1.0, 1.0];
for index = 1:2
    filename = fullfile(dataDir,"resource_landscape_1", "Evolution", "Evolutionary Stable State", "a0 = "+num2str(a0)+", s = "+num2str(s(index),'%.1f')+", alpha = 0.5,beta = 0.1, gamma = 1.0, delta = 0.0001.xlsx");
    spatial_plot(filename);
end

function spatial_plot(filename)
    % Check if savefolder is provided, otherwise set save to false
    
    % Read data from the specified Excel file
    simulation  = readtable(filename,"Sheet","Consumer Densities"); % Getting consumer densities
    inflow      = readtable(filename,"Sheet","Inflow Rates"); % Getting resource densities
    traits      = readtable(filename,"Sheet","Traits"); % Getting the traits of the evolved communities
    resource    = readtable(filename,"Sheet","Resource Densities"); % Getting resource densities
    
    % Small epsilon value to avoid division by zero or log(0)
    eps = 1e-10;
    
    % Calculate various ecological metrics
    abundance        = table2array(simulation);% Getting consumer densities in matrix form
    abundance_limits = [min(min(abundance)) max(max(abundance))];  % limits of abundances for the colorbar
    inflow_abundance = table2array(inflow); % Getting resource supply densities in matrix form
    res_abundance    = table2array(resource);% Getting resource densities in matrix form
    
    % Sort traits and rearrange abundance accordingly
    [traits.traits,ind] = sort(traits.traits); % sorting the traits;
    abundance           = abundance(:,ind);% sorting the abundance;
    
    % Define grid size and interpolation points for spatial plotting
    grid_size_updated   = 128;%sqrt(length(tot_abundance)); % change to the grid_size you need
    interp_points       = 1/grid_size_updated:1/grid_size_updated:1; % the interpolated grid (one dimension)
    
    % Threshold for plotting species abundance
    threshold = 1e-4;
    
    % Define a custom colormap for the the plots
    customMap = [255,255,255
                237,248,177
                199,233,180
                127,205,187
                65,182,196
                29,145,192
                34,94,168
                37,52,148
                8,29,88];
    customMap = custom_spectrum(customMap);
    customMap = [1,1,1
               customMap];  
    
    % Plot the resource supplies
    for index = 1:2
        fig = figure;
        clf
        set(gcf,"Position",[0 0 1400 1400])
        colormap(customMap);
        
        % Reshape and interpolate resource supply densities
        inflow_interp    = reshape(inflow_abundance(:,index),sqrt(length(inflow_abundance(:,index))),[]);
        inflow_interp    = periodic_interpolation_2D(inflow_interp,interp_points,1);
        
        % Create the image plot
        imagesc(inflow_interp)
        axis("square")
        xticks([])
        yticks([])
        colorbar
        fontsize(40,"pixels")
        clim([0 1.7])
    end
                
    % Plot the residual resource densities
    for index = 1:2
        fig = figure;
        clf
        set(gcf,"Position",[0 0 1400 1400])
        colormap(customMap);
        
        % Reshape and interpolate resource supply densities
        resource_interp    = reshape(res_abundance(:,index),sqrt(length(res_abundance(:,index))),[]);
        resource_interp    = periodic_interpolation_2D(resource_interp,interp_points,1);
        
        % Create the image plot
        imagesc(resource_interp)
        axis("square")
        xticks([])
        yticks([])
        colorbar
        fontsize(40,"pixels")
        clim([0 1.7])
    end

    % Plot the consumer biomass of individual species
    for index = 1:length(traits.traits)
        fig = figure;
        set(gcf,"Position",[0 0 1400 1400])
        colormap(customMap);
        % Reshape, interpolate, and threshold species abundance
        species_abundance = reshape(max(abundance(:,index),0),sqrt(length(abundance)),[]);
        species_abundance = periodic_interpolation_2D(species_abundance,interp_points,1);
        species_abundance(species_abundance < threshold) = 0;
        
        % Create the image plot
        imagesc(species_abundance,abundance_limits)
        clim([0 1.7])
        axis("square")
        title("\chi = "+num2str(traits.traits(index)))
        fontsize(40,"pixels")
        xlim([1 grid_size_updated])
        ylim([1 grid_size_updated])
        xticks([])
        yticks([])
    end
    
end

% Function to perform 2D periodic interpolation
function FF = periodic_interpolation_2D( F , xx , L )
    if nargin < 3
        L = 2*pi;
    end
    
    xx = xx(:);
    
    N = size( F , 1 );
    h = L/N;
    x = h*(1:N)';
    
    Lsf = (2*pi)/L;
    
    Xd = xx*ones( 1 , length(x) ) - ones( length(xx) , 1 ) * x';
    
    S = sin( Lsf*Xd*N*0.5 ) ./ ( N * ( tan( Lsf*Xd*0.5 ) ) );
    S( Xd == 0 ) = 1;
    
    FF = S*F*(S');
end