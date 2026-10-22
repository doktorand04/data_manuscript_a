clear; close all; clc
codeDir   = fileparts(mfilename("fullpath"));   % .../matlab_code
addpath(fullfile(codeDir,"crameri"))
parentDir = fileparts(codeDir);                 % shared parent folder
dataDir   = fullfile(parentDir,"affinity_tradeoff_substitutability");

% Combined results matrix
grid_x = 54;
grid_y = 18;
combined_abundance        = zeros(grid_y,grid_x);
combined_abundance_gen    = zeros(grid_y,grid_x);
combined_abundance_hom    = zeros(grid_y,grid_x);
combined_resource_1       = zeros(grid_y,grid_x);
combined_resource_2       = zeros(grid_y,grid_x);
combined_biodiversity     = zeros(grid_y,grid_x);
combined_rao_entropy      = zeros(grid_y,grid_x);
combined_trait_width      = zeros(grid_y,grid_x);
combined_niche_alignment  = zeros(grid_y,grid_x);

landscape_vector = 1 : 4;
for landscape = landscape_vector

    % Get the results from the files
    results = readtable(fullfile(dataDir,"resource_landscape_"+num2str(landscape),"results_combined"));

    % Get the necessary parameters from the table.
    trade_off        = reshape(results.trade_off,grid_y,grid_x);
    substitutability = reshape(results.substitutability,grid_y,grid_x);
    abundance        = reshape(results.consumer_biomass,grid_y,grid_x);
    abundance_gen    = reshape(results.consumer_biomass_gen,grid_y,grid_x);
    abundance_hom    = reshape(results.consumer_biomass_hom,grid_y,grid_x);
    resource_1       = reshape(results.resource_1_biomass,grid_y,grid_x);
    resource_2       = reshape(results.resource_2_biomass,grid_y,grid_x);
    biodiversity     = reshape(results.evolved_species,grid_y,grid_x);
    rao_entropy      = reshape(results.rao_entropy,grid_y,grid_x);
    trait_width      = reshape(results.trait_width,grid_y,grid_x);
    niche_alignment  = reshape(results.resource_alignment,grid_y,grid_x);
    % Combining the results
    combined_abundance        = combined_abundance + abundance;
    combined_abundance_gen    = combined_abundance_gen + abundance_gen;
    combined_abundance_hom    = combined_abundance_hom + abundance_hom;
    combined_resource_1       = combined_resource_1 + resource_1;
    combined_resource_2       = combined_resource_2 + resource_2;
    combined_biodiversity     = combined_biodiversity + biodiversity;
    combined_rao_entropy      = combined_rao_entropy + rao_entropy;
    combined_trait_width      = combined_trait_width + trait_width;
    combined_niche_alignment  = combined_niche_alignment + niche_alignment;
end

% Combined results matrix
num_of_landscapes         = numel(landscape_vector);
combined_abundance        = combined_abundance/num_of_landscapes;
combined_abundance_gen    = combined_abundance_gen/num_of_landscapes;
combined_abundance_hom    = combined_abundance_hom/num_of_landscapes;
combined_resource_1       = combined_resource_1/num_of_landscapes;
combined_resource_2       = combined_resource_2/num_of_landscapes;
combined_biodiversity     = combined_biodiversity/num_of_landscapes;
combined_rao_entropy      = combined_rao_entropy/num_of_landscapes;
combined_trait_width      = combined_trait_width/num_of_landscapes;
combined_niche_alignment  = combined_niche_alignment/num_of_landscapes;
curvature_trd    = trd_crv(trade_off);
curvature_sub    = sub_crv(substitutability);

%% Plots

%% Figure 1 : Biodiversity in the evolved community
colors        = crameri("lajolla",256);
fig1          = figure(1);
hold on
surf(curvature_trd,curvature_sub,combined_biodiversity,EdgeColor="none",FaceColor="texturemap");
plot(linspace(0,2),linspace(0,2))
colormap(colors)

clim([1,12]);
view(2)
fig_prop(fig1)
hold off

%% Figure 2 : Biodiversity in the evolved community (Rao Q)
colors = crameri("lajolla",256);
fig2 = figure(2);
surf(curvature_trd,curvature_sub,combined_rao_entropy,EdgeColor="none",FaceColor="texturemap");
colormap(colors)
view(2)
fig_prop(fig2)

%% Figure 3 : Biodiversity in the evolved community (Trait Width)
colors = crameri("lajolla",256);
fig3 = figure(3);
surf(curvature_trd,curvature_sub,combined_trait_width,EdgeColor="none",FaceColor="texturemap");
colormap(colors)
view(2)
fig_prop(fig3)

%% Figure 4 : Average consumer biomass in the evolved community

colors = [255,255,229
          247,252,185
          217,240,163
          173,221,142
          120,198,121
          65,171,93
          35,132,67
          0,104,55
          0,69,41]; % Generate a set of distinct colors
spectrum = custom_spectrum(colors);
fig4 = figure(4);
surf(curvature_trd,curvature_sub,combined_abundance,EdgeColor="none",FaceColor="texturemap");
colormap(spectrum);
clim([1.3,1.9])
view(2)
fig_prop(fig4)

%% Figure 5 : Resource abundances in the evolved community
colors = [255,255,229
          247,252,185
          217,240,163
          173,221,142
          120,198,121
          65,171,93
          35,132,67
          0,104,55
          0,69,41]; % Generate a set of distinct colors
spectrum = custom_spectrum(colors);
fig5 = figure(5);
surf(curvature_trd,curvature_sub,combined_resource_1+combined_resource_2,EdgeColor="none",FaceColor="texturemap");
colormap(spectrum);
clim([0.1,0.7])
view(2)
fig_prop(fig5)

%% Figure 6 : Consumer biomass of the evolved community in heterogeneous environment w.r.t homogeneous environment
colors = crameri("vik",255)*255;
spectrum = custom_spectrum_diverging(colors,[-18,0.5]);
fig6 = figure(6);
surf(curvature_trd,curvature_sub,combined_abundance_hom*100,EdgeColor="none",FaceColor="texturemap")
colormap(spectrum);
clim([-18,0.25])
view(2)
fig_prop(fig6);

%% Figure 7 : Consumer biomass of the evolved community in heterogeneous environment w.r.t the best monoculture
colors = crameri("vik",255)*255;
spectrum = custom_spectrum_diverging(colors,[-5,16]);
fig7 = figure(7);
surf(curvature_trd,curvature_sub,combined_abundance_gen*100,EdgeColor="none",FaceColor="texturemap")
colormap(spectrum);
clim([-5,16])
view(2)
fig_prop(fig7)

%% Figure 8 : Resource alignment
colors = [178,24,43
        214,96,77
        244,165,130
        253,219,199
        255,255,255
        224,224,224
        186,186,186
        135,135,135
        77,77,77]; % Generate a set of distinct colors
spectrum = custom_spectrum_diverging(colors,[-0.4,0.4]);
fig8 = figure(8);
surf(curvature_trd,curvature_sub,combined_niche_alignment,EdgeColor="none",FaceColor="texturemap")
colormap(spectrum);
clim([-0.4,0.4])
view(2)
fig_prop(fig8)

%% Necessary functions
function curvature_trd = trd_crv(trade_off)
    b = 4;
    curvature_trd = b * ((1+b)*trade_off - b*log(0.5*(exp(trade_off)+1)))/(b+2)^2;
end
function curvature_sub = sub_crv(substitutability)
    curvature_sub = (1-substitutability);
end

function fig_prop(fig)
    set(gcf,'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.8, 0.8]);
    daspect([1 1 1])
    xticks([ ])
    yticks([ ])
    box(gca,"on")
    set(gca,"LineWidth",5,"InnerPosition",[0.13,0.11,0.775,0.815])
    fontsize(30,"points")
    fontname("Bahnschrift")
    c=colorbar;
    colorbar(gca,'LineWidth',3,Position=[0.917868522638404,0.190783410138249,0.015611544073166,0.653456221198157])
end