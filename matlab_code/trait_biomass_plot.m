clear; close all; clc
codeDir   = fileparts(mfilename("fullpath"));   % .../matlab_code
addpath(fullfile(codeDir,"crameri"))
parentDir = fileparts(codeDir);                 % shared parent folder
dataDir   = fullfile(parentDir,"trait_biomass_data","cross_corr_0");

numbins = 27;
max_files = 160;
color = crameri('lajolla', max_files + 1);

a0_values = [-5.95244, -0.16526, 6.91064, 24.22741];
s_values  = [-1.0, 1.0, 0.05882];
bin_centers = linspace(-4.5 + 1/6, 4.5 - 1/6, numbins);
index = 0;

for a0 = a0_values
    for s = s_values
        index = index+1;
        f = figure(index);
        %f.Visible = 'off';
        f.Position = [100 100 1300 650];
        if s == 1 || s == -1
            shortname = sprintf("a0 = %.5f, s = %.1f", a0, s);
        else
            shortname = sprintf("a0 = %.5f, s = %.5f", a0, s);
        end
        filename = shortname + ", alpha = 0.5,beta = 0.1, gamma = 1.0, delta = 0.0001.xlsx";
        
        % Pre-allocate cell arrays for file-level data gathering
        aff_cell = cell(max_files, 1);
        dens_cell = cell(max_files, 1);
        param_sample = [];
        valid_count = 0;
        
        for ind = 1:max_files
            filepath = fullfile(dataDir, sprintf("resource_landscape_%d", ind), ...
                                "Evolution", "Evolutionary Stable State", filename);
            if ~isfile(filepath)
                continue;
            end
            
            % Read data efficiently
            param              = readtable(filepath, 'Sheet', 'Parameters');
            consumer_density   = readmatrix(filepath, 'Sheet', 'Consumer Densities');
            traits             = readmatrix(filepath, 'Sheet', 'Traits');
            
            if isempty(param_sample)
                param_sample = param; % Retain first valid parameter set for fill calculation
            end
            
            consumer_abundance = mean(consumer_density, 1)';
            affinities         = trait_to_affinity(param.a0, param.b, traits);
            distance_metric    = sign(traits) .* sqrt((affinities(:,1) - 1).^2 + (affinities(:,2) - 1).^2);
            trait_bin          = ceil((distance_metric + 4.5) * numbins / 9);
            
            % Aggregate biomass by bin for this file using accumarray
            uni_bins = unique(trait_bin);
            biomass_bin = accumarray(trait_bin, consumer_abundance, [numbins 1]);
            
            valid_count = valid_count + 1;
            aff_cell{valid_count}  = uni_bins;
            dens_cell{valid_count} = biomass_bin(uni_bins);
        end
        
        if valid_count == 0
            continue;
        end
        
        % Combine pre-allocated data
        affinity_bin         = vertcat(aff_cell{1:valid_count});
        consumer_density_bin = vertcat(dens_cell{1:valid_count});
        
        % Bin statistics using accumarray
        incidence    = accumarray(affinity_bin, 1, [numbins 1]);
        mean_biomass = accumarray(affinity_bin, consumer_density_bin, [numbins 1], @mean);
        std_biomass  = accumarray(affinity_bin, consumer_density_bin, [numbins 1], @std);
        
        % Normalize biomass metrics
        norm_factor  = sum(mean_biomass, 'omitnan');
        if norm_factor > 0
            mean_biomass = mean_biomass / norm_factor;
            std_biomass  = std_biomass / norm_factor;
        end
        
        % Feasible region boundaries
        aff_specialists      = trait_to_affinity(param_sample.a0, param_sample.b, [-1; 1]);
        distance_specialists = sign([-1; 1]) .* sqrt((aff_specialists(:,1) - 1).^2 + (aff_specialists(:,2) - 1).^2);
        
        fill_x_1 = [-4.7, -4.7, distance_specialists(1), distance_specialists(1)];
        fill_x_2 = [distance_specialists(2), distance_specialists(2), 4.7, 4.7];
        fill_y   = [0, 1.7, 1.7, 0];
        
        % Plotting
        ax = axes(f);
        hold(ax, 'on');
        
        trbio = bar(ax, bin_centers, mean_biomass, 'FaceColor', 'flat', 'LineWidth', 3);
        trbio.CData = color(incidence + 1, :);
        
        errorbar(ax, bin_centers, mean_biomass, std_biomass, 'k', 'LineStyle', 'none', 'LineWidth', 3);
        fill(ax, fill_x_1, fill_y, 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        fill(ax, fill_x_2, fill_y, 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        
        xticks(ax, [-4.5, -2.25, 0, 2.25, 4.5]);
        xlim(ax, [-4.7, 4.7]);
        ylim(ax, [0, 1.05]);
        box(ax, 'on');
        set(ax, 'FontSize', 50, 'LineWidth', 3, 'FontWeight', 'bold');
        hold(ax, 'off');
        
    end
end

%% Creating the colorbar
total = 160;
color = crameri("lajolla",total+1);
figure;
set(gcf,"Position",[0 0 1400 1400])
% Create invisible axes just to host the colorbar
ax = axes('Visible','off');

% Attach a colorbar to those axes
c = colorbar(ax);

% Set the limits for the colorbar
c.Limits = [0 1];   % example limits

% Optional: set ticks and labels
c.Ticks = 0:.2:1;
c.TickLabels = arrayfun(@num2str, c.Ticks*100, 'UniformOutput', false) + "%";
set(gca,'FontSize',30,'LineWidth',3, 'FontWeight', 'bold');
c.Position([1,1,1,1])

% Optional: choose a colormap
colormap(color);  % or jet, hot, etc.


%% The function to calculate affinity from traits
function affinity = trait_to_affinity(trd_off,b,traits)
    affinity = [1+(b.*log(0.5*(exp(-trd_off.*traits)+1))-trd_off.*traits)/(trd_off-b*log(0.5*(exp(-trd_off)+1))),1+(b.*log(0.5*(exp(trd_off.*traits)+1))+trd_off.*traits)/(trd_off-b*log(0.5*(exp(-trd_off)+1)))];
end