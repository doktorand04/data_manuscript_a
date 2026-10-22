clear; close all; clc

trd_off = [-5.95244, 24.22741]; % The extreme generalist-favoring and specialist favoring trade-off
b = 4; % The shape parameter

% Colors for the plots
colors = [0,0,0]/256; 
ax_colors = [213,62,79
            50,136,189]/256;
scatter_colors = [230,97,1
                  94,60,153]/256;

% Plotting affinities against each other
for index = 1:length(trd_off)
    figure(index)
    hold on
    plot_affinity(trd_off(index),b,colors,scatter_colors)
    
    xlim([0 6])
    ylim([0 6])
    xticks(linspace(0,6,7))
    yticks(linspace(0,6,7))
    
    axis("square")
    set(gcf,"OuterPosition",[100 100 1000 850])
    ax = gca;
    ax.XColor = ax_colors(1,:);
    ax.YColor = ax_colors(2,:);
    hold off
end

% Plotting affinities as functions of trait.
for index = 1:length(trd_off)
    figure(index + length(trd_off))
    hold on
    trait_vs_affinity(trd_off(index),b,ax_colors,scatter_colors)
    
    axis("square")
    set(gcf,"OuterPosition",[100 100 1000 850])
    hold off
end

%% Necessary functions
function plot_affinity(trd_off,b,colors,scatter_colors)
    traitpool = linspace(-1,1,101);
    affinitypool = [1+(b.*log(0.5*(exp(-trd_off.*traitpool)+1))-trd_off.*traitpool)/(trd_off-b*log(0.5*(exp(-trd_off)+1)));1+(b.*log(0.5*(exp(trd_off.*traitpool)+1))+trd_off.*traitpool)/(trd_off-b*log(0.5*(exp(-trd_off)+1)))];
    hold on
    plot(affinitypool(1,:),affinitypool(2,:),"Color",colors,"LineWidth",6)
    scatter(affinitypool(1,[15,86]),affinitypool(2,[15,86]),800,scatter_colors,"filled")
    plot([0,1],[1,1],LineStyle="--",LineWidth=4,Color="black",HandleVisibility="off")
    plot([1,1],[0,1],LineStyle="--",LineWidth=4,Color="black",HandleVisibility="off")
    box(gca,'on');
    set(gca,'FontSize',40,'LineWidth',5, 'FontWeight', 'bold');
end

function trait_vs_affinity(trd_off,b,colors,scatter_colors)
    traitpool = linspace(-1,1,101);
    affinitypool = [1+(b.*log(0.5*(exp(-trd_off.*traitpool)+1))-trd_off.*traitpool)/(trd_off-b*log(0.5*(exp(-trd_off)+1)));1+(b.*log(0.5*(exp(trd_off.*traitpool)+1))+trd_off.*traitpool)/(trd_off-b*log(0.5*(exp(-trd_off)+1)))];
    cla
    hold on 
    yyaxis left
    plot(traitpool,affinitypool(1,:),"LineWidth",5,Color=colors(1,:))
    ylim([0 6])
    yticks(linspace(0,6,7))
    ax = gca;
    ax.YColor = colors(1,:);
    yyaxis right
    plot(traitpool,affinitypool(2,:),"LineWidth",5,Color=colors(2,:))
    ylim([0 6])
    yticks(linspace(0,6,7))
    ax = gca;
    ax.YColor = colors(2,:);
    scatter(traitpool([15,86]),affinitypool(1,[15 86]),800,scatter_colors,"filled")
    scatter(traitpool([15,86]),affinitypool(2,[15 86]),800,scatter_colors,"filled")
    plot([0,0],[0,1],LineStyle="--",LineWidth=4,Color="black",HandleVisibility="off",Marker="none")
    plot([-1,0],[1,1],LineStyle="--",LineWidth=4,Color="black",HandleVisibility="off",Marker="none")
    box(gca,'on');
    set(gca,'FontSize',40,'LineWidth',5, 'FontWeight', 'bold');
end

