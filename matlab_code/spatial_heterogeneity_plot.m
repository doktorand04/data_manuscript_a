% Import all necessary metrics from excel

richness = [3.997175141	4.836614173	2.030612245	6.057142857;
            3.033898305	4.155511811	2.025510204	4.857142857;
            1.098870056	1.408464567	2	1.778571429];

biomass  = [1.592072532	1.370310968	1.757486469	1.386812763;
            1.596200905	1.462332647	1.757796054	1.477292809;
            1.603212205	1.595710624	1.756990368	1.616502476];

resource = [0.407927469	0.629691402	0.24251353	0.61320372;
            0.403799096	0.53766845	0.242203954	0.522715582;
            0.396787795	0.404289378	0.243009634	0.383497582];

rel_biomass =  [0.002683217	-0.00655255	0.100209509	-0.015284756;
                0.001376661	-0.005252753	0.096128259	-0.013200936;
                4.7816E-06	-0.00037074	0.089775089	-0.001222925];

hom_biomass =  [-0.004954667	-0.143652037	0.001409115	-0.145074012;
                -0.002374434	-0.086145316	0.001586607	-0.089296265;
                0.002007628	-0.002794488	0.001115788	-0.003479111];

func_rich = [0.790031311	1.066154699	5.879680319	3.512353069;
             0.617695785	0.89779184	5.853867002	3.196456011;
             0.040233755	0.196156677	5.8135471	1.164948671];

trait_align = [0.283107979	-0.304569468	0.229638492	-0.406466667;
               0.129482066	-0.204334195	0.120092561	-0.293003861;
               0.000917229	-0.013499599	0.005569852	-0.020800778];


% Do barchart for each variable
for ind = 1:4
    if ind == 4 || ind == 2
        figure("Color",[253,219,199]/256)
    else
        figure("Color",[209,229,240]/256)
    end
    b = bar(richness(:,ind),LineWidth = 3.5);
    b.FaceColor = 'flat';
    b.CData(1,:) = [128,0,38]/256;   % red
    b.CData(2,:) = [253,141,60]/256;   % green
    b.CData(3,:) = [255,255,204]/256;   % blue
    set(gca,'FontSize',80,'LineWidth', 4, 'FontWeight', 'bold');
    ylim([0 7])
    xticks([])
end

%%
for ind = 1:4
    if ind == 4 || ind == 2
        figure("Color",[253,219,199]/256)
    else
        figure("Color",[209,229,240]/256)
    end
    b = bar(biomass(:,ind),LineWidth = 3.5);
    b.FaceColor = 'flat';
    b.CData(1,:) = [128,0,38]/256;   % red
    b.CData(2,:) = [253,141,60]/256;   % green
    b.CData(3,:) = [255,255,204]/256;   % blue
    set(gca,'FontSize',80,'LineWidth', 4, 'FontWeight', 'bold');
    ylim([1.2 2])
    yticks([1.2 1.6 2.0])
    xticks([])
end

%%
for ind = 1:4
    if ind == 4 || ind == 2
        figure("Color",[253,219,199]/256)
    else
        figure("Color",[209,229,240]/256)
    end
    b = bar(resource(:,ind),LineWidth = 3.5);
    b.FaceColor = 'flat';
    b.CData(1,:) = [128,0,38]/256;   % red
    b.CData(2,:) = [253,141,60]/256;   % green
    b.CData(3,:) = [255,255,204]/256;   % blue
    set(gca,'FontSize',80,'LineWidth', 4, 'FontWeight', 'bold');
    ylim([0 0.8])
    yticks([0 0.4 0.8])
    xticks([])
end

%%
for ind = 1:4
    if ind == 4 || ind == 2
        figure("Color",[253,219,199]/256)
    else
        figure("Color",[209,229,240]/256)
    end
    b = bar(rel_biomass(:,ind)*100,LineWidth = 3.5);
    b.FaceColor = 'flat';
    b.CData(1,:) = [128,0,38]/256;   % red
    b.CData(2,:) = [253,141,60]/256;   % green
    b.CData(3,:) = [255,255,204]/256;   % blue
    set(gca,'FontSize',80,'LineWidth', 4, 'FontWeight', 'bold');
    ylim([-2 12])
    ytickformat('percentage');
    xticks([])
end

%%
for ind = 1:4
    if ind == 4 || ind == 2
        figure("Color",[253,219,199]/256)
    else
        figure("Color",[209,229,240]/256)
    end
    b = bar(hom_biomass(:,ind)*100,LineWidth = 3.5);
    b.FaceColor = 'flat';
    b.CData(1,:) = [128,0,38]/256;   % red
    b.CData(2,:) = [253,141,60]/256;   % green
    b.CData(3,:) = [255,255,204]/256;   % blue
    set(gca,'FontSize',80,'LineWidth', 4, 'FontWeight', 'bold');
    ylim([-16 1])
    ytickformat('percentage');
    xticks([])
end

%%
for ind = 1:4
    if ind == 4 || ind == 2
        figure("Color",[253,219,199]/256)
    else
        figure("Color",[209,229,240]/256)
    end
    b = bar(func_rich(:,ind),LineWidth = 3.5);
    b.FaceColor = 'flat';
    b.CData(1,:) = [128,0,38]/256;   % red
    b.CData(2,:) = [253,141,60]/256;   % green
    b.CData(3,:) = [255,255,204]/256;   % blue
    set(gca,'FontSize',80,'LineWidth', 4, 'FontWeight', 'bold');
    ylim([0 7])
    xticks([])
end

%%
for ind = 1:4
    if ind == 4 || ind == 2
        figure("Color",[253,219,199]/256)
    else
        figure("Color",[209,229,240]/256)
    end
    b = bar(trait_align(:,ind),LineWidth = 3.5);
    b.FaceColor = 'flat';
    b.CData(1,:) = [128,0,38]/256;   % red
    b.CData(2,:) = [253,141,60]/256;   % green
    b.CData(3,:) = [255,255,204]/256;   % blue
    set(gca,'FontSize',80,'LineWidth', 4, 'FontWeight', 'bold');
    ylim([-0.45 0.3])
    xticks([])
end