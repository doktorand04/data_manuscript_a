alpha = 0.5;
beta = 0.1;
colors = [140,81,10
          255,255,191
          1,102,94]/256;
fillcolors = [251,180,174
              179,205,227
              204,235,197]/256;
plot_generalist_zngi(alpha,beta,colors,fillcolors);

function plot_generalist_zngi(alpha,beta,colors,fillcolors)
    x = linspace(0,0.5,1001);
    y = linspace(0,0.5,1001);
    [X,Y] = meshgrid(x,y);
    Z1 = (X+Y)/2;
    Z2 = min(X,Y);
    Z3 = 2*X.*Y./(X+Y);
    Z = cat(3,Z1,Z2,Z3);
    for index = 1:3
        ZNGI = figure(index);
        cla
        hold on
        contourf(gca,x,y,Z(:,:,index),[0 beta*(1-alpha)/alpha^2],LineWidth=5)
        colormap(fillcolors)
        xticks([0,0.1,0.2,0.3,0.4,0.5])
        yticks([0,0.1,0.2,0.3,0.4,0.5])
        hold off
        axis("square")
        set(gcf,"OuterPosition",[0 0 1000 1000])
    end
end