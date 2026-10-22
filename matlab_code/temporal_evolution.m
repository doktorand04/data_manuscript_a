clear; close all; clc
% Load the results into MATLAB
codeDir   = fileparts(mfilename("fullpath"));   % .../matlab_code
parentDir = fileparts(codeDir);                 % shared parent folder
dataDir   = fullfile(parentDir, "affinity_tradeoff_substitutability", "temporal_evolution");
filename   = fullfile(dataDir,"Evolution", "Evolutionary Stable State", "a0 = -0.16526, s = -1.0, alpha = 0.5,beta = 0.1, gamma = 1.0, delta = 0.0001.xlsx");
addpath(fullfile(codeDir,"crameri"))
colors        = crameri("bamako");

result     = readtable(filename,"Sheet","Traits");
consumers  = length(result.traits);

init_time = 0;
clf
f = figure(1);
hold on
max_all_density = 0;
for cons = 1:consumers
    new_filename = fullfile(dataDir, "a0 = -0.16526, s = -1.0, alpha = 0.5,beta = 0.1, gamma = 1.0, delta = 0.0001, consumer = "+num2str(cons)+".xlsx");
    timetable            = readtable(new_filename,"Sheet","Time");
    traittable           = readtable(new_filename,"Sheet","Traits");
    consumerdensitytable = readtable(new_filename,"Sheet","Consumer Densities");      
    time      = timetable.time+init_time;
    init_time = time(end);
    trait           = table2array(traittable);
    consumerdensity = table2array(consumerdensitytable);
    max_density = max(max(consumerdensity));
    max_all_density = max(max_density,max_all_density);
    for ind = 1: cons
        scatter(time(1),trait(1,ind),25,'black',"filled")
        plot(time,trait,'Color','black',LineWidth=5,LineJoin='chamfer')
        scatter(time(end),trait(end,ind),25,'black',"filled")
    end
end

spectrum = [255,255,229
            247,252,185
            217,240,163
            173,221,142
            120,198,121
            65,171,93
            35,132,67
            0,104,55
            0,69,41];
f.Position = [100 100 1200 650];
ylim([-1 1])
xlim([0,ceil(max(time)/10^9)*10^9])
box(gca,'on');
set(gca,'FontSize',40,'LineWidth',5, 'FontWeight', 'bold');
ax = gca;
ax.XAxis.Exponent = 0;
xtickformat('%.0f')
hold off
