clear; close all; clc
% When the consumers reach an evolutionar singular point, we take a look at
% the fitness landscape. If we have a possible mutant with positive
% fitness, we let the mutant with the highest fitness invade the system.

% Load the results into MATLAB
codeDir   = fileparts(mfilename("fullpath"));   % .../matlab_code
parentDir = fileparts(codeDir);                 % shared parent folder
dataDir   = fullfile(parentDir, "affinity_tradeoff_substitutability", "temporal_evolution");
filename   = fullfile(dataDir, "Evolution", "a0 = -0.16526, s = -1.0, alpha = 0.5,beta = 0.1, gamma = 1.0, delta = 0.0001, consumer = 4.xlsx");
parametertable        = readtable(filename,Sheet="Parameters");
traittable            = readtable(filename,"Sheet","Traits");
consumerdensitytable  = readtable(filename,"Sheet","Consumer Densities");
resourcedensitytable = readtable(filename,"Sheet","Resource Densities");
traits                 = table2array(traittable);
consumerdensity       = table2array(consumerdensitytable);
resourcedensity1      = resourcedensitytable.R1;
resourcedensity2      = resourcedensitytable.R2;

gridsize              = sqrt(length(consumerdensity));

%% Grid generation
h = 2*pi/gridsize;
x=1/gridsize.*(1:gridsize); % creating equally spaced points between [-pi,pi]
[xx,yy]=meshgrid(x); % Converting 1D mesh to a 2D mesh.

%% Laplacian Operator
col = [-pi^2/(3*h^2)-1/6  -(-1).^(1:gridsize-1)./(2*sin((1:gridsize-1)*h/2).^2)];
Lap = (2*pi)^2 * toeplitz(col); % The Laplacian operator in 1 dimension, domain = [0,1]
I = eye(gridsize);
Lap2D = kron(I,Lap)+kron(Lap,I); % The Laplacian operator in 2 Dimension.

%% Calculate the fitness landscape


% Here, we create a function that can help us find the mutant with trait
% value at the fitness maximum.

% Separating the constants into different variables

NN = gridsize^2;
alpha = parametertable.alpha;          % Scaled maximal per-capita growth.
beta  = parametertable.beta;           % Scaled resource-dependent growth parameter.
dC    = parametertable.deltaC;             % Scaled consumer diffusion rate.
a0    = parametertable.a0;             % Curvature parameter of the affinity trade-off 
b     = parametertable.b;              % Shape parameter of the affinity trade-off.
s     = parametertable.s;              % Control Parameter for the resource substitutability



%% Seperating the structure containing matrices into different variables

% We need to find the trait value with the highest positive dominant
% eigenvalue. So we need two variables to store the trait value and
% dominant eigenvalue.
Xk = [];
r = 201;
% At first, we create a for loop to collect all the dominant eigenvalues of
% the matrix F to see if any of them are greater than 0.
D = zeros(r,1); trait = linspace(-1,1,r); 
for i = 1:r
    
    
    Xr = trait(i);
    
    % First, we should define the resource affinity of the consumer towards
    % resource 1 and 2. This will help in reducing complexity in the final
    % equation. Much easier to make changes and make corrections if necessary.
   
    a1 = 1+ (b.*log(0.5.*(exp(-a0.*Xr)+1))-a0.*Xr)/(a0-b*log(0.5*(exp(-a0)+1))); % Affinity towards resource 1
    a2 = 1+ (b.*log(0.5.*(exp(+a0.*Xr)+1))+a0.*Xr)/(a0-b*log(0.5*(exp(-a0)+1))); % Affinity towards resource 2


    % Now that we have defined the affinities for both the resources, we can
    % define the growth rate G of the new mutant consumer under selection.
    if s==0
        T = sqrt(a1*resourcedensity1*a2.*resourcedensity2);
    else
        T = 0.5^(1/s).*((a1*resourcedensity1).^(s)+(a2*resourcedensity2).^(s)).^(1/s);
    end
    G = 1./(alpha + beta./T);

    % G in this form is a vector and cannot be added to or subtracted 
    % from a matrix. Inorder to define the matrix F, we need to
    % convert G into its appropriate matrix form. This can be achieved
    % by defining a matrix whose diagonal elements are the elements of G.
    Gdiag = diag(G);

    % m is a scalar and to convert it into its appropriate matrix from, we
    % have to multiply it with an identity matrix of appropriate
    % dimensions(N^2 x N^2).
    I=eye(NN);

    % Now, we have all the ingredients to create the matrix form of F. 
    F = Gdiag + (1-1/alpha).*I + dC*Lap2D;



    % Here, we compute the eigenvalue of the matrix F and sort it from
    % highest to lowest to pick out the dominant eigenvalue.  
    eigenvalues = sort(eig(F),ComparisonMethod="real"); 
    dominant_eigenvalue = eigenvalues(NN);
    D(i)=dominant_eigenvalue; 
end

Eigenval = zeros(size(traits));
for i = 1:length(traits)
    
    
    Xr = traits(i);
    
    % First, we should define the resource affinity of the consumer towards
    % resource 1 and 2. This will help in reducing complexity in the final
    % equation. Much easier to make changes and make corrections if necessary.
   
    a1 = 1+ (b.*log(0.5.*(exp(-a0.*Xr)+1))-a0.*Xr)/(a0-b*log(0.5*(exp(-a0)+1))); % Affinity towards resource 1
    a2 = 1+ (b.*log(0.5.*(exp(+a0.*Xr)+1))+a0.*Xr)/(a0-b*log(0.5*(exp(-a0)+1))); % Affinity towards resource 2


    % Now that we have defined the affinities for both the resources, we can
    % define the growth rate G of the new mutant consumer under selection.
    if s==0
        T = sqrt(a1*resourcedensity1*a2.*resourcedensity2);
    else
        T = 0.5^(1/s).*((a1*resourcedensity1).^(s)+(a2*resourcedensity2).^(s)).^(1/s);
    end
    G = 1./(alpha + beta./T);

    % G in this form is a vector and cannot be added to or subtracted 
    % from a matrix. Inorder to define the matrix F, we need to
    % convert G into its appropriate matrix form. This can be achieved
    % by defining a matrix whose diagonal elements are the elements of G.
    Gdiag = diag(G);

    % m is a scalar and to convert it into its appropriate matrix from, we
    % have to multiply it with an identity matrix of appropriate
    % dimensions(N^2 x N^2).
    I=eye(NN);

    % Now, we have all the ingredients to create the matrix form of F. 
    F = Gdiag + (1-1/alpha).*I + dC*Lap2D;



    % Here, we compute the eigenvalue of the matrix F and sort it from
    % highest to lowest to pick out the dominant eigenvalue.  
    eigenvalues = sort(eig(F),ComparisonMethod="real"); 
    dominant_eigenvalue = eigenvalues(NN);
    Eigenval(i)=dominant_eigenvalue; 
end

% Plotting the fitness landscape
fig = figure(1);
fig.Position = [100 100 800 650];
hold on
yline(0,LineWidth=3,Color='r',LineStyle='-')
plot(trait,D,LineWidth=6,Color=[107,174,214]/256)
scatter(traits,Eigenval,200,"black","filled")
ylim([-0.001,0.001])
%yticks([-0.01,0,0.01])
xlim([-1, 1])
xticks([-1 -0.5 0 0.5 1])
box(gca,'on');
set(gca,'FontSize',40,'LineWidth',5, 'FontWeight', 'bold');
hold off