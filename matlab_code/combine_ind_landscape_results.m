% Combining results of the simulations (evolutionary dynamics + trade-off + substitutability) 
codeDir   = fileparts(mfilename("fullpath"));   % .../matlab_code
parentDir = fileparts(codeDir);                 % shared parent folder
dataDir   = fullfile(parentDir,"affinity_tradeoff_substitutability");
% Getting data from land
for land = 4:4
land_index = land;
landpath = fullfile(dataDir,"resource_landscape_"+num2str(land_index),"ResourceLandscape_grid_size=32.xlsx");
land_param = readtable(landpath,'Sheet','Sheet2');

% Getting the results file from the folder path
folderpath = fullfile(dataDir,"resource_landscape_"+num2str(land_index),"Evolution","Evolutionary Stable State");
fin_filepath = folderpath+"\*.xlsx";
fin_filename = dir(fin_filepath);

% Getting the results file from the folder path for generalist
gen_folderpath = fullfile(dataDir,"resource_landscape_"+num2str(land_index),"Ecological Dynamics");
gen_filepath = gen_folderpath + "\*.xlsx";
gen_filename = dir(gen_filepath);

% Getting the results file from the folder path for homogeneous resources
hom_folderpath = fullfile(dataDir,"resource_landscape_0","Evolution","Evolutionary Stable State");

%%
% Creating a structure to store the results (combined)
trade_off            = zeros(972,1);
substitutability     = zeros(972,1);
consumer_biomass     = zeros(972,1);
consumer_biomass_gen = zeros(972,1);
consumer_biomass_hom = zeros(972,1);
resource_biomass     = zeros(972,2);
species_richness     = zeros(972,1);
rao_entropy          = zeros(972,1);
trait_width          = zeros(972,1);
gen_biomass          = zeros(972,1);
hom_biomass          = zeros(972,1);
resource_alignment   = zeros(972,1);

for num = 1:length(fin_filename)
    ind_filepath         = fullfile(folderpath,fin_filename(num).name);
    gen_filepath         = fullfile(gen_folderpath,gen_filename(num).name);
    hom_filepath         = fullfile(hom_folderpath,fin_filename(num).name);
    param                = readtable(ind_filepath,'Sheet','Parameters');
    consumer_density     = readtable(ind_filepath,'Sheet','Consumer Densities');
    resource_density     = readtable(ind_filepath,'Sheet','Resource Densities');
    inflow_rates         = readtable(ind_filepath,'Sheet','Inflow Rates');
    traits               = readtable(ind_filepath,'Sheet','Traits');
    gen_consumer_density = readtable(gen_filepath,'Sheet','Consumer Densities');
    hom_consumer_density = readtable(hom_filepath,'Sheet','Consumer Densities');
    ind_consumer_biomass = sum(consumer_density{1:end,:},1)/length(consumer_density.x1);
    sorted_consumer_biomass = sort(ind_consumer_biomass);
    
    trade_off(num)            = param.a0;
    substitutability(num)     = param.s;
    affinities                                     = trait_to_affinity(trade_off(num),param.b,traits.traits);


    consumer_biomass(num)     = sum(consumer_density{1:end,:},'all')/length(consumer_density.x1);
    gen_biomass(num)          = sum(gen_consumer_density{1:end,:},'all')/length(consumer_density.x1);
    hom_biomass(num)          = sum(hom_consumer_density{1:end,:},'all')/length(consumer_density.x1);
    consumer_biomass_gen(num) = (consumer_biomass(num)-gen_biomass(num))/gen_biomass(num);
    consumer_biomass_hom(num) = (consumer_biomass(num)-hom_biomass(num))/hom_biomass(num);
    resource_biomass(num,1)   = sum(resource_density.R1)/length(resource_density.R1);  
    resource_biomass(num,2)   = sum(resource_density.R2)/length(resource_density.R2);  
    species_richness(num)     = length(traits.traits);
    if length(traits.traits)~=1
        for index_i = 1:length(traits.traits)-1
            for index_j = index_i:length(traits.traits)
                rao_entropy(num) = rao_entropy(num)+sqrt((affinities(index_i,1)-affinities(index_j,1))^2+(affinities(index_i,2)-affinities(index_j,2))^2)*ind_consumer_biomass(index_i)*ind_consumer_biomass(index_j)/consumer_biomass(num)^2;
            end
        end
    end
    
    for trait_index = 1 : length(traits.traits)
        niche_resource_1      = sum(consumer_density{:,trait_index}.*inflow_rates.K1)/sum(consumer_density{:,trait_index});
        niche_resource_2      = sum(consumer_density{:,trait_index}.*inflow_rates.K2)/sum(consumer_density{:,trait_index});
        if traits.traits(trait_index) < 0
            resource_alignment(num)   = resource_alignment(num) + ind_consumer_biomass(trait_index)/consumer_biomass(num)*log(niche_resource_1/niche_resource_2);
        else
            resource_alignment(num)   = resource_alignment(num) + ind_consumer_biomass(trait_index)/consumer_biomass(num)*log(niche_resource_2/niche_resource_1);
        end
    end
    sorted_traits             = sort(traits.traits);
    sorted_affinities         = trait_to_affinity(trade_off(num),param.b,sorted_traits);
    trait_width(num)          = sqrt((sorted_affinities(1,1)-sorted_affinities(end,1))^2+(sorted_affinities(1,2)-sorted_affinities(end,2))^2);   
end

%% Creating the table
[trade_off,ind]      = sort(trade_off);
substitutability     = substitutability(ind);
consumer_biomass     = consumer_biomass(ind);
consumer_biomass_gen = consumer_biomass_gen(ind);
consumer_biomass_hom = consumer_biomass_hom(ind);
resource_biomass     = resource_biomass(ind,:);
species_richness     = species_richness(ind);
rao_entropy          = rao_entropy(ind);
trait_width          = trait_width(ind);
resource_alignment   = resource_alignment(ind);

[substitutability_mat,ind] = sort(reshape(substitutability,18,54),1,"descend");
consumer_abundance_mat     = reshape(consumer_biomass,18,54);
resource_1_abundance_mat   = reshape(resource_biomass(:,1),18,54);
resource_2_abundance_mat   = reshape(resource_biomass(:,2),18,54);
consumer_biomass_gen_mat   = reshape(consumer_biomass_gen,18,54);
consumer_biomass_hom_mat   = reshape(consumer_biomass_hom,18,54);
species_richness_mat       = reshape(species_richness,18,54);
rao_entropy_mat            = reshape(rao_entropy,18,54);
trait_width_mat            = reshape(trait_width,18,54);
resource_alignment_mat     = reshape(resource_alignment,18,54);

for col_ind = 1:54
    consumer_abundance_mat(:,col_ind)   = consumer_abundance_mat(ind(:,col_ind),col_ind);
    resource_1_abundance_mat(:,col_ind) = resource_1_abundance_mat(ind(:,col_ind),col_ind);
    resource_2_abundance_mat(:,col_ind) = resource_2_abundance_mat(ind(:,col_ind),col_ind);
    species_richness_mat(:,col_ind)     = species_richness_mat(ind(:,col_ind),col_ind);
    rao_entropy_mat(:,col_ind)          = rao_entropy_mat(ind(:,col_ind),col_ind);
    trait_width_mat(:,col_ind)          = trait_width_mat(ind(:,col_ind),col_ind);
    consumer_biomass_gen_mat(:,col_ind) = consumer_biomass_gen_mat(ind(:,col_ind),col_ind);
    consumer_biomass_hom_mat(:,col_ind) = consumer_biomass_hom_mat(ind(:,col_ind),col_ind);
    resource_alignment_mat(:,col_ind)   = resource_alignment_mat(ind(:,col_ind),col_ind);
end

%%
substitutability      = substitutability_mat(:);
consumer_biomass      = consumer_abundance_mat(:);
resource_biomass(:,1) = resource_1_abundance_mat(:);
resource_biomass(:,2) = resource_2_abundance_mat(:);
species_richness      = species_richness_mat(:);
rao_entropy           = rao_entropy_mat(:);
trait_width           = trait_width_mat(:);
consumer_biomass_gen  = consumer_biomass_gen_mat(:);
consumer_biomass_hom  = consumer_biomass_hom_mat(:);
resource_alignment    = resource_alignment_mat(:);

%% Filling the results structure
result = struct("corr_coef",[],"trade_off",[],"substitutability",[],"consumer_biomass",[],"resource_1_biomass",[],"resource_2_biomass",[],"consumer_biomass_gen",[],"consumer_biomass_hom",[],"evolved_species",[],"rao_entropy",[],"trait_width",[],"resource_alignment",[]);

result.corr_coef            = repmat(land_param.rho,972,1);
result.trade_off            = trade_off;
result.substitutability     = substitutability;
result.consumer_biomass     = consumer_biomass;
result.consumer_biomass_gen = consumer_biomass_gen;
result.consumer_biomass_hom = consumer_biomass_hom;
result.resource_1_biomass   = resource_biomass(:,1); 
result.resource_2_biomass   = resource_biomass(:,2); 
result.evolved_species      = species_richness;
result.rao_entropy          = rao_entropy;
result.trait_width          = trait_width;
result.resource_alignment   = resource_alignment;
T                           = struct2table(result);
writetable(T,fullfile(dataDir,"resource_landscape_"+num2str(land_index),"results_combined.xlsx"));
end

function affinity = trait_to_affinity(trd_off,b,traits)
    affinity = [1+(b.*log(0.5*(exp(-trd_off.*traits)+1))-trd_off.*traits)/(trd_off-b*log(0.5*(exp(-trd_off)+1))),1+(b.*log(0.5*(exp(trd_off.*traits)+1))+trd_off.*traits)/(trd_off-b*log(0.5*(exp(-trd_off)+1)))];
end