% fileprefix: prefix that is part of the config file name
% directory: directory of the config file
% tstop: final time point of simulation
% dt: time step
% n_R: number of receptors
% n_phe_1: the lowest pheromone concentration in the linear gradient (nM)
% n_phe_2: the highest pheromone concentration in the linear gradient (nM)
% samplingrate: how many time steps between each record
% (here we record every 10 secs which is 2000 time steps)
% random seed: random number generator

function smoldyn_cfg(fileprefix,directory,tstop,dt,n_R,n_phe_1,n_phe_2,samplingrate,random_seed)

%% Basic file and parameter setup %%
% Name the file
cfg_name = [directory '/' fileprefix '.cfg']; 
xyz_name = [fileprefix '.xyz']; 
fid=fopen(cfg_name,'w');
fprintf(fid, 'random_seed %d\n',random_seed);

fprintf(fid, 'variable L = 4\n'); % length of the simulation domain

% receptor deactivation rate
kri = 0.1;
fprintf(fid, 'variable kri = %g\n',kri); 
fprintf(fid, '\n');

% Set up domain boundaries (um)
fprintf(fid, 'dim 3\n');
fprintf(fid, 'boundaries x -L L\n');
fprintf(fid, 'boundaries y -3.01 3.01\n');
fprintf(fid, 'boundaries z -3.01 3.01\n');
fprintf(fid, '\n');

% Set up species
fprintf(fid, 'species Ram Rim Ga Gi complex_Gi_Ra complex_Ga_Ri\n'); fprintf(fid, '\n');
% Define diffusion rates for species (um^2/s)
fprintf(fid, 'difc Ga(all) %g\n',0.002);
fprintf(fid, 'difc Gi(all) %g\n',0.002);

% Set up a cell with radius = 2.5 um
fprintf(fid, 'start_surface cell\n');
fprintf(fid, 'polygon both edge\n');
fprintf(fid, 'action both all(all) reflect\n');
fprintf(fid, 'panel sphere 0 0 0 2.5 10 10 panel_sphere\n');
fprintf(fid, 'end_surface\n');
fprintf(fid, '\n');

% Define the compartment outside the cell. Not essential to simulations
fprintf(fid, 'start_compartment full_domain\n');
fprintf(fid, 'surface cell\n');
fprintf(fid, 'point 0 0 0\n');
fprintf(fid, 'end_compartment\n');
fprintf(fid, '\n');

fprintf(fid, 'start_surface outer_walls_1\n');
fprintf(fid, 'action both all(all) absorb\n');
fprintf(fid, 'polygon both edge\n');
fprintf(fid, 'panel rect +x 3.5 3.01 -3.01 -6.02 6.02\n');
fprintf(fid, 'panel rect -x -3.5 -3.01 -3.01 6.02 6.02\n');
fprintf(fid, 'end_surface\n');
fprintf(fid, '\n');

fprintf(fid, 'start_surface outer_walls_2\n');
fprintf(fid, 'action both all(all) reflect\n');
fprintf(fid, 'polygon both edge\n');
fprintf(fid, 'panel rect +y -3.5 3.01 -3.01 7 6.02\n');
fprintf(fid, 'panel rect -y 3.5 -3.01 -3.01 -7 6.02\n');
fprintf(fid, 'panel rect +z -3.5 3.01 3.01 7 -6.02\n');
fprintf(fid, 'panel rect -z -3.5 3.01 -3.01 7 -6.02\n');
fprintf(fid, 'end_surface\n');
fprintf(fid, '\n');

%% Set up phermone gradients %%
% 11 gradient sections
cell_radius = 2.5; N = 11; interval = 2*cell_radius/N; 
for i = 1:N-1
    fprintf(fid, 'start_surface wall_%d\n',i);
    fprintf(fid, 'polygon both edge\n');
    x_position = -cell_radius + i*interval;
    fprintf(fid, 'panel rect +x %g 3.01 -3.01 -6.02 6.02\n',x_position);
    fprintf(fid, 'end_surface\n');
end
fprintf(fid, '\n');

for i = 1:N
    fprintf(fid, 'start_compartment cmpt_%d\n',i);
    if i == 1 % For the first and last section we just need to add one boundary, as the other boundary is just the cell surface
        fprintf(fid, 'surface wall_%d\n',i);
        fprintf(fid, 'surface cell\n');
        wall_position_1 = -cell_radius + i*interval;
        wall_position_2 = -cell_radius;
        cmpt_position = (wall_position_1 + wall_position_2)/2;
        cmpt_height_1 = sqrt(cell_radius^2-(cell_radius-abs(wall_position_1))^2);
        cmpt_height_2 = 0;
        mid_height = (cmpt_height_1 + cmpt_height_2)/2;
        fprintf(fid, 'point %g 0 0\n',cmpt_position);
        fprintf(fid, 'point %g 0 %g\n',cmpt_position,mid_height);
        fprintf(fid, 'point %g 0 %g\n',cmpt_position,-mid_height);
        fprintf(fid, 'point %g 0 %g\n',cmpt_position,mid_height/2);
        fprintf(fid, 'point %g 0 %g\n',cmpt_position,-mid_height/2);
        fprintf(fid, 'end_compartment\n');
        fprintf(fid, '\n');
    elseif i == N
        fprintf(fid, 'surface wall_%d\n',i-1);
        fprintf(fid, 'surface cell\n');
        wall_position_1 = -cell_radius + (i-1)*interval;
        wall_position_2 = cell_radius;
        cmpt_position = (wall_position_1 + wall_position_2)/2;
        cmpt_height_1 = sqrt(cell_radius^2-(cell_radius-abs(wall_position_1))^2);
        cmpt_height_2 = 0;
        mid_height = (cmpt_height_1 + cmpt_height_2)/2;
        fprintf(fid, 'point %g 0 0\n',cmpt_position);
        fprintf(fid, 'point %g 0 %g\n',cmpt_position,mid_height/2);
        fprintf(fid, 'point %g 0 %g\n',cmpt_position,-mid_height/2);
        fprintf(fid, 'end_compartment\n');
    else
        fprintf(fid, 'surface wall_%d\n',i-1);
        fprintf(fid, 'surface wall_%d\n',i);
        fprintf(fid, 'surface cell\n');
        wall_position_1 = -cell_radius + i*interval;
        wall_position_2 = -cell_radius + (i-1)*interval;
        cmpt_position = (wall_position_1 + wall_position_2)/2;
        cmpt_height_1 = sqrt(cell_radius^2-(cell_radius-abs(wall_position_1))^2);
        cmpt_height_2 = sqrt(cell_radius^2-(cell_radius-abs(wall_position_2))^2);
        mid_height = (cmpt_height_1 + cmpt_height_2)/2;
        fprintf(fid, 'point %g 0 0\n',cmpt_position);
        fprintf(fid, 'point %g 0 %g\n',cmpt_position,mid_height/2);
        fprintf(fid, 'point %g 0 %g\n',cmpt_position,-mid_height/2);
        fprintf(fid, 'end_compartment\n');
        fprintf(fid, '\n');
    end
    diff = (n_phe_2 - n_phe_1)/(N-1);
    phe_conc(i) = n_phe_1 + (i-1)*diff; % pheromone concentration in each section
    activation_rate = kri/6*phe_conc(i); % calculate receptor activation rate. 
                                         % Kd = 6 nM, receptor activation rate 
                                         % = inactivation rate / Kd * pheromone conentration in the current section
    % Write the reaction Ri (inactive receptor) -> Ra (active receptor) with the activation rate in each section
    fprintf(fid, 'reaction compartment=cmpt_%g Receptor_activation_%d Rim(down) -> Ram(down) %g\n',i,i,activation_rate);
    fprintf(fid, '\n');
end

% Define molecule list for speeding up simulations in Smoldyn
fprintf(fid, 'molecule_lists list1 list2 list3 list4\n');
fprintf(fid, 'mol_list Rim list1\n');
fprintf(fid, 'mol_list Ram list2\n');
fprintf(fid, 'mol_list Ga list3\n');
fprintf(fid, 'mol_list Gi list4\n');

% Receptor activation: Ri (inactive receptor) -> Ra (active receptor)
fprintf(fid, 'reaction Receptor_deactivation Ram(down) -> Rim(down) kri\n');

% G protein activation: Gi (inactive G protein) + Ra (active receptor) ->
% Ga (active G protein) + Ra (active receptor)
fprintf(fid, 'reaction Gprotein_activation_1 Gi(down) + Ram(down) -> complex_Gi_Ra(down)\n');
fprintf(fid, 'reaction Gprotein_activation_2 complex_Gi_Ra(down) -> Ram(down) + Ga(down)\n');
fprintf(fid, 'reaction_probability Gprotein_activation_1 1\n');
fprintf(fid, 'binding_radius Gprotein_activation_1 0.005\n'); % interaction radius, 
                                                              % if reactants are within the radius they will react (bind)
fprintf(fid, 'reaction_probability Gprotein_activation_2 1\n');
fprintf(fid, 'product_placement Gprotein_activation_2 unbindrad 0.005+0.00005\n'); % unbinding radius, 
                                                                                   % added a very small distance to avoid repeated binding

% G protein inactivation: Ga (active G protein) + Ri (inactive receptor) ->
% Gi (inactive G protein) + Ri (inactive receptor)
fprintf(fid, 'reaction Gprotein_deactivation_1 Ga(down) + Rim(down) -> complex_Ga_Ri(down)\n');
fprintf(fid, 'reaction Gprotein_deactivation_2 complex_Ga_Ri(down) -> Rim(down) + Gi(down)\n');
fprintf(fid, 'reaction_probability Gprotein_deactivation_1 1\n');
fprintf(fid, 'binding_radius Gprotein_deactivation_1 0.005\n');
fprintf(fid, 'reaction_probability Gprotein_deactivation_2 1\n');
fprintf(fid, 'product_placement Gprotein_deactivation_2 unbindrad 0.005+0.00005\n');

%% Simulation setup %%
fprintf(fid, 'time_start 0\n'); % start time = 0
fprintf(fid, 'time_stop %g\n',tstop); % stop time
fprintf(fid, 'time_step %g\n',dt); % time step
fprintf(fid, 'surface_mol %d Rim(down) cell sphere panel_sphere\n',n_R); % add how many inactive receptors on the cell surface at t=0
fprintf(fid, 'surface_mol %d Gi(down) cell sphere panel_sphere\n',2500); % add how many inactive G proteins on the cell surface at t=0

% Record the poitions of molecules every 10 secs
fprintf(fid, 'output_files %s\n', xyz_name);
fprintf(fid, 'cmd I 1 %d %d molpos Ram(all) %s\n', tstop/dt, samplingrate, xyz_name);
fprintf(fid, 'cmd I 1 %d %d molpos Rim(all) %s\n', tstop/dt, samplingrate, xyz_name);
fprintf(fid, 'cmd I 1 %d %d molpos Ga(all) %s\n', tstop/dt, samplingrate, xyz_name);
fprintf(fid, 'cmd I 1 %d %d molpos Gi(all) %s\n', tstop/dt, samplingrate, xyz_name);
end

