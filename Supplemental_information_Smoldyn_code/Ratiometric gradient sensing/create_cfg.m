% Create the directory that stores the simulations
directory = 'Smoldyn_simulations';
filebase = '231201.8'; % file version control, can change to any identifiers

mkdir(directory);

% End point of the simulation. Units in seconds.
tstop = 20000;

% Generate 10 realizations
random_seeds = 1:10;
dt = 0.005; % time step
samplingrate = round(10/dt); % how many time steps to record molecule positions (we record every 10 secs)
n_phe_2 = [8]; % highest pheromone concentration in the linear gradient
n_phe_1 = [4]; % lowest pheromone concentration in the linear gradient

n_R = [375]; % Total number of receptors. 

% Write the smoldyn file
fid=fopen(sprintf('%s/run.sh',directory),'w');
fprintf(fid,'#!/bin/bash\n\n');

for j = random_seeds
    for m = 1:numel(n_R)
        for h = 1:numel(n_phe_1)
            curr_seed  = j;
            curr_fileprefix = sprintf('%s-n_R_%g-n_phe_1_%g-n_phe_2_%g-seeds_%d',filebase,n_R(m),n_phe_1(h),n_phe_2(h),j);
            smoldyn_cfg(curr_fileprefix,directory,tstop,dt,n_R(m),n_phe_1(h),n_phe_2(h),samplingrate,j);
            cfg_name = sprintf('%s.cfg',curr_fileprefix);

            % optional: you can use a HPC cluster
            %fprintf(fid,'sbatch -p general -N 1 -J Smoldyn -t 168:00:00 --mem=2g --wrap="smoldyn %s"\n',cfg_name);
        end
    end
end
fclose(fid);