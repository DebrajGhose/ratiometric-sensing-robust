%% ratiometric receptor with G-protein simulation

clear all
close all

load('seed.mat')
rng(seed);

%% simulation parameters

totaltime = 600; %seconds
cellradius = 2.5; %in microns

kri = 0.1; % receptor deactivation rate
kD = 6; %dissociation constant in nM
kra = kri/kD;
kpgi = 1; % probability of g protein deactivation when active g protein encounters inactive receptor 
kgi = 0.1; % base G protein deactivation rate

reaction_radius = 0.005; %or you can think of this as being receptor size

NR = 6000; %number of receptors
Rprot = zeros(NR,4);% each row is a receptor protein, first Rp(:,1:3) are x y z positions, Rp(:,4) are states.
Rprot(:,1:3) = scatter_on_sphere(NR,cellradius);
newrecstates = zeros(NR,1); %vector to store new rec states

Pher = Pheromone_gradient(8,4,Rprot,cellradius); %create pheromone gradient (high, low, receptor positions,...)

Rprot(:,4) = rand(NR,1) < (Pher./(Pher+kD)); % set receptor states at steady state

NG = 2500; %number of G-proteins
Gprot = zeros(NG,4);% each row is a g protein, first Gp(:,1:3) are x y z positions, Gp(:,4) are states.
Gprot(:,1:3) = scatter_on_sphere(NG,cellradius); % randomly distribute G proteins
Gprot(:,4) = (rand(NG,1)>0.5);% assign G protein states

Dm = 0.002; %micron^2/s
r6Dt = sqrt(2*3*0.002*0.00005); %I want to keep r6Dt consistent with what it was for Dm = 0.002 and dt = 0.00005. 2*3 because diffusion is in 3D and only later projected into 2D; this is probably fine on short distance scales, but for accurate diffusion, use code from Ghose and Lew. MBoC 2020.
dt = round(r6Dt^2/(6*Dm),10); %seconds. I am changing dt to keep r6Dt the same as what it was when dt = 0.00005. Rounding is to prevent numerical errors

% NOTE: Make sure that sqrt(2*4*Dm*dt) is several fold less than
% reaction_radius. The "4" is because diffusion happens essentially 2
% dimension on the spherical cell.

%% data collection matrices

capturefreq = 1; %how many seconds per capture

allRprot = zeros(NR,4,totaltime/capturefreq);

allGprot = zeros(NG,4,totaltime/capturefreq);

%% create spatial partition or receptors

RKdtree = KDTreeSearcher(Rprot(:,1:3)); % generate model to perform partitioned search later

%% simulate

for tt = 1:totaltime/dt


%% diffuse g proteins on membrane

Gprot = diffuse_membrane(Gprot,r6Dt,cellradius);

%% receptors react to field of pheromone

dieroll = rand(NR,1);

newrecstates(Rprot(:,4)==0) = (dieroll(Rprot(:,4)==0) < Pher(Rprot(:,4)==0)*kra*dt);  %activate
newrecstates(Rprot(:,4)==1) = ~(dieroll(Rprot(:,4)==1) < kri*dt); %deactivate

Rprot(:,4) = newrecstates;

%% using range search to find all particles within some distance

%range search uses kdtree searching to generate indices of neighbors
%within some distance

R_close_IDs = rangesearch(RKdtree,Gprot(:,1:3),reaction_radius); %identify receptors that are close to g proteins

for ii = 1:size(Gprot,1)
if ~isempty(R_close_IDs{ii})
%The chance that any G protein becomes active is determined by
%the number of active/total receptors within the reaction
%radius. If there only active receptors, the g protein will
%become active, if there are only inactive receptors, the g
%protein will become inactive.

prob_react = sum(Rprot(R_close_IDs{ii},4))/numel(R_close_IDs{ii}); % probability of reaction to occur

if rand()<prob_react %toss a die for reaction
Gprot(ii,4)=1;
else
Gprot(ii,4)=0;
end

end
end

%% collect data every second

if mod(tt*dt,capturefreq) == 0

allRprot(:,:,tt*dt) = Rprot;

allGprot(:,:,tt*dt) = Gprot;

end

%% visualize simulations
%{
if mod(tt*dt,1)==0
scatter3( Gprot(:,1), Gprot(:,2), Gprot(:,3),30, Gprot(:,4),'filled' );
hold on
view([0 90])
axis equal
axis off
title(num2str(dt*tt))
drawnow
hold off
end
%}


end

save 'output.mat' %just save the entire workspace

function [Pher] = Pheromone_gradient(high,low,Rprot,cellradius) % generate pheromone gradient

slope = (high-low)/(cellradius*2); offset = (high+low)/2;

Pher = Rprot(:,1)*slope + offset;

end

function [molloc] = diffuse_membrane(molloc,r6Dt,r) %naively diffuse particles on a sphere. For more accurate random walk on a sphere, use (time-consuming) function from https://github.com/DebrajGhose/Actin_driven_patch_movement/blob/master/Agent%20based%20model%20and%20analysis/MeasurePersistence_improved.m

mover = normrnd(0,r6Dt, [size(molloc,1),1] ); %amount of movement
diffm = scatter_on_sphere(size(molloc,1),mover); % acquire vector that is basically small movement in a random direction in 3 dimensions

molloc(:,1:3) = molloc(:,1:3)+diffm;

[azm,elem,rm] = cart2sph(molloc(:,1),molloc(:,2),molloc(:,3)); %convert to spherical coodinates
rm(:)=r; %confine within radius
[x,y,z]=sph2cart(azm,elem,rm);
molloc(:,1:3) = [x,y,z];

end

function [particle_location] = scatter_on_sphere(N,cellradius) %N is number of particles, r is radius of sphere

az = 2*pi*rand(N,1); %generate random coordinates
ele = asin(2*rand(N,1)-1);
radius = cellradius.*ones(N,1);
[particle_location(:,1),particle_location(:,2),particle_location(:,3)] = sph2cart(az,ele,radius);

end

