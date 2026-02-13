%%Set Parameters
NR = 375; %number of receptors
Dm = 0.002; %micron^2/sec
cellradius = 2.5; %micron
koff = 0.008;  %rate constant For G protein deactivation, 1/sec.
reaction_radius = 0.005; %micron, reaction radius

Receptor_area_fraction = NR*(reaction_radius/cellradius)^2

NG = 10000;  %number of G proteins

%%Distribute NR Receptors uniformly on sphere. Locations are fixed.
Rprot = zeros(NR,4);%each row is a receptor protein, first Rp(:,1:3) are x y z positions, Rp(:,4) is a flag variable.
Rprot(:,1:3) = scatter_on_sphere2(NR,cellradius);


%%Create spatial partition For receptors
RKdtree = KDTreeSearcher(Rprot(:,1:3)); 


%%Generate independent Toff times For each Gprotein, exponentially distributed with mean 1/koff
Gtoffs = -(1/koff)*log(1-rand(1,NG));  

%%initial locations For NR Gproteins
Gprot = zeros(NG,3);
Gprot = scatter_on_sphere(NG,cellradius);

Gencounters = zeros(NG,1); %record number of encounters
Genctimes = zeros(NG,1); %record time to (first?) encounter
GprotLast = zeros(NG,3); %%final positions For Gproteins
Rnearest = -1.0 + zeros(NG,1); %%distance to first encounter location
RnearestEuc = -1.0 + zeros(NG,1); %%distance to first encounter location
totalBMsteps = zeros(NG,1);
Tmax = 1000;

%%loop through G proteins.
%%For each protein: evolve until toff has occured, compute number of distinct receptors encountered, and their locations.
for gg=1:NG

    cgloc = Gprot(gg,1:3); %current Gprotein location
    Toff = Gtoffs(gg); %time of state reset; exponentially distributed
    Tcum = 0; %accumulated time

    %we will remove rows after each encounter, to exclude proteins already visited.    
    RprotTemp = Rprot;

    %%compute distinct receptor encounters up to time Toff 
    while (Tcum < Toff)

        tlimit = Toff - Tcum; %time limit For computing next encounter

        if (Gencounters(gg,1) < 1)
            tlimit = Tmax; %allow longer computation time to compute first encounter.
        end

        RKdtree = KDTreeSearcher(RprotTemp(:,1:3)); 

        %%%GetNextEncounter to compute the following:
        %%Tenc, time until next encounter, 
        %%loc, location of next encounter,
        %%idx, index of the receptor encountered.
        %%bsteps, debugging variable, number of brownian steps

        ycurrent = cgloc; %%current G protein location
        steps=0;
        R_nearest_loc = [0,0,0];
        R_nearest_id = -1;
        tcount = 0;
        arrived = false; %whether or not Gprotein has arrived at a receptor
        ind = -1; %this will stay -1 if no encounter occurs.

        %keep stepping until arrived at a receptor or until time expires 
        while (arrived == false) && (tcount < tlimit)
            %%get distance to the nearest receptor, from current G protein location
            R_nearest_id = knnsearch(RKdtree,ycurrent);  
            R_nearest_loc = RprotTemp(R_nearest_id,1:3);
            rstep = spheredist(R_nearest_loc,ycurrent,cellradius); %%sphere distance to nearest receptor
            rstep = min(rstep,(pi/20)*cellradius);  %%limit the step size, since this is based on local linearization of shere
            if (rstep <= reaction_radius)
                %%encounter has occured.
                arrived = true;
                ind = R_nearest_id;
            else
                %%not arrived yet...take a Brownian step out to distance rstep
                [ynew,tstep] = BmotionStep(ycurrent,rstep,cellradius,Dm);
                ycurrent = ynew;  %%update G protein location
                tcount = tcount + tstep; %%update time
                steps = steps+1;
            end
        end
        
        loc = ycurrent;
        Tenc = tcount; 
        Tcum = Tcum + Tenc; %update accumulated time

         if (arrived==true)&&(Gencounters(gg,1) < 1)
             Genctimes(gg,1) = Tenc; %record time of first encounter
             Rnearest(gg,1) = spheredist(R_nearest_loc,cgloc,cellradius); %record distance traveled to first encounter
              RnearestEuc(gg,1) = eucldist(R_nearest_loc,cgloc); %record euclidean distance traveled to first encounter
         end

        totalBMsteps(gg) = totalBMsteps(gg) + steps;


        if (arrived==true) && (Tcum < Toff)
            %%this encounter happens before Toff, so we make a record of it
            cgloc = loc;  %%update current g protein location
            Gencounters(gg,1) = Gencounters(gg,1) + 1; %%increment record of encounters
            %Rprot(ind,4) = 1; %%flag most recent receptor encountered.
            RprotTemp(ind,:) = [];  %remove this receptor from list
        end
    end


    GprotLast(gg,1:3) = cgloc;
    %%record location last encounter (if any) (not same as position at time Toff)

end



%%ANALYSIS, DATA DISPLAY

meanNumberEncounters = mean(Gencounters)   %average number of encounters per G protein

kvals = [0:1:9];

geop = 1/(meanNumberEncounters + 1);  %set parameter for geometric distribution, for comparison
geoDist = geop*((1-geop).^kvals); %Geometric(p) distribution
edges = [-0.5:1:9.5]; %bins for histogram counts
hg = histcounts(Gencounters,edges);
probZeroEncounters = hg(1)/NG  %frequency of no encounter

%% Define custom colors
color_gold  = [255, 195, 57] / 255;   % #FFC339
color_cyan  = [11, 176, 213] / 255;   % #0BB0D5
color_green = [0, 224, 112] / 255;    % #00E070

kvalsoff = kvals + 0.2;
figure(1); hold on;
b1 = bar(kvalsoff, geoDist, 'FaceColor', color_gold);
b2 = bar(kvals, hg/NG, 'FaceColor', color_cyan);
xlabel('Distinct Encounters');
ylabel('Frequency');
legend([b1, b2], 'Geometric Distribution', 'Simulation', 'Location', 'best');

MeanFirstEncTime = mean(Genctimes(:,1))

figure(2); hold on;
bins = [0:30:600];
hg2 = histcounts(Genctimes,bins);
b3 = bar(bins(1:end-1), hg2/NG, 'FaceColor', color_green, 'FaceAlpha', 0.4);
xlabel('Time to first encounter (sec)');
ylabel('Frequency');
legend(b3, 'First Encounter Time', 'Location', 'best');






%%Given current Gprotein location y, and max_time interval
%%Computes time to next receptor encounter, location of next encounter, and index of receptor encountered.
function [ enc_time, enc_loc, enc_ind, bsteps] = getNextEncounter(receptors,kdtree,y,max_time,cellrad,dm,rstar)

    ycurrent = y; %%current G protein location
    steps=0;
    R_nearest_loc = y;
    R_nearest_id = 1;
    tcount = 0;
    arrived = false;

    while (arrived == false) && (tcount < max_time)
        %%get distance to the nearest receptor, from current G protein location
        R_nearest_id = knnsearch(kdtree,ycurrent);  
        R_nearest_loc = receptors(R_nearest_id,1:3);
        rstep = spheredist(R_nearest_loc,ycurrent,cellrad); %%sphere distance to nearest receptor

        if (rstep <= rstar)
            %%encounter has occured.
            arrived = true;
        else
            %%not arrived yet...take a Brownian step out to distance rstep
            [ynew,tstep] = BmotionStep(ycurrent,rstep,cellrad,dm);
            ycurrent = ynew;  %%update G protein location
            tcount = tcount + tstep; %%update time
            steps = steps+1;
        end
    end

    enc_ind = R_nearest_id;
    enc_loc = ycurrent;
    enc_time = tcount;
    bsteps = steps;
end




%%%% compute geodesic distance between two points on sphere of radius cellrad
%%posA and posB are in euclidean coordinates
function [sphdist] = spheredist(posA,posB,cellrad)

    rA = sqrt(posA(1)^2 + posA(2)^2 + posA(3)^2);
    rB = sqrt(posB(1)^2 + posB(2)^2 + posB(3)^2);

    %%inner product of normalized vectors
    scaledot = (posA(1)*posB(1) + posA(2)*posB(2) + posA(3)*posB(3))/(rA*rB);

    scaledot = min(scaledot,1);
    scaledot = max(scaledot,-1);

    %angular distance between the two vectors    
    theta = acos(scaledot);  

    sphdist = cellrad*theta;
end



%%%% compute euclidean distance between two points im R3
%%posA and posB are in euclidean coordinates
function [edist] = eucldist(posA,posB)
     
    edist = sqrt((posA(1) - posB(1))^2 + (posA(2) - posB(2))^2 + (posA(3) - posB(3))^2);

end



%%%%Evolve brownian motion on sphere, starting from initial position pos until is has moved distance stepsize
%%pos should be a point on the sphere of radius cellradius
%%stepsize should be in the interval [0,(pi/4)*cellradius]
%%the time calculation is based on linearization of the sphere about y, so it assumes (stepsize/cellradius) should be small.
%%returns location at time Bmotion reaches target distance, and the time elapsed until that event occurs 
function [loc,time] = BmotionStep(pos, stepsize, cellrad, dv)

    angle = stepsize/cellrad;  %angular distance to move on sphere
    vsame = pos*cos(angle); %new position vector projected onto the vector pos

    posnorm = pos/(cellrad); %unit vector in direction of the vector pos.

    %%Now we compute the projection of the new position vector onto subspace orthogonal to pos vector
    %%get two basis vectors for the subspace orthogonal to posnorm vector
    v1 = [0,0,0];
    v2 = [0,0,0];
    for j=1:3
        if (posnorm(j) < 3/4)  %there must be at least one component of posnorm that is less than 3/4
            k = 1 + mod(j,3); %k,q are the other two indices
            q = 1 + mod(j+1,3);  

            %project e_j onto orthogonal subspace: v1 = (I - posnorm*posnorm^T)e_j,   v1 is orthog to pos
            v1(j) = 1 - posnorm(j)*posnorm(j); 
            v1(k) = - posnorm(k)*posnorm(j); 
            v1(q) = - posnorm(q)*posnorm(j); 
            sz = sqrt((v1(1))^2 + (v1(2))^2 + (v1(3))^2); %length of v1

            %normalize v1 to be a unit vector.
            v1 = v1/sz;

            %compute cross product of v1 and posnorm; this gives a vector v2 orthogonal to both v1 and posnorm
            v2(1) = posnorm(2)*v1(3) - posnorm(3)*v1(2);
            v2(2) = - (posnorm(1)*v1(3) - posnorm(3)*v1(1));
            v2(3) = posnorm(1)*v1(2) - posnorm(2)*v1(1);  

            break;
        end    
    end
    
    %choose random rotation angle in the plane orthogonal to pos
    uval = 2*pi*rand(1,1);
    vperp = (v1*cos(uval) + v2*sin(uval))*sin(angle)*cellrad;  %component of new position in plane orthogonal to pos

    %%compute new location, uniform on circle (in sphere metric) of radius stepsize, centered at pos.
    loc = vsame + vperp;

    %%loc should have length equal to cellradius, but we'll normalize to make sure.
    %sz = sqrt((loc(1))^2 + (loc(2))^2 + (loc(3))^2);    
    %loc = cellrad*loc/sz;
    
    %%generate time traveled.  
    %T1 = time required For a standard Bmotion to go distance 1, in two dimensions (euclidean space)
    T1 = - (1/4)*log(1-rand(1,1)); %exponential random variable with mean 1/4; mean time to exit unit disc in 2d is 1/4.  
    %%We approximate this by an exponential with the correct mean; the exact distribution is not exponential.

    time = (stepsize^2/dv)*T1; %Brownian scaling in two dimensions (euclidean).
end





%%Debraj's scatter function.
function [particle_location] = scatter_on_sphere(N,cellradius) %N is number of particles, r is radius of sphere

    az = 2*pi*rand(N,1); %generate random coordinates
    ele = asin(2*rand(N,1)-1);
    radius = cellradius.*ones(N,1);
    [particle_location(:,1),particle_location(:,2),particle_location(:,3)] = sph2cart(az,ele,radius);
end


%%alternative scatter function.
function [particle_location] = scatter_on_sphere2(N,cellradius) %N is number of particles, r is radius of sphere
    v = randn(N,3);
    lengths = sqrt(v(:,1).^2 + v(:,2).^2 + v(:,3).^2);
    for i=1:N
       particle_location(i,1)=cellradius*v(i,1)/lengths(i);
       particle_location(i,2)=cellradius*v(i,2)/lengths(i);
       particle_location(i,3)=cellradius*v(i,3)/lengths(i);
    end
end