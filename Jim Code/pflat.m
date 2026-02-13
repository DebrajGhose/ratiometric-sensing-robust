close all
clear all

x = [-1:0.01:1]; %position
c0 = 0.5;b = 0.1; 
cvals = c0*( 1 + b*x); %concentration gradient
kd = c0;

prx = cvals./(cvals+kd); %probability of active receptor p_R(x)

%plot(x,cvals);

nq = 10;
qvals = [1:nq]/nq;

figure(1);hold on;
for j = 1:nq
    q = qvals(j);
    ps = [0.0:0.01:1.0];
    pg = 1 - q./(1- (1-ps)*(1-q));
    plot(ps,pg,'LineWidth',2);
    
end
plot(ps,ps,'k');
xlabel('p_R');
ylabel('p_G^c');




nq = 5;
nvals = (2.^[0:nq-1])-1;

figure(3);hold on;
for j = 1:nq
    n = nvals(j);
    ps = [0.0:0.01:1.0];
    pg = ps*n./(1 + ps*n);
    plot(ps,pg,'LineWidth',2);
    %if n > 1
    %x = (1/sqrt(n)) - (1/n);
    %y=x*n/(1 + x*n);
    %plot(x,y,'ko','LineWidth',2);
    %end
end

%joff = 0.8*[0:10]/10;
%for j=1:11
%plot(ps,ps+joff(j),'k');
%end

plot(ps,ps,'k');
xlabel('p_R');
ylabel('p_G^c');
