clear all
close all
NR=1000;
NG=1000;
ep=0.1;
nrf = 1000;
Rfvals = binornd(NR,0.5+ep,nrf,1);
Rbvals = binornd(NR,0.5-ep,nrf,1);

%% Define custom colors
color_gold  = [255, 195, 57] / 255;   % #FFC339
color_cyan  = [11, 176, 213] / 255;   % #0BB0D5
color_green = [0, 224, 112] / 255;    % #00E070
color_teal  = [11, 128, 131] / 255;   % #0B8083

%% Histogram bin settings
binWidth_fig2 = 0.005;   % common bin width for Figure 2
binWidth_fig3 = 0.005;   % common bin width for Figure 3

%ratiometric model
%
alpha = 1.0;
actF = 0.5 + alpha*( (Rfvals/NR) - 0.5);
actB = 0.5 + alpha*( (Rbvals/NR) - 0.5);
Gfvals = zeros(1,nrf);
Gbvals = zeros(1,nrf);
for i=1:nrf
	Gfvals(i) = binornd(NG,actF(i));
	Gbvals(i) = binornd(NG,actB(i));
end
Gdifs = Gfvals-Gbvals;  %difference between front and back
figure(1);
scatter(Rfvals/NR,Gfvals/NG,20,color_gold,'filled','MarkerFaceAlpha',0.5);
hold on;
scatter(Rbvals/NR,Gbvals/NG,20,color_cyan,'filled','MarkerFaceAlpha',0.5);
xlabel('Active Receptor Fraction')
ylabel('Active G Protein Fraction')
title('Empirical Joint Distribution')
figure(2);
h1=histogram(Gfvals/NG,'BinWidth',binWidth_fig2,'Normalization','pdf','FaceColor',color_gold,'FaceAlpha',0.5);
hold on;
h2=histogram(Gbvals/NG,'BinWidth',binWidth_fig2,'Normalization','pdf','FaceColor',color_cyan,'FaceAlpha',0.5);
xlabel('Active G Protein Fraction')
title('Histogram')
cv2r = var(Gfvals)/(mean(Gfvals)^2)
figure(3);
hd1 = histogram(Gdifs/NG,'BinWidth',binWidth_fig3,'Normalization','pdf');
hd1.FaceColor = color_cyan;
hd1.FaceAlpha = 0.5;
hold on;
%
%classical model
%
alpha = 0.50;
actF = 0.5 + alpha*( (Rfvals/NR) - 0.5);
actB = 0.5 + alpha*( (Rbvals/NR) - 0.5);
Gfvals = zeros(1,nrf);
Gbvals = zeros(1,nrf);
for i=1:nrf
	Gfvals(i) = binornd(NG,actF(i));
	Gbvals(i) = binornd(NG,actB(i));
end
Gdifs = Gfvals-Gbvals;
figure(1);
scatter(Rfvals/NR,Gfvals/NG,20,color_green,'filled','MarkerFaceAlpha',0.5);
hold on;
scatter(Rbvals/NR,Gbvals/NG,20,color_teal,'filled','MarkerFaceAlpha',0.5);
xv = [0.30:0.001:0.7];
plot(xv,xv,'Color',color_gold,'LineWidth',1.5);
nl = xv./(xv + 0.5);
plot(xv,nl,'Color',color_cyan,'LineWidth',1.5);
lnl = 0.5 + 0.5*(xv - 0.5);
plot(xv,lnl,'Color',color_teal,'LineWidth',1.5);
legend('Ratiometric Front','Ratiometric Back','Classical Front','Classical Back',...
       'y = x','Nonlinear','Linear','Location','best');
figure(2);
h3=histogram(Gfvals/NG,'BinWidth',binWidth_fig2,'Normalization','pdf','FaceColor',color_green,'FaceAlpha',0.5);
hold on;
h4=histogram(Gbvals/NG,'BinWidth',binWidth_fig2,'Normalization','pdf','FaceColor',color_teal,'FaceAlpha',0.5);
cv2c = var(Gfvals)/(mean(Gfvals)^2)
legend('Ratiometric Front','Ratiometric Back','Classical Front','Classical Back','Location','best');
figure(3);
hd2 = histogram(Gdifs/NG,'BinWidth',binWidth_fig3,'Normalization','pdf');
hd2.FaceColor = color_gold;
hd2.FaceAlpha = 0.5;
hold on;
title('Histogram of G Protein gradient')
legend('Ratiometric','Classical','Location','best');