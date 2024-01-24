%%%
%%% plotUCflowTimeSeries.m
%%%
%%% Plots a time series of simulated undercurrent flow speed.
%%%

%%% Load pre-computed undercurrent flow time series
load('undercurrent_LWF_tau0.05_Nc4_Yc25_UCflow.mat');
tt_wind = tt;
uu_wind = u_UC;
load('undercurrent_RF_F0.75_Nc4_Yc25_UCflow.mat');
tt_rand = tt;
uu_rand = u_UC;
load('undercurrent_LWF_tau0.05_Nc4_Yc25_Fbaro0.025_UCflow.mat');
tt_baro = tt;
uu_baro = u_UC;

fontsize = 18;
linewidth = 1.5;
t1day = 86400;
t1year = t1day*365;

figure(5);
clf;
set(gcf,'Position',[479         513        1364         716]);
plot(tt_wind/t1day,uu_wind,'LineWidth',linewidth);
hold on;
plot(tt_rand/t1day,uu_rand,'LineWidth',linewidth);
plot(tt_baro/t1day,uu_baro,'LineWidth',linewidth);
hold off;
legend('Wind-forced','Randomly-forced','Wind- and APF-forced','Location','SouthEast');
set(gca,'FontSize',fontsize);
set(gca,'XLim',[0 365]);
set(gca,'Position',[0.07 0.1 0.9 0.85]);
xlabel('Time (days)');
ylabel('Undercurrent speed (m/s)');

figure(6);
clf;
set(gcf,'Position',[479         513        1364         716]);
plot(tt_wind/t1year,uu_wind,'LineWidth',linewidth);
hold on;
plot(tt_rand/t1year,uu_rand,'LineWidth',linewidth);
plot(tt_baro/t1year,uu_baro,'LineWidth',linewidth);
hold off;
legend('Wind-forced','Randomly-forced','Wind- and APF-forced','Location','SouthEast');
set(gca,'FontSize',fontsize);
set(gca,'XLim',[0 10]);
set(gca,'Position',[0.07 0.1 0.9 0.85]);
xlabel('Time (years)');
ylabel('Undercurrent speed (m/s)');