% script to reproduce Figure 2 from Diekman and Wei (2021)
% Circadian Rhythms of Early Afterdepolarizations and Ventricular Arrhythmias in a Cardiomyocyte Model
% Biophysical Journal, VOLUME 120, ISSUE 2, P319-333
% https://doi.org/10.1016/j.bpj.2020.11.2264

clear all; close all

% set parameters

global g_ca g_k Iapp Eca Ek tauf taux
global dhalf dslope fhalf fslope xhalf xslope

g_k = 0.05;

dhalf=7.3; dslope=8.6; fhalf=13.3; fslope=11.9;

xhalf=40; xslope=5;

Eca=60;
Ek=-80;
tauf=80;
taux=300;

Iapp=0;

inits1 = [-80 0.9963 0.0003354]; % f_inf(-80) and x_inf(-80)

options = odeset('AbsTol',1e-9,'RelTol',1e-9);

t0=0;
tf1=10;
tin1=t0:.1:tf1;

tf2=1000;
tin2=tf1:.1:tf2;


%% ZT 15 (Fig 2 blue trace)

g_ca = 0.15; 

[t1,u1] = ode15s('diekman_wei_model',tin1,inits1,options);

inits2=u1(end,:);
inits2(1)=0;

[t2,u2] = ode15s('diekman_wei_model',tin2,inits2,options);

t_zt15=[t1; t2];
V_zt15=[u1(:,1); u2(:,1)];

%% ZT 3 (Fig 2 red trace)

g_ca = 0.3; 

[t1,u1] = ode15s('diekman_wei_model',tin1,inits1,options);

inits2=u1(end,:);
inits2(1)=0;

[t2,u2] = ode15s('diekman_wei_model',tin2,inits2,options);

t_zt3=[t1; t2];
V_zt3=[u1(:,1); u2(:,1)];

%% make plot

set(0,'DefaultAxesFontSize',24)

figure(1)
hold on
plot(t_zt15,V_zt15,'b','linewidth',3)
plot(t_zt3,V_zt3,'r--','linewidth',3)
xlim([0 700])
ylim([-90 60])
set(gca,'XTick',0:100:700,'YTick',-80:40:40)
xlabel('t (ms)','Interpreter','latex')
ylabel('V (mV)','Interpreter','latex')
legend('ZT 15','ZT 3')
legend('boxoff')




