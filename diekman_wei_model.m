function z = diekman_wei_model(~,u)

% equations for modified Sato et al 2010 model 

global g_ca g_k Iapp Eca Ek tauf taux
global dhalf dslope fhalf fslope xhalf xslope

V = u(1);
f = u(2);
x = u(3);

dinf = 1/(1+exp(-(V+dhalf)/dslope));
finf = 1/(1+exp((V+fhalf)/fslope));
xinf = 1/(1+exp(-(V+xhalf)/xslope));

Ica = g_ca*dinf*f*(V-Eca);
Ik = g_k*x*(V-Ek);

z(1) = Iapp - Ica - Ik;
z(2) = (finf - f)/tauf;
z(3) = (xinf - x)/taux;

z=z';



