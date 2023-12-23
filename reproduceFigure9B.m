clear all; close('all')

set(0,'DefaultAxesFontSize',24)

StimProtocol=2;

% This script is based on Peter Hammer's code "Spiral waves in monodomain
% reaction-diffusion-model" provided on the MATLAB Central File Exchange:
% https://uk.mathworks.com/matlabcentral/fileexchange/22492-spiral-waves-in-monodomain-reaction-diffusion-model

% Matlab implementation a monodomain reaction-diffusion model in 2-D. 
% The model equations are a variant of the Fitzhugh-Nagumo equations
% modified to simulate the cardiac action potential. The progression
% of the two normalized state variables, membrane voltage (v) and recovery 
% (r), is computed across a 128 x 128 spatial domain and across time. This 
% function simulates spiral waves, which are hypothesized to underlie 
% reentrant tachycardia. The spiral waves can be initiated by two different 
% cardiac pacing methods: 
%
% (1) two-point stimulation where a point stimulus is delivered in the 
% center of the domain followed by another point stimulus on the partially
% refractory wake of the first wave of excitation.
%
% (2) cross-field stimulation where a stimulus is applied to the left 
% domain boundary causing a plane wave. As this wave travels across the
% domain, a second stimulus is applied to the bottom boundary of the domain.
%
% This function accepts only one input argument, StimProtocol, which can 
% be either the numerical values '1' (for two-point stimulation) or '2' 
% (for cross-field stimulation). As the simulation runs, the activation 
% state of the individual units comprising the domain is mapped to color 
% and plotted in a figure window. A count of time steps is displayed at 
% the top of the plot along with the values of v and r for the upper left 
% element of the domain. 
% 
% Model equations are solved using a finite difference method for spatial 
% derivatives and explicit Euler integration for time derivatives. Newman
% boundary conditions are assumed. Model parameters are taken from two 
% journal articles: [1] Rogers JM et al. "A collocation-Galerkin finite 
% element model of cardiac action potential propagation". IEEE TBME;41:
% 743-757,1994. [2] Pertsov AM et al. "Spiral waves of excitation underlie 
% reentrant activity in isolated cardiac muscle". Circulation Research;72:
% 631-650, 1993. 
%
% Following the simulated spiral waves, a movie (AVI) is generated, and 
% the user is given the option to save the movie to disk. One simulation 
% takes about 160 seconds on a 2.33 GHz Intel Dual Core 64-bit laptop 
% (3.3 GB RAM).
%
% This function was written by Peter E. Hammer (hammer@usa.com) in 
% October 2006 as part of an academic project associated with PhD studies 
% in biomedical engineering. Have fun!

ncols=128;                               % Number of columns in domain
nrows=128;                               % Number of rows in domain
dur=60000;                               % Number of time steps
h=0.01;                                 % Grid size
h2=h^2;
dt=0.1;                                     % Time step                                
Iex=500;                                % Amplitude of external current
mu=1.0;                                  % Anisotropy
Gx=0.000025; Gy=Gx/mu;                  % Conductances               
v=-80*ones(nrows,ncols);                % Initialize voltage array
fgate=ones(nrows,ncols);                % Initialize gatingg variable array
xgate=zeros(nrows,ncols);               % Initialize gatingg variable array

% we save the potential at several points for plotting
% this is the index to plot the potential at (center of domain)
rowInd=round(nrows/2);
colInd=round(ncols/2);

% Parameters for Diekman and Wei model
Iapp=0;
Eca=60;
Ek=-80;
%g_ca=0.15;
g_ca=0.3;
g_k=0.1;
%g_k=0.05;

gk_mat=g_k*ones(nrows,ncols);
gk_mat(:,1:102)=0.05;                 % heterogeniety in potassium channel conductance

dhalf=7.3; dslope=8.6; fhalf=13.3; fslope=11.9;
xhalf = 40; xslope = 5;

tauf=80;
taux=300;

% Set initial stim current and pattern
iex=zeros(nrows,ncols);
if StimProtocol==1
    iex(62:67,62:67)=Iex;
else
    iex(:,1)=Iex;
end

clim=[-80 60];
% Setup image
ih=imagesc(v,clim); %set(ih,'cdatamapping','direct')
colormap(hot); axis image off; th=title('');
colorbar

set(gcf, 'PaperPosition', [0.25 2.5 8.0 6.0]);


n=0;                 % Counter for time loop
k=0;                 % Counter for movie frames
done=0;              % Flag for while loop

n1e=20;              % Step at which to end 1st stimulus
switch StimProtocol
    case 1           % Two-point stimulation
        n2b=3800;    % Step at which to begin 2nd stimulus
        n2e=3900;    % Step at which to end 2nd stimulus
    case 2           % Cross-field stimulation
        n2b=8100;     % Step at which to begin 2nd stimulus
        n2e=8130;     % Step at which to end 2nd stimulus
end

n3b=10000;
n3e=10030;

while ~done          % Time loop
    
    if n == n1e      % End 1st stimulus
        iex=zeros(nrows,ncols);
    end
    
    if n == n2b      % Begin 2nd stimulus
        switch StimProtocol
            case 1
                iex(62:67,49:54)=Iex;
            case 2
                iex(end,:)=Iex;
        end
    end
    
    if n == n2e      % End 2nd stimulus
        iex=zeros(nrows,ncols);
    end
    
    % Create padded v matrix to incorporate Newman boundary conditions 
    vv=[[0 v(2,:) 0];[v(:,2) v v(:,end-1)];[0 v(end-1,:) 0]];
    
    % Update v
    vxx=(vv(2:end-1,1:end-2) + vv(2:end-1,3:end) -2*v)/h2; 
    vyy=(vv(1:end-2,2:end-1) + vv(3:end,2:end-1) -2*v)/h2;
        
    dinf = 1./(1+exp(-(v+dhalf)/dslope));
    finf = 1./(1+exp((v+fhalf)/fslope));
    
    Ica = g_ca*dinf.*fgate.*(v-Eca);
        
    xinf = 1./(1+exp(-(v+xhalf)/xslope));
    %Ik = g_k*xgate.*(v-Ek);
    Ik = gk_mat.*xgate.*(v-Ek);
    
    dvdt = Iapp - Ica - Ik + iex + Gx*vxx + Gy*vyy;

    v_new = v + dvdt*dt;
    
    dfdt=(finf-fgate)/tauf;
    fgate=fgate+dfdt*dt;
    dxdt=(xinf-xgate)/taux;
    xgate=xgate+dxdt*dt;
    v=v_new; clear v_new
    
    VsaveMiddleCol(n+1,:)=v(:,64);
    VsaveMiddleRow(n+1,:)=v(64,:);
    
    % Update image and text 
    set(ih,'cdata',v);
    set(th,'string',sprintf('%d  %0.2f   %0.2f    %0.2f',n,v(rowInd,colInd),fgate(rowInd,colInd),xgate(rowInd,colInd)))
    drawnow
        
    n=n+1;
    done=(n > dur);

end


%% make plot similar to Figure 7D

t=(0:(n-1))*dt;

figure(2)
hold on
plot(t,VsaveMiddleRow(:,1),'b','linewidth',3)
plot(t,VsaveMiddleRow(:,64),'r','linewidth',3)
plot(t,VsaveMiddleRow(:,128),'k','linewidth',3)
ylim([-90 60])
set(gca,'XTick',0:1000:5000,'YTick',[-80:40:80])
xlabel('$t$ (ms)','interpreter','latex')
ylabel('$V$ (mV)','interpreter','latex')



