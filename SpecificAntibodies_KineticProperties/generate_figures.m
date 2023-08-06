function generate_figures(x)

% x is the parameter vector consisting of: x= [Pup, Pcl,  koff, KD, kint ], i.e.:
% the mass transfer coefficient of specific Ab during uptake
% the mass transfer coefficient of specific Ab during clearance
% the dissociation rate constant, koff
% the dissociation equilibrium constant, KD
% the internalization rate constant, kint.

% Load data
T=readtable('Specific_ab_Trastuzumab_Uptake and Clearance_Antibody.xlsx','Sheet','Sheet2');

 
% number of nodes 
np = 100; 
% Discretize domain
xpt=zeros(np,1); 
dx=1/np; 
xpt(1)=dx/2; 
for i=2:np
    xpt(i)=xpt(i-1)+dx;
end
% Load parameters from non-specific Ab simulations 
load nonspecific_data.mat 
Deff=xnew(1)*(205/200)^2; 
% 205 um is the radius of spheroids where fitting was performed for transport properties of non-specific Ab and 200 um is the radius of spheroids 
% where the fitting for specific Ab is now performed
     
% load values for concentration of Ab (in nM) in the solution (Abbulk)
% and Rt the concentration of antibody receptors (in nM)
load experimental_Abbulk_Rt              

hup=x(1); 
hcl=x(2);
koff=x(3); 
KD=x(4); 
kr=x(5);

% Uptake experiments
[xpt,ytheor_up]=specantibody_uptake([Deff hup koff koff/KD Rt Abbulk kr], [0 1 2 4 6 12 16 24]);
% Pick the last solution that corresponds to 24 hrs incubation and set
% it as initial condition for the clearance experiments 
y0=ytheor_up(end,:); y0=y0';

% Clearance experiment
tspan=[0, 1, 2, 4]; % 0 corresponds to 24 hrs incubation 
[xpt,ytheor_cl]=antibodies_clearance([Deff hcl koff koff/KD Rt Abbulk kr],y0, tspan); 
% Fishing simulation
Ytheor=zeros(length(xpt),4);
j=1:np;
for i=1:4
          y0=ytheor_cl(i,:); y0=y0';
          tspan=[0,0.25/60,0.5/60];
          [xpt,ytheor]=antibodies_clearance([Deff 10000 koff koff/KD Rt Abbulk kr],y0, tspan); 
          y1=ytheor(end,3*(j-1)+1);
         y2=ytheor(end,3*(j-1)+2);
         y3=ytheor(end,3*(j-1)+3);
         Ytheor(:,i)=y1+y2*Rt/Abbulk+y3; % this includes all forms of antibodies (interstitium, bound, and internalized forms)
end


subplot(2,2,1),errorbar(T.Var1,T.Avg,T.StdDev),hold on
subplot(2,2,1),plot(xpt*max(T.Var1),Ytheor(:,1)*0.06); axis([0 205 0 0.06])
xlabel('distance from spheroid center (um)'), ylabel('Ab concentration (uM)'), title('24 hrs incubation')
 

subplot(2,2,2),errorbar(T.Var1,T.Avg_2,T.StdDev_2),hold on
subplot(2,2,2),plot(xpt*max(T.Var1),Ytheor(:,2)*0.06); axis([0 205 0 0.06])
xlabel('distance from spheroid center (um)'), ylabel('Ab concentration (uM)'), title('1 hr clearance')

subplot(2,2,3),errorbar(T.Var1,T.Avg_3,T.StdDev_3),hold on
subplot(2,2,3),plot(xpt*max(T.Var1),Ytheor(:,3)*0.06); axis([0 205 0 0.06])
xlabel('distance from spheroid center (um)'), ylabel('Ab concentration (uM)'), title('2 hrs clearance')

subplot(2,2,4),errorbar(T.Var1,T.Avg_4,T.StdDev_4),hold on
subplot(2,2,4),plot(xpt*max(T.Var1),Ytheor(:,4)*0.06); axis([0 205 0 0.06])
xlabel('distance from spheroid center (um)'), ylabel('Ab concentration (uM)'), title('4 hrs clearance')
