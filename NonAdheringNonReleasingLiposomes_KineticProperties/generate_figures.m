function generate_figures(x)

% x is the parameter vector consisting of: x= [ Deff, hd, hcl ], i.e.:
% the effective diffusivity of liposomes, Deff and hd,hcl the mass transfer
% coefficients during uptake and clearance experiments

T=readtable('Non-Adhering Liposomes_Uptake and Clearance_Final Concentration.xlsx');
 
% number of nodes 
np = 100; 
% Discretize domain
xpt=zeros(np,1); 
dx=1/np; 
xpt(1)=dx/2; 
for i=2:np
    xpt(i)=xpt(i-1)+dx;
end
              
% Uptake simulation
[xpt,ytheor_up]=nonspecliposome_uptake([x(1) x(2)]);

i=4; % end of incubation time
y0=ytheor_up(i,:); y0=y0';
tspan=[0,0.25/60,0.5/60]; % Assume fishing lasts half a minute
[xpt,ytheor]=liposomes_clearance([x(1) 1000],y0,tspan); 

Y6 = ytheor(end,:);

subplot(2,2,1),errorbar(T.Var1,T.Average_3,T.StdDev_3),hold on
subplot(2,2,1),plot(xpt*max(T.Var1),Y6*0.5); axis([0 205 0 0.5])
xlabel('distance from spheroid center (um)'), ylabel('Liposome concentration (mM)'), title('6 hrs incubation')
 
% Clearance simulation 
tspan=[0, 0.5, 1, 2];  % hrs
[xpt,ytheor_cl]=liposomes_clearance([x(1) x(3)],y0, tspan); 
% Fishing simulation
Ytheor=[];
for i=1:3
          y0=ytheor_cl(i,:); y0=y0';
          tspan=[0,0.25/60,0.5/60];
          [xpt,ytheor]=liposomes_clearance([x(1) 1000],y0,tspan);          
          Ytheor(:,end+1)=ytheor(end,:);
end

subplot(2,2,2),errorbar(T.Var1,T.Average_4,T.StdDev_4),hold on
subplot(2,2,2),plot(xpt*max(T.Var1),Ytheor(:,1)*0.5); axis([0 205 0 0.5])
xlabel('distance from spheroid center (um)'), ylabel('Liposome concentration (mM)'), title('0.5 hrs clearance')

subplot(2,2,3),errorbar(T.Var1,T.Average_5,T.StdDev_5),hold on
subplot(2,2,3),plot(xpt*max(T.Var1),Ytheor(:,2)*0.5); axis([0 205 0 0.5])
xlabel('distance from spheroid center (um)'), ylabel('Liposome concentration (mM)'), title('1 hr clearance')

subplot(2,2,4),errorbar(T.Var1,T.Average_6,T.StdDev_6),hold on
subplot(2,2,4),plot(xpt*max(T.Var1),Ytheor(:,3)*0.5); axis([0 205 0 0.5])
xlabel('distance from spheroid center (um)'), ylabel('Liposome concentration (mM)'), title('2 hrs clearance')
