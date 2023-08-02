function generate_figures(x)

% x is the parameter vector consisting of: x= [ Deff, hd, hcl ], i.e.:
% the effective diffusivity of antibodies, Deff and hd, hcl the mass transfer
% coefficients during uptake and clearance experiments, respectively. 

% First read uptake
Tup=readtable('Non-specific_Ab_Rituximab_Uptake and Clearance_Antibody.xlsx','Sheet','Sheet1');
% Read sheet with clearance results 
Tcl=readtable('Non-specific_Ab_Rituximab_Uptake and Clearance_Antibody.xlsx','Sheet','Sheet2');

 
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
[xpt,ytheor_up]=nonspecantibody_uptake([x(1) x(2)]);

i=6; % end of incubation time
y0=ytheor_up(i,:); y0=y0';
tspan=[0,0.25/60,0.5/60]; % Assume fishing lasts half a minute
[xpt,ytheor]=antibodies_clearance([x(1) 10000],y0,tspan); 

Y24 = ytheor(end,:);

subplot(2,2,1),errorbar(Tup.Var1,Tup.Average_5,Tup.StdDev_5),hold on
subplot(2,2,1),plot(xpt*max(Tup.Var1),Y24*0.06); axis([0 205 0 0.06])
xlabel('distance from spheroid center (um)'), ylabel('Ab concentration (uM)'), title('24 hrs incubation')
 
% Clearance simulation 
tspan=[0, 1, 2, 4];  % hrs
[xpt,ytheor_cl]=antibodies_clearance([x(1) x(3)],y0,tspan); 
% Fishing simulation
Ytheor=[];
for i=1:3
          y0=ytheor_cl(i,:); y0=y0';
          tspan=[0,0.25/60,0.5/60];
          [xpt,ytheor]=antibodies_clearance([x(1) 1000],y0,tspan);  
          Ytheor(:,end+1)=ytheor(end,:);
end

subplot(2,2,2),errorbar(Tup.Var1,Tcl.Average_1,Tcl.StdDev_1),hold on
subplot(2,2,2),plot(xpt*max(Tup.Var1),Ytheor(:,1)*0.06); axis([0 205 0 0.06])
xlabel('distance from spheroid center (um)'), ylabel('Ab concentration (uM)'), title('1 hr clearance')

subplot(2,2,3),errorbar(Tup.Var1,Tcl.Average_2,Tcl.StdDev_2),hold on
subplot(2,2,3),plot(xpt*max(Tup.Var1),Ytheor(:,2)*0.06); axis([0 205 0 0.06])
xlabel('distance from spheroid center (um)'), ylabel('Ab concentration (uM)'), title('2 hrs clearance')

subplot(2,2,4),errorbar(Tup.Var1,Tcl.Average_3,Tcl.StdDev_3),hold on
subplot(2,2,4),plot(xpt*max(Tup.Var1),Ytheor(:,3)*0.06); axis([0 205 0 0.06])
xlabel('distance from spheroid center (um)'), ylabel('Ab concentration (uM)'), title('4 hrs clearance')
