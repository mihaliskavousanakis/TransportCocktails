% number of nodes 
np = 100; 
% Discretize domain
xpt=zeros(np,1); 
dx=1/np; 
xpt(1)=dx/2; 
for i=2:np
    xpt(i)=xpt(i-1)+dx;
end
              
% Setup parameters for simulation
x(1)=1.5e-13/(300e-6)^2*3600;  % effective diffusivity of liposomes
x(2)=1.9e-9/(300e-6)*3600; % mass transfer coefficient during uptake
x(3)=5.8e-9/(300e-6)*3600; % mass transfer coefficient during clearance
Dd =2.e-11/(300e-6)^2*3600; % effective drug diffusivity 
hdup=1.7e-8/(300e-6)*3600; % drug mass transfer coefficient during uptake
hdcl=2.9e-8/(300e-6)*3600; % drug mass transfer coefficient during clearance

% uptake simulation
[xpt,ytheor_up]=specliposome_uptake([x(1) x(2) Dd hdup]);
Ytheorup=[];
for met=1:13
    Ytheorup(end+1,:)=ytheor_up(met,:); 
end

% Clearance simulation
y0=ytheor_up(end,:); y0=y0';
tspan1=[0:0.5:8];
[xpt,ytheor_cl]=specliposome_clearance([x(1) x(3) Dd hdcl],y0,tspan1);
Ytheorcl=[];
for met=1:length(tspan1)
    Ytheorcl(end+1,:)=ytheor_cl(met,:); 
end

