% number of nodes 
np = 100; 
% Discretize domain
xpt=zeros(np,1); 
dx=1/np; 
xpt(1)=dx/2; 
for i=2:np
    xpt(i)=xpt(i-1)+dx;
end

% Parameters for Ab simulation
Deff=6e-12/(300e-6)^2*3600; % effective diffusivity of Ab in the spheroids
hup=2.5e-10/(300e-6)*3600; % mass transfer coefficient during uptake
hcl=8e-7/(300e-6)*3600; % mass transfer coefficient during clearance
koff=1e-2*3600; % dissociation rate constant
KD=4.5; % dissociation equilbrium constant
kr=5.5e-6*3600; % internalization rate constant


% Uptake experiment
ti=0:0.5:24;
[xpt,ytheor_up]=specantibody_uptake([Deff hup koff koff/KD Rt Abbulk kr], ti);
y0=ytheor_up(end,:); y0=y0';
Ytheorup=[];
for met=1:length(ti)
    Ytheorup(end+1,:)=ytheor_up(met,:); 
end

% Clearance experiment
ti=[0:0.5:8];
j=1:np; 
[xpt,ytheor_cl]=antibodies_clearance([Deff hcl koff koff/KD Rt Abbulk kr],y0,xpt, ti); 
Ytheorcl=[];
for i=1:length(ti)
    Ytheorcl(end+1,:)=ytheor_cl(i,:);
end
    
