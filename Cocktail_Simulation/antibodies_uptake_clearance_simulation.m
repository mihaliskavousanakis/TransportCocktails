% number of nodes 
np = 100; 
% Discretize domain
xpt=zeros(np,1); 
dx=1/np; 
xpt(1)=dx/2; 
for i=2:np
    xpt(i)=xpt(i-1)+dx;
end

Deff=6e-12/(200e-6)^2*3600;
hup=2.5e-10/(200e-6)*3600;
hcl=2.5e-7/(200e-6)*3600;
koff=4e-3*3600;
KD=5;
kr=1.4e-5*3600;

ti=0:0.5:24;
[xpt,ytheor_up]=specantibody_uptake([Deff hup koff koff/KD Rt Abbulk kr], ti);
y0=ytheor_up(end,:); y0=y0';
tspan=[0,0.25/60,0.5/60];
Ytheorup=[];
for met=1:length(ti)
    ynew=ytheor_up(met,:); ynew=ynew';
    [xpt,ytheorup]=antibodies_clearance([Deff 10000 koff koff/KD Rt Abbulk kr],ynew,xpt, tspan);     
    Ytheorup(end+1,:)=ytheorup(1,:); 
end

ti=[0:0.5:4*24];
j=1:np; 
[xpt,ytheor_cl]=antibodies_clearance([Deff hcl koff koff/KD Rt Abbulk kr],y0,xpt, ti); 
Ytheorcl=[];
for i=1:length(ti)
    y0=ytheor_cl(i,:); y0=y0';
    tspan=[0,0.25/60,0.5/60];
    [xpt,ytheor]=antibodies_clearance([Deff 10000 koff koff/KD Rt Abbulk kr],y0,xpt, tspan);     
    Ytheorcl(end+1,:)=ytheor(1,:);
end
    
