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
x(1)=1.5e-13/(200e-6)^2*3600; x(2)=1.9e-9/(200e-6)*3600; x(3)=5.8e-9/(200e-6)*3600;
Dd = 8e-4/(200e-6)^2*3600; 

[xpt,ytheor_up]=specliposome_uptake([x(1) x(2) Dd]);

tspan=[0,0.25/60,0.5/60];  
Ytheorup=[];
for met=1:13
    y0=ytheor_up(met,:); y0=y0';
    [xpt,ytheorup]=specliposome_clearance([x(1) 1000 Dd],y0,tspan);
    %Ytheorup(end+1,:)=ytheorup(end,:); 
    Ytheorup(end+1,:)=ytheorup(1,:); 
end

% Clearance
y0=ytheor_up(end,:); y0=y0';
tspan1=[0:0.5:4*24];
[xpt,ytheor_cl]=specliposome_clearance([x(1) x(3) Dd],y0,tspan1);
tspan=[0,0.25/60,0.5/60];  
Ytheorcl=[];
for met=1:length(tspan1)
    y0=ytheor_cl(met,:); y0=y0';
    [xpt,ytheorcl]=specliposome_clearance([x(1) 1000 Dd],y0,tspan);
    %Ytheorcl(end+1,:)=ytheorcl(end,:); 
    Ytheorcl(end+1,:)=ytheorcl(1,:); 
end

