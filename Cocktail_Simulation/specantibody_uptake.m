function [xpt,ytheor] = specantibody_uptake(x,tspan)

% input parameters 
Deff = x(1); % effective diffusivity um^2/hrs
hd=x(2); % mass transfer coefficient um/hrs
koff = x(3); 
kon = x(4); 
Rt = x(5); 
Abbulk=x(6); 
kr=x(7); 


% number of nodes 
np = 100; 

% Discretize domain
xpt=zeros(np,1); 
dx=1/np; 

xpt(1)=dx/2; 
for i=2:np
    xpt(i)=xpt(i-1)+dx;
end


% Load porosity value
%load porosity_value
load porosity_Ab
phi=xpor(1)*xpt.^xpor(2)+xpor(3);

Mass = sparse(6*np,6*np); 
i=1:np; 

xpt2=xpt(i).*xpt(i);
Mass = Mass + sparse(6*(i-1)+1,6*(i-1)+1,xpt2.*phi,6*np,6*np);
Mass = Mass + sparse(6*(i-1)+2,6*(i-1)+2,xpt2.*phi,6*np,6*np);
Mass = Mass + sparse(6*(i-1)+3,6*(i-1)+3,1,6*np,6*np);
Mass = Mass + sparse(6*(i-1)+4,6*(i-1)+4,1,6*np,6*np);
Mass = Mass + sparse(6*(i-1)+5,6*(i-1)+5,1,6*np,6*np);
Mass = Mass + sparse(6*(i-1)+6,6*(i-1)+6,1,6*np,6*np);


y0=zeros(6*np,1); % initial condition
    
%options = odeset('Mass',Mass,'RelTol',1e-3,'AbsTol',1e-6,'Jacobian',@eleven_kinds_plusDox_Jacobian);
options = odeset('Mass',Mass);

% tspan: measured in hrs
[t,y] = ode15s(@(t,y) antibodies(t,y,xpt,Deff,hd,koff,kon,Rt,Abbulk,kr,phi,1,xpor), tspan, y0, options);


ytheor=y(1:end,:);  j=1:np; 
for i=1:length(tspan)
    %ytheor(i,3*(j-1)+1)=ytheor(i,3*(j-1)+1)*phi;
    ytheor(i,6*(j-1)+1)=ytheor(i,6*(j-1)+1).*phi';
    ytheor(i,6*(j-1)+2)=ytheor(i,6*(j-1)+2).*phi';
end

end