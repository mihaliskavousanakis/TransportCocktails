function [xpt,ytheor] = antibodies_clearance(x,y0,R,tspan)

% input parameters 
Deff = x(1); % effective diffusivity 
hd=x(2); % mass transfer coefficient 
koff = x(3); % dissociation rate constant
kon = x(4); % association rate constant
Rt = x(5); % Ab receptors concentration 
Abbulk=x(6); % Ab bulk concentration
kr=x(7); % internalization rate constant 


% number of nodes 
np = 100; 

% Discretize domain
xpt=zeros(np,1); 
dx=1/np; 

xpt(1)=dx/2; 
for i=2:np
    xpt(i)=xpt(i-1)+dx;
end


% Load porosity functional
load porosity_Ab.mat
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

% initial condition
y0(6*(i-1)+1)=y0(6*(i-1)+1)./phi;    
y0(6*(i-1)+2)=y0(6*(i-1)+2)./phi;    
options = odeset('Mass',Mass);

% tspan: measured in hrs

[t,y] = ode15s(@(t,y) antibodies(t,y,xpt,Deff,hd,koff,kon,Rt,Abbulk,kr,phi,0,xpor), tspan, y0, options);
ytheor=y; j=1:np;
for i=1:length(tspan)
    ytheor(i,6*(j-1)+1)=ytheor(i,6*(j-1)+1).*phi';
    ytheor(i,6*(j-1)+2)=ytheor(i,6*(j-1)+2).*phi';
end

end