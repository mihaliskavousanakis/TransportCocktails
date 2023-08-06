function [xpt,ytheor] = specantibody_uptake(x,tspan)

% input parameters 
Deff = x(1); % effective diffusivity 
hd=x(2); % mass transfer coefficient 
koff = x(3); % dissociation rate constant
kon = x(4);  % association rate constant (KD=koff/kon)
Rt = x(5); % Concentration of antibody receptors in the spheroid (in nM)
Abbulk=x(6); % Concentration of Ab in the bulk solution (in nM)
kr=x(7);  % Internalization rate constant of bound Ab


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
load porosity_new
phi=xpor(1)*xpt.^xpor(2)+xpor(3);

Mass = sparse(3*np,3*np); 
i=1:np; 

xpt2=xpt(i).*xpt(i);
Mass = Mass + sparse(3*(i-1)+1,3*(i-1)+1,xpt2.*phi,3*np,3*np);
Mass = Mass + sparse(3*(i-1)+2,3*(i-1)+2,1,3*np,3*np);
Mass = Mass + sparse(3*(i-1)+3,3*(i-1)+3,1,3*np,3*np);


y0=zeros(3*np,1); % initial condition
    
options = odeset('Mass',Mass);

% tspan: measured in hrs
[t,y] = ode15s(@(t,y) antibodies(t,y,xpt,Deff,hd,koff,kon,Rt,Abbulk,kr,phi,1,xpor), tspan, y0, options);


ytheor=y(2:end,:);  j=1:np; 
for i=1:length(tspan)-1
    ytheor(i,3*(j-1)+1)=ytheor(i,3*(j-1)+1).*phi';
end

end