function [xpt,ytheor] = nonspecliposome_uptake(x)

% input parameters 
Deff = x(1); % effective diffusivity um^2/hrs
hd=x(2); % mass transfer coefficient um/hrs


% number of nodes 
np = 100; 

% Discretize domain
xpt=zeros(np,1); 
dx=1/np; 

xpt(1)=dx/2; 
for i=2:np
    xpt(i)=xpt(i-1)+dx;
end


% Load porosity coefficients
load porosity
phi=xpor(1)*xpt.^xpor(2)+xpor(3); 

Mass = sparse(np,np); 
i=1:np; 

xpt2=xpt(i).*xpt(i);
Mass = Mass + sparse(i,i,xpt2.*phi,np,np);



y0=zeros(np,1); % initial condition
    
options = odeset('Mass',Mass);

% tspan: measured in hrs
tspan=[0, 1, 2, 4, 6];  
[t,y] = ode15s(@(t,y) liposomes(t,y,xpt,Deff,hd,phi,1,xpor), tspan, y0, options);

ytheor=y(2:end,:); 
for i=1:length(tspan)-1
    ytheor(i,:)=ytheor(i,:).*phi';
end

end