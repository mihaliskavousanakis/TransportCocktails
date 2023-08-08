function [xpt,ytheor] = specliposome_uptake(x)

% input parameters 
DL = x(1); % liposome effective diffusivity
hd = x(2); % mass transfer coefficient during uptake experiments
Dd = x(3); % drug diffusivity


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
phiL=xpor(1)*xpt.^xpor(2)+xpor(3); 
% Porosity for drug
phiD= 0.7313;

Mass = sparse(26*np,26*np); 
i=1:np; i=i';

xpt2=xpt(i).*xpt(i);

% species L
for met=1:22
    j=26*(i-1)+met; 
    Mass = Mass +sparse(j,j,xpt2.*phiL,26*np,26*np); 
end
% Drug unspent
j=26*(i-1)+23; 
Mass = Mass + sparse(j,j,xpt2.*phiD,26*np,26*np);
% Drug spent
j=26*(i-1)+24; 
Mass = Mass + sparse(j,j,xpt2.*phiD,26*np,26*np);

% Drug unspent internalized
j=26*(i-1)+25; 
Mass = Mass + sparse(j,j,1,26*np,26*np);
% Drug spent internalized
j=26*(i-1)+26; 
Mass = Mass + sparse(j,j,1,26*np,26*np);
    
y0=zeros(26*np,1); % initial condition
    
options = odeset('Mass',Mass,'Jacobian',@(t,y) liposomes_jac(t,y,xpt,DL,hd,xpor,Dd,0.7,1));

% tspan: measured in hrs
tspan=0:0.5:6;  
[t,y] = ode15s(@(t,y) liposomes(t,y,xpt,DL,hd,xpor,Dd,0.7,1), tspan, y0, options);

for tmet=1:length(tspan)
    for met=1:22
        j=26*(i-1)+met;
        ytheor(tmet,j)=y(tmet,j).*phiL';
    end
    j=26*(i-1)+23; 
    ytheor(tmet,j)=y(tmet,j).*phiD';
    j=26*(i-1)+24; 
    ytheor(tmet,j)=y(tmet,j).*phiD';
    j=26*(i-1)+25; 
    ytheor(tmet,j)=y(tmet,j);
    j=26*(i-1)+26; 
    ytheor(tmet,j)=y(tmet,j);
end

end
