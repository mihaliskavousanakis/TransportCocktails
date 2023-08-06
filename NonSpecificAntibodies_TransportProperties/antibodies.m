function dydt=antibodies(t,y,xpt,Deff,hd,phi,cinf,xpor)

np=length(y); 
dydt=zeros(np,1); 
dx=xpt(2)-xpt(1); 

phi=@(r) xpor(1)*r.^xpor(2)+xpor(3);

% Symmetry boundary 
i=1;

% West node
xw=0; 
fw=0; 
%phiw=phi;
phiw=phi(xw);

% East node
xe=(xpt(i)+xpt(i+1))/2;
fe=(y(i+1)-y(i))/dx; 
%phie=phi;
phie=phi(xe);
dydt(i)=Deff/dx*(xe^2*phie*fe-xw^2*phiw*fe);

% inner nodes 
for i=2:np-1
        
    % West node
    xw = (xpt(i)+xpt(i-1))/2; 
    %phiw=phi;
    phiw=phi(xw); 
    % Total flux
    fw = (y(i)-y(i-1))/dx; 
    
    % East node
    xe=(xpt(i)+xpt(i+1))/2;
    %phie=phi;
    phie=phi(xe); 
    % Total flux
    fe = (y(i+1)-y(i))/dx;
     
    
    dydt(i)=Deff/dx*(xe^2*phie*fe-xw^2*phiw*fw);
    
end

% Right boundary - spheroid radius
i=np; 
% West node
xw = (xpt(i)+xpt(i-1))/2; 
%phiw=phi;
phiw=phi(xw); 
% Total flux
fw = (y(i)-y(i-1))/dx; 

% East node
xe = 1;
%phie=phi;
phie=phi(xe); 

if (hd<1000)
ce=(2*Deff*phie*y(i)/dx+hd*cinf)/(2*Deff/dx*phie+hd); 
else
    ce=cinf;
end

%ce=cinf;
fe = (ce-y(i))/(dx/2);

dydt(i)=Deff/dx*(xe^2*phie*fe-xw^2*phiw*fw);
