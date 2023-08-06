function dydt=antibodies(t,y,xpt,Deff,hd,koff,kon,Rt,Abbulk,kr,phi,cinf,xpor)

np=length(xpt); 
dydt=zeros(3*np,1); 
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
ia = 3*(i-1)+1; 
itheta = ia+1; 
iint=itheta+1; 
iap1 = ia+3; 
fe=(y(iap1)-y(ia))/dx; 
phie=phi(xe); 
dydt(ia)=Deff/dx*(xe^2*phie*fe-xw^2*phiw*fe)-kon*xpt(i).^2*Rt*y(ia).*(1-y(itheta))+koff*xpt(i)^2*Rt/Abbulk*y(itheta);
dydt(itheta)=kon*Abbulk*y(ia).*(1-y(itheta))-koff*y(itheta)-kr*y(itheta);
dydt(iint)=kr*Rt/Abbulk*y(itheta);

% inner nodes 
for i=2:np-1
    
    ia=3*(i-1)+1; 
    iap1=ia+3; 
    iam1=ia-3; 
    itheta=ia+1;
    iint=itheta+1; 
    % West node
    xw = (xpt(i)+xpt(i-1))/2; 
    phiw=phi(xw); 
    % Total flux
    fw = (y(ia)-y(iam1))/dx; 
    
    % East node
    xe=(xpt(i)+xpt(i+1))/2;
    phie=phi(xe); 
    % Total flux
    fe = (y(iap1)-y(ia))/dx;
     
    
    dydt(ia)=Deff/dx*(xe^2*phie*fe-xw^2*phiw*fw)-kon*xpt(i).^2*Rt*y(ia).*(1-y(itheta))+koff*xpt(i)^2*Rt/Abbulk*y(itheta);
    dydt(itheta)=kon*Abbulk*y(ia).*(1-y(itheta))-koff*y(itheta)-kr*y(itheta);
    dydt(iint)=kr*Rt/Abbulk*y(itheta);
end

% Right boundary - spheroid radius
i=np; 
ia=3*(i-1)+1; 
iam1=ia-3; 
itheta=ia+1; 
iint=itheta+1;
% West node
xw = (xpt(i)+xpt(i-1))/2;
phiw=phi(xw); 
% Total flux
fw = (y(ia)-y(iam1))/dx; 

% East node
xe = 1;
phie=phi(xe); 

if hd<1000
ce=(2*Deff*phie*y(ia)/dx+hd*cinf)/(2*Deff/dx*phie+hd); 
else
ce=cinf;
end
fe = (ce-y(ia))/(dx/2);

dydt(ia)=Deff/dx*(xe^2*phie*fe-xw^2*phiw*fw)-kon*xpt(i).^2*Rt*y(ia).*(1-y(itheta))+koff*xpt(i)^2*Rt/Abbulk*y(itheta);
dydt(itheta)=kon*Abbulk*y(ia).*(1-y(itheta))-koff*y(itheta)-kr*y(itheta);
dydt(iint)=kr*Rt/Abbulk*y(itheta);
