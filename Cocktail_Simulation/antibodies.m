function dydt=antibodies(t,y,xpt,Deff,hd,koff,kon,Rt,Abbulk,kr,phi,cinf,xpor)

% kspent
thalf=10*24; %hrs
kspent = log(2)/thalf;

np=length(xpt); 
dydt=zeros(3*np,1); 
dx=xpt(2)-xpt(1); 

% Porosity for Antibodies
phi=@(r) xpor(1)*r.^xpor(2)+xpor(3);

% Symmetry boundary 
i=1;

% West node
xw=0; 
fw=0; 
phiw=phi(xw); 

% East node
xe=(xpt(i)+xpt(i+1))/2;
ia = 6*(i-1)+1; 
iaB=ia+1;
itheta = iaB+1; 
ithetaB=itheta+1;
iint=ithetaB+1; 
iintB=iint+1;
iap1 = ia+6;
iaBp1 = iaB+6;
fe=(y(iap1)-y(ia))/dx; 
feB=(y(iaBp1)-y(iaB))/dx; 
phie=phi(xe); 
dydt(ia)=Deff/dx*(xe^2*phie*fe-xw^2*phiw*fe)-kon*xpt(i).^2*Rt*y(ia).*(1-y(itheta)-y(ithetaB))+koff*xpt(i)^2*Rt/Abbulk*y(itheta)-kspent*y(ia)*xpt(i)^2; % unspent interstitium Ab
dydt(iaB)=Deff/dx*(xe^2*phie*feB-xw^2*phiw*feB)-kon*xpt(i).^2*Rt*y(iaB).*(1-y(itheta)-y(ithetaB))+koff*xpt(i)^2*Rt/Abbulk*y(ithetaB)+kspent*y(ia)*xpt(i)^2; % spent interstitium Ab
dydt(itheta)=kon*Abbulk*y(ia).*(1-y(itheta)-y(ithetaB))-koff*y(itheta)-kr*y(itheta)-kspent*y(itheta); % unspent bound Ab
dydt(ithetaB)=kon*Abbulk*y(iaB).*(1-y(ithetaB)-y(itheta))-koff*y(ithetaB)-kr*y(ithetaB)+kspent*y(itheta); % spent bound Ab
dydt(iint)=kr*Rt/Abbulk*y(itheta)-kspent*y(iint); % unspent internalized Ab
dydt(iintB)=kr*Rt/Abbulk*y(ithetaB)+kspent*y(iint); % spent internalized Ab

% inner nodes 
for i=2:np-1
    
    ia=6*(i-1)+1; 
    iaB=ia+1;
    iap1=ia+6;
    iaBp1=iaB+6;
    iam1=ia-6; 
    iaBm1=iaB-6; 
    itheta=iaB+1;
    ithetaB=itheta+1;
    iint=ithetaB+1; 
    iintB=iint+1;
    % West node
    xw = (xpt(i)+xpt(i-1))/2; 
    phiw=phi(xw); 
    % Total flux
    fw = (y(ia)-y(iam1))/dx; 
    fwB = (y(iaB)-y(iaBm1))/dx; 
    
    % East node
    xe=(xpt(i)+xpt(i+1))/2;
    phie=phi(xe); 
    % Total flux
    fe = (y(iap1)-y(ia))/dx;
    feB = (y(iaBp1)-y(iaB))/dx;
     
    
    dydt(ia)=Deff/dx*(xe^2*phie*fe-xw^2*phiw*fw)-kon*xpt(i).^2*Rt*y(ia).*(1-y(itheta)-y(ithetaB))+koff*xpt(i)^2*Rt/Abbulk*y(itheta)-kspent*y(ia)*xpt(i)^2; % unspent interstitium Ab
    dydt(iaB)=Deff/dx*(xe^2*phie*feB-xw^2*phiw*fwB)-kon*xpt(i).^2*Rt*y(iaB).*(1-y(itheta)-y(ithetaB))+koff*xpt(i)^2*Rt/Abbulk*y(ithetaB)+kspent*y(ia)*xpt(i)^2; % spent interstitium Ab
    dydt(itheta)=kon*Abbulk*y(ia).*(1-y(itheta)-y(ithetaB))-koff*y(itheta)-kr*y(itheta)-kspent*y(itheta); % unspent bound Ab
    dydt(ithetaB)=kon*Abbulk*y(iaB).*(1-y(itheta)-y(ithetaB))-koff*y(ithetaB)-kr*y(ithetaB)+kspent*y(itheta); % spent bound Ab
    dydt(iint)=kr*Rt/Abbulk*y(itheta)-kspent*y(iint); % unspent internalized Ab
    dydt(iintB)=kr*Rt/Abbulk*y(ithetaB)+kspent*y(iint); % spent internalized Ab
end

% Right boundary - spheroid radius
i=np; 
ia=6*(i-1)+1; 
iaB=ia+1;
iam1=ia-6; 
iaBm1=iaB-6; 
itheta=iaB+1; 
ithetaB=itheta+1;
iint=ithetaB+1;
iintB=iint+1;
% West node
xw = (xpt(i)+xpt(i-1))/2; 
phiw=phi(xw); 
% Total flux
fw = (y(ia)-y(iam1))/dx; 
fwB = (y(iaB)-y(iaBm1))/dx; 

% East node
xe = 1;
phie=phi(xe); 
cinfR=cinf*exp(-kspent*t);
cinfB=cinf-cinfR;
if hd<1000
ce=(2*Deff*phie*y(ia)/dx+hd*cinfR)/(2*Deff/dx*phie+hd); 
ceB=(2*Deff*phie*y(iaB)/dx+hd*cinfB)/(2*Deff/dx*phie+hd); 
else
ce=cinfR;
ceB=cinfB;
end
fe = (ce-y(ia))/(dx/2);
feB = (ceB-y(iaB))/(dx/2);

dydt(ia)=Deff/dx*(xe^2*phie*fe-xw^2*phiw*fw)-kon*xpt(i).^2*Rt*y(ia).*(1-y(itheta)-y(ithetaB))+koff*xpt(i)^2*Rt/Abbulk*y(itheta)-kspent*y(itheta)*xpt(i)^2; % unspent interstitium Ab
dydt(iaB)=Deff/dx*(xe^2*phie*feB-xw^2*phiw*fwB)-kon*xpt(i).^2*Rt*y(iaB).*(1-y(itheta)-y(ithetaB))+koff*xpt(i)^2*Rt/Abbulk*y(ithetaB)+kspent*y(itheta)*xpt(i)^2; % spent interstitium Ab
dydt(itheta)=kon*Abbulk*y(ia).*(1-y(itheta)-y(ithetaB))-koff*y(itheta)-kr*y(itheta)-kspent*y(itheta); % unspent bound Ab
dydt(ithetaB)=kon*Abbulk*y(iaB).*(1-y(itheta)-y(ithetaB))-koff*y(ithetaB)-kr*y(ithetaB)+kspent*y(itheta); % spent bound Ab
dydt(iint)=kr*Rt/Abbulk*y(itheta)-kspent*y(iint); % unspent internalized Ab
dydt(iintB)=kr*Rt/Abbulk*y(ithetaB)+kspent*y(iint); % spent internalized Ab