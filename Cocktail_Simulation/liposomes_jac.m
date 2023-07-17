%function ajac=liposomes_jac(t,y,xpt,Deff,hd,kd,xpor,Dd,phiD,cinf)
function ajac=liposomes_jac(t,y,xpt,Deff,hd,xpor,Dd,phiD,cinf)
% kspent
thalf=10*24; %hrs
kspent = log(2)/thalf;

np=length(xpt);
npt=length(y); 
ajac=sparse(npt,npt); 
dx=xpt(2)-xpt(1); 

phiL=@(r) xpor(1)*r.^xpor(2)+xpor(3); 
%load porosity_drug_ne
%phiD=@(r) xporD(1)*r.^xporD(2)+xporD(3); 
phiD=@(r) 0.7313;
pH =@(x) (7.4929-6.5778)/2*tanh((x-0.6240)*7.8671)+(7.4929+6.5778)/2;
%pH = @(x) 7.4;
kd=@(x) -3.8195*pH(x) + 28.8336;
kup=@(x) 0.0052*pH(x) -0.0307;

i=1:np; i=i';
for met=1:26
    j(met,:)=26*(i-1)+met; 
end


cnt=1:-0.1:0; % content percentage
for i=1:length(cnt)
    content(2*(i-1)+(1:2))=cnt(i); % this accounts for spent/unspent species
end

% Symmetry boundary 
i=1;
% West node
xw=0; 
fw=0; 
phiw=phiL(xw);

% East node
xe=(xpt(i)+xpt(i+1))/2;
phie=phiL(xe);

met=1;
nd=j(met,i); 
ndp=j(met,i+1);

ajac=ajac+sparse(nd,nd,-Deff/dx^2*xe^2*phie-kd(xpt(i))/0.1.*phiL(xpt(i))*xpt(i)^2*content(met)-kspent*xpt(i)^2.*phiL(xpt(i)),npt,npt);
ajac=ajac+sparse(nd,ndp,Deff/dx^2*xe^2*phie,npt,npt);
% Spent
metB=met+1;
nd=j(metB,i); 
ndp=j(metB,i+1);
ajac=ajac+sparse(nd,nd,-Deff/dx^2*xe^2*phie-kd(xpt(i))/0.1.*phiL(xpt(i))*xpt(i)^2*content(metB),npt,npt);
ajac=ajac+sparse(nd,ndp,Deff/dx^2*xe^2*phie,npt,npt);
ajac=ajac+sparse(nd,j(met,i),kspent*xpt(i)^2.*phiL(xpt(i)),npt,npt);

for met=3:2:19
    nd=j(met,i); 
    ndp=j(met,i+1);        
    ajac=ajac+sparse(nd,nd,-Deff/dx^2*xe^2*phie-kd(xpt(i))/0.1.*phiL(xpt(i))*xpt(i)^2*content(met)-kspent*xpt(i)^2.*phiL(xpt(i)),npt,npt);
    ajac=ajac+sparse(nd,ndp,Deff/dx^2*xe^2*phie,npt,npt);
    ajac=ajac+sparse(nd,j(met-2,i),kd(xpt(i))/0.1*phiL(xpt(i)).*xpt(i)^2*content(met-2),npt,npt);
    % Spent
    metB=met+1;
    nd=j(metB,i); 
    ndp=j(metB,i+1);        
    ajac=ajac+sparse(nd,nd,-Deff/dx^2*xe^2*phie-kd(xpt(i))/0.1.*phiL(xpt(i))*xpt(i)^2*content(metB),npt,npt);
    ajac=ajac+sparse(nd,ndp,Deff/dx^2*xe^2*phie,npt,npt);
    ajac=ajac+sparse(nd,j(metB-2,i),kd(xpt(i))/0.1*phiL(xpt(i)).*xpt(i)^2*content(metB-2),npt,npt);
    ajac=ajac+sparse(nd,j(met,i),kspent*xpt(i)^2.*phiL(xpt(i)),npt,npt);
end
met=21; 
nd=j(met,i); 
ndp=j(met,i+1); 
ajac=ajac+sparse(nd,nd,-Deff/dx^2*xe^2*phie-kspent*xpt(i)^2.*phiL(xpt(i)),npt,npt);
ajac=ajac+sparse(nd,ndp,Deff/dx^2*xe^2*phie,npt,npt);
ajac=ajac+sparse(nd,j(met-2,i),kd(xpt(i))/0.1.*phiL(xpt(i))*xpt(i)^2*content(met-2),npt,npt);
% Spent
metB=met+1;
nd=j(metB,i); 
ndp=j(metB,i+1); 
ajac=ajac+sparse(nd,nd,-Deff/dx^2*xe^2*phie,npt,npt);
ajac=ajac+sparse(nd,ndp,Deff/dx^2*xe^2*phie,npt,npt);
ajac=ajac+sparse(nd,j(metB-2,i),kd(xpt(i))/0.1.*phiL(xpt(i))*xpt(i)^2*content(metB-2),npt,npt);
ajac=ajac+sparse(nd,j(met,i),kspent*xpt(i)^2.*phiL(xpt(i)),npt,npt);

met=23; 
nd=j(met,i); 
ndp=j(met,i+1); 
phie=phiD(xe); phiw=phiD(xw);
ajac=ajac+sparse(nd,nd,-Dd/dx^2*xe^2*phie,npt,npt);
ajac=ajac+sparse(nd,ndp,Dd/dx^2*xe^2*phie,npt,npt);
for count=1:2:19    
    ajac=ajac+sparse(nd,j(count,i),kd(xpt(i)).*phiD(xpt(i))*xpt(i)^2*content(count),npt,npt);
end
ajac=ajac+sparse(nd,nd,-kup(xpt(i)).*xpt(i).^2*phiD(xpt(i))-kspent*xpt(i)^2.*phiD(xpt(i)),npt,npt);
% Spent
metB=met+1; 
nd=j(metB,i); 
ndp=j(metB,i+1); 
phie=phiD(xe); phiw=phiD(xw);
ajac=ajac+sparse(nd,nd,-Dd/dx^2*xe^2*phie,npt,npt);
ajac=ajac+sparse(nd,ndp,Dd/dx^2*xe^2*phie,npt,npt);
for count=2:2:20    
    ajac=ajac+sparse(nd,j(count,i),kd(xpt(i)).*phiD(xpt(i))*xpt(i)^2*content(count),npt,npt);
end
ajac=ajac+sparse(nd,nd,-kup(xpt(i)).*xpt(i).^2*phiD(xpt(i)),npt,npt);
ajac=ajac+sparse(nd,j(met,i),kspent*xpt(i)^2.*phiD(xpt(i)),npt,npt);

met=25;
nd=j(met,i);
ajac=ajac+sparse(nd,j(met-2,i),kup(xpt(i)).*phiD(xpt(i)),npt,npt);
ajac=ajac+sparse(nd,nd,-kspent,npt,npt);
% Spent
metB=met+1;
nd=j(metB,i);
ajac=ajac+sparse(nd,j(metB-2,i),kup(xpt(i)).*phiD(xpt(i)),npt,npt);
ajac=ajac+sparse(nd,j(met,i),kspent,npt,npt);


% inner nodes 
i=2:np-1; i=i';
% West node
xw = (xpt(i)+xpt(i-1))/2; 
phiw=phiL(xw);
% East node
xe=(xpt(i)+xpt(i+1))/2;
phie=phiL(xe);
     
met=1;
nd=j(met,i)'; 
ndp=j(met,i+1)'; 
ndm=j(met,i-1)'; 
ajac=ajac+sparse(nd,nd,-Deff/dx^2*xe.^2.*phie-Deff/dx^2*xw.^2.*phiw-kd(xpt(i))/0.1.*phiL(xpt(i)).*xpt(i).^2*content(met)-kspent*phiL(xpt(i)).*xpt(i).^2,npt,npt);
ajac=ajac+sparse(nd,ndm,Deff/dx^2*xw.^2.*phiw,npt,npt);
ajac=ajac+sparse(nd,ndp,Deff/dx^2*xe.^2.*phie,npt,npt);
% Spent
metB=met+1;
nd=j(metB,i)'; 
ndp=j(metB,i+1)'; 
ndm=j(metB,i-1)'; 
ajac=ajac+sparse(nd,nd,-Deff/dx^2*xe.^2.*phie-Deff/dx^2*xw.^2.*phiw-kd(xpt(i))/0.1.*phiL(xpt(i)).*xpt(i).^2*content(metB),npt,npt);
ajac=ajac+sparse(nd,ndm,Deff/dx^2*xw.^2.*phiw,npt,npt);
ajac=ajac+sparse(nd,ndp,Deff/dx^2*xe.^2.*phie,npt,npt);
ajac=ajac+sparse(nd,j(met,i)',kspent*phiL(xpt(i)).*xpt(i).^2,npt,npt);
for met=3:2:19
    nd=j(met,i)'; 
    ndp=j(met,i+1)';
    ndm=j(met,i-1)';     
    ajac=ajac+sparse(nd,nd,-Deff/dx^2*xe.^2.*phie-Deff/dx^2*xw.^2.*phiw-kd(xpt(i))/0.1.*phiL(xpt(i)).*xpt(i).^2*content(met)-kspent*phiL(xpt(i)).*xpt(i).^2,npt,npt);
    ajac=ajac+sparse(nd,ndm,Deff/dx^2*xw.^2.*phiw,npt,npt);
    ajac=ajac+sparse(nd,ndp,Deff/dx^2*xe.^2.*phie,npt,npt);
    ajac=ajac+sparse(nd,j(met-2,i),kd(xpt(i))/0.1.*phiL(xpt(i)).*xpt(i).^2*content(met-2),npt,npt);
    % Spent
    metB=met+1;
    nd=j(metB,i)'; 
    ndp=j(metB,i+1)';
    ndm=j(metB,i-1)';     
    ajac=ajac+sparse(nd,nd,-Deff/dx^2*xe.^2.*phie-Deff/dx^2*xw.^2.*phiw-kd(xpt(i))/0.1.*phiL(xpt(i)).*xpt(i).^2*content(metB),npt,npt);
    ajac=ajac+sparse(nd,ndm,Deff/dx^2*xw.^2.*phiw,npt,npt);
    ajac=ajac+sparse(nd,ndp,Deff/dx^2*xe.^2.*phie,npt,npt);
    ajac=ajac+sparse(nd,j(metB-2,i),kd(xpt(i))/0.1.*phiL(xpt(i)).*xpt(i).^2*content(metB-2),npt,npt);
    ajac=ajac+sparse(nd,j(met,i)',kspent*phiL(xpt(i)).*xpt(i).^2,npt,npt);
end
met=21; 
nd=j(met,i)'; 
ndp=j(met,i+1)'; 
ndm=j(met,i-1)'; 
ajac=ajac+sparse(nd,nd,-Deff/dx^2*xe.^2.*phie-Deff/dx^2*xw.^2.*phiw-kspent*phiL(xpt(i)).*xpt(i).^2,npt,npt);
ajac=ajac+sparse(nd,ndm,Deff/dx^2*xw.^2.*phiw,npt,npt);
ajac=ajac+sparse(nd,ndp,Deff/dx^2*xe.^2.*phie,npt,npt);
ajac=ajac+sparse(nd,j(met-2,i),kd(xpt(i))/0.1.*phiL(xpt(i)).*xpt(i).^2*content(met-2),npt,npt);
% Spent
metB=met+1; 
nd=j(metB,i)'; 
ndp=j(metB,i+1)'; 
ndm=j(metB,i-1)'; 
ajac=ajac+sparse(nd,nd,-Deff/dx^2*xe.^2.*phie-Deff/dx^2*xw.^2.*phiw,npt,npt);
ajac=ajac+sparse(nd,ndm,Deff/dx^2*xw.^2.*phiw,npt,npt);
ajac=ajac+sparse(nd,ndp,Deff/dx^2*xe.^2.*phie,npt,npt);
ajac=ajac+sparse(nd,j(metB-2,i),kd(xpt(i))/0.1.*phiL(xpt(i)).*xpt(i).^2*content(metB-2),npt,npt);
ajac=ajac+sparse(nd,j(met,i)',kspent*phiL(xpt(i)).*xpt(i).^2,npt,npt);

met=23; 
nd=j(met,i)'; 
ndp=j(met,i+1)'; 
ndm=j(met,i-1)'; 
phie=phiD(xe); phiw=phiD(xw);
ajac=ajac+sparse(nd,nd,-Dd/dx^2*xe.^2.*phie-Dd/dx^2*xw.^2.*phiw,npt,npt);
ajac=ajac+sparse(nd,ndm,Dd/dx^2*xw.^2.*phiw,npt,npt);
ajac=ajac+sparse(nd,ndp,Dd/dx^2*xe.^2.*phie,npt,npt);
for count=1:2:19
    ajac=ajac+sparse(nd,j(count,i),kd(xpt(i)).*phiD(xpt(i)).*xpt(i).^2.*content(count),npt,npt);
end
ajac=ajac+sparse(nd,nd,-kup(xpt(i)).*xpt(i).^2*phiD(xpt(i))-kspent*xpt(i).^2.*phiD(xpt(i)),npt,npt);
% Spent
metB=met+1; 
nd=j(metB,i)'; 
ndp=j(metB,i+1)'; 
ndm=j(metB,i-1)'; 
phie=phiD(xe); phiw=phiD(xw);
ajac=ajac+sparse(nd,nd,-Dd/dx^2*xe.^2.*phie-Dd/dx^2*xw.^2.*phiw,npt,npt);
ajac=ajac+sparse(nd,ndm,Dd/dx^2*xw.^2.*phiw,npt,npt);
ajac=ajac+sparse(nd,ndp,Dd/dx^2*xe.^2.*phie,npt,npt);
for count=2:2:20
    ajac=ajac+sparse(nd,j(count,i),kd(xpt(i)).*phiD(xpt(i)).*xpt(i).^2.*content(count),npt,npt);
end
ajac=ajac+sparse(nd,nd,-kup(xpt(i)).*xpt(i).^2*phiD(xpt(i)),npt,npt);
ajac=ajac+sparse(nd,j(met,i)',kspent.*xpt(i).^2*phiD(xpt(i)),npt,npt);

met=25;
nd=j(met,i);
ajac=ajac+sparse(nd,j(met-2,i),kup(xpt(i)).*phiD(xpt(i)),npt,npt);
ajac=ajac+sparse(nd,nd,-kspent,npt,npt);
% Spent
metB=met+1;
nd=j(metB,i);
ajac=ajac+sparse(nd,j(metB-2,i),kup(xpt(i)).*phiD(xpt(i)),npt,npt);
ajac=ajac+sparse(nd,j(met,i),kspent,npt,npt);

% Right boundary - spheroid radius
i=np; 
% West node
xw = (xpt(i)+xpt(i-1))/2; 
phiw=phiL(xw);
% East node
xe = 1;
phie=phiL(xe);


met=1;
nd=j(met,i); 
ndm=j(met,i-1); 
ajac=ajac+sparse(nd,nd,- kd(xpt(i))/0.1*phiL(xpt(i))*xpt(i)^2 - (dx*(2*Deff*hd*phie*xe^2 + Deff*hd*phiw*xw^2) + 2*Deff^2*phie*phiw*xw^2)/(dx^2*(2*Deff*phie + dx*hd))-kspent*xpt(i).^2.*phiL(xpt(i)),npt,npt);
ajac=ajac+sparse(nd,ndm,(Deff*phiw*xw^2)/dx^2,npt,npt);
% Spent
metB=met+1;
nd=j(metB,i); 
ndm=j(metB,i-1); 
ajac=ajac+sparse(nd,nd,- kd(xpt(i))/0.1*phiL(xpt(i))*xpt(i)^2 - (dx*(2*Deff*hd*phie*xe^2 + Deff*hd*phiw*xw^2) + 2*Deff^2*phie*phiw*xw^2)/(dx^2*(2*Deff*phie + dx*hd)),npt,npt);
ajac=ajac+sparse(nd,ndm,(Deff*phiw*xw^2)/dx^2,npt,npt);
ajac=ajac+sparse(nd,j(met,i),kspent*xpt(i).^2.*phiL(xpt(i)),npt,npt);
for met=3:2:19
    nd=j(met,i); 
    ndm=j(met,i-1);      
    ajac=ajac+sparse(nd,nd,- content(met)*kd(xpt(i))/0.1*phiL(xpt(i))*xpt(i)^2 - (dx*(2*Deff*hd*phie*xe^2 + Deff*hd*phiw*xw^2) + 2*Deff^2*phie*phiw*xw^2)/(dx^2*(2*Deff*phie + dx*hd))-kspent*xpt(i).^2.*phiL(xpt(i)),npt,npt);
    ajac=ajac+sparse(nd,ndm,(Deff*phiw*xw^2)/dx^2,npt,npt);
    ajac=ajac+sparse(nd,j(met-2,i),content(met-2)*kd(xpt(i))/0.1*phiL(xpt(i))*xpt(i)^2,npt,npt);
    % Spent
    metB=met+1;
    nd=j(metB,i); 
    ndm=j(metB,i-1);      
    ajac=ajac+sparse(nd,nd,- content(metB)*kd(xpt(i))/0.1*phiL(xpt(i))*xpt(i)^2 - (dx*(2*Deff*hd*phie*xe^2 + Deff*hd*phiw*xw^2) + 2*Deff^2*phie*phiw*xw^2)/(dx^2*(2*Deff*phie + dx*hd)),npt,npt);
    ajac=ajac+sparse(nd,ndm,(Deff*phiw*xw^2)/dx^2,npt,npt);
    ajac=ajac+sparse(nd,j(metB-2,i),content(metB-2)*kd(xpt(i))/0.1*phiL(xpt(i))*xpt(i)^2,npt,npt);
    ajac=ajac+sparse(nd,j(met,i),kspent*xpt(i).^2.*phiL(xpt(i)),npt,npt);
end
met=21; 
nd=j(met,i); 
ndm=j(met,i-1); 
ajac=ajac+sparse(nd,nd,-(Deff*dx*(2*hd*phie*xe^2 + hd*phiw*xw^2) + 2*Deff^2*phie*phiw*xw^2)/(dx^2*(2*Deff*phie + dx*hd))-kspent*xpt(i).^2.*phiL(xpt(i)),npt,npt);
ajac=ajac+sparse(nd,ndm,(Deff*phiw*xw^2)/dx^2,npt,npt);
ajac=ajac+sparse(nd,j(met-2,i),content(met-2)*kd(xpt(i))/0.1*phiL(xpt(i))*xpt(i)^2,npt,npt);
% Spent
metB=met+1; 
nd=j(metB,i); 
ndm=j(metB,i-1); 
ajac=ajac+sparse(nd,nd,-(Deff*dx*(2*hd*phie*xe^2 + hd*phiw*xw^2) + 2*Deff^2*phie*phiw*xw^2)/(dx^2*(2*Deff*phie + dx*hd)),npt,npt);
ajac=ajac+sparse(nd,ndm,(Deff*phiw*xw^2)/dx^2,npt,npt);
ajac=ajac+sparse(nd,j(metB-2,i),content(metB-2)*kd(xpt(i))/0.1*phiL(xpt(i))*xpt(i)^2,npt,npt);
ajac=ajac+sparse(nd,j(met,i),kspent*xpt(i).^2.*phiL(xpt(i)),npt,npt);

met=23; 
nd=j(met,i); 
ndm=j(met,i-1); 
hdd=1*hd;
%if hd>10
 %   hdd=10;
%end
phie=phiD(xe); phiw=phiD(xw);
ajac=ajac+sparse(nd,nd,- (Dd*phiw*xw^2)/dx^2 - (2*Dd*dx*hdd*phie*xe^2)/(hdd*dx^3 + 2*Dd*phie*dx^2),npt,npt);
ajac=ajac+sparse(nd,ndm,(Dd*phiw*xw^2)/dx^2,npt,npt);
for count=1:10    
    ajac=ajac+sparse(nd,j(count,i),kd(xpt(i))*phiD(xpt(i)).*xpt(i).^2.*(content(count)),npt,npt);
end
ajac=ajac+sparse(nd,nd,-kup(xpt(i)).*xpt(i).^2*phiD(xpt(i))-kspent*xpt(i).^2.*phiD(xpt(i)),npt,npt);
% Spent
metB=met+1; 
nd=j(metB,i); 
ndm=j(metB,i-1); 
hdd=1*hd;
%if hd>10
 %   hdd=10;
%end
phie=phiD(xe); phiw=phiD(xw);
ajac=ajac+sparse(nd,nd,- (Dd*phiw*xw^2)/dx^2 - (2*Dd*dx*hdd*phie*xe^2)/(hdd*dx^3 + 2*Dd*phie*dx^2),npt,npt);
ajac=ajac+sparse(nd,ndm,(Dd*phiw*xw^2)/dx^2,npt,npt);
for count=1:10    
    ajac=ajac+sparse(nd,j(count,i),kd(xpt(i))*phiD(xpt(i)).*xpt(i).^2.*(content(count)),npt,npt);
end
ajac=ajac+sparse(nd,nd,-kup(xpt(i)).*xpt(i).^2*phiD(xpt(i)),npt,npt);
ajac=ajac+sparse(nd,j(met,i),kspent*xpt(i).^2.*phiD(xpt(i)),npt,npt);

met=25;
nd=j(met,i);
ajac=ajac+sparse(nd,j(met-2,i),kup(xpt(i)).*phiD(xpt(i)),npt,npt);
ajac=ajac+sparse(nd,nd,-kspent,npt,npt);
% Spent
metB=met+1;
nd=j(metB,i);
ajac=ajac+sparse(nd,j(metB-2,i),kup(xpt(i)).*phiD(xpt(i)),npt,npt);
ajac=ajac+sparse(nd,j(met,i),kspent,npt,npt);
