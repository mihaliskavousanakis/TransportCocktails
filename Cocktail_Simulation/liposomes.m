function dydt=liposomes(t,y,xpt,Deff,hd,xpor,Dd,hdd,phiD,cinf)
% kspent
thalf=10*24; %hrs
kspent = log(2)/thalf;

np=length(xpt);
npt=length(y); 
dydt=zeros(npt,1); 
dx=xpt(2)-xpt(1); 

phiL=@(r) xpor(1)*r.^xpor(2)+xpor(3); 
% porosity for drug
phiD=@(r) 0.7313;

% pH spatial distribution in the spheroid
pH =@(x) (7.4929-6.5778)/2*tanh((x-0.6240)*7.8671)+(7.4929+6.5778)/2;

% functional for drug release and drug uptake rate constant
kd=@(x) -3.8195*pH(x) + 28.8336;
kup=@(x) 0*(0.0052*pH(x) -0.0307);


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
fe=(y(ndp)-y(nd))/dx;
dydt(nd)=Deff/dx*(xe^2*phie*fe-xw^2*phiw*fw)-kd(xpt(i))/0.1.*phiL(xpt(i))*xpt(i)^2*y(nd)*content(met)-kspent*y(nd)*xpt(i)^2.*phiL(xpt(i));
% Spent
metB=met+1;
nd=j(metB,i); 
ndp=j(metB,i+1);
fe=(y(ndp)-y(nd))/dx;
dydt(nd)=Deff/dx*(xe^2*phie*fe-xw^2*phiw*fw)-kd(xpt(i))/0.1.*phiL(xpt(i))*xpt(i)^2*y(nd)*content(metB)+kspent*y(j(met,i))*xpt(i)^2.*phiL(xpt(i));

for met=3:2:19
    nd=j(met,i); 
    ndp=j(met,i+1);     
    fe=(y(ndp)-y(nd))/dx;    
    dydt(nd)=Deff/dx*(xe^2*phie*fe-xw^2*phiw*fw)+kd(xpt(i))/0.1.*phiL(xpt(i))*xpt(i)^2*(y(j(met-2,i))*content(met-2)-y(nd)*content(met))-kspent*y(nd)*xpt(i)^2.*phiL(xpt(i));
    % Spent
    metB=met+1;
    nd=j(metB,i); 
    ndp=j(metB,i+1);     
    fe=(y(ndp)-y(nd))/dx;    
    dydt(nd)=Deff/dx*(xe^2*phie*fe-xw^2*phiw*fw)+kd(xpt(i))/0.1.*phiL(xpt(i))*xpt(i)^2*(y(j(metB-2,i))*content(metB-2)-y(nd)*content(metB))+kspent*y(j(met,i))*xpt(i)^2.*phiL(xpt(i));
end
met=21; 
nd=j(met,i); 
ndp=j(met,i+1); 
fe=(y(ndp)-y(nd))/dx;
dydt(nd)=Deff/dx*(xe^2*phie*fe-xw^2*phiw*fw)+kd(xpt(i))/0.1.*phiL(xpt(i))*xpt(i)^2*(y(j(met-2,i))*content(met-2))-kspent*y(nd)*xpt(i)^2.*phiL(xpt(i));
% Spent
metB=met+1; 
nd=j(metB,i); 
ndp=j(metB,i+1); 
fe=(y(ndp)-y(nd))/dx;
dydt(nd)=Deff/dx*(xe^2*phie*fe-xw^2*phiw*fw)+kd(xpt(i))/0.1.*phiL(xpt(i))*xpt(i)^2*(y(j(metB-2,i))*content(metB-2))+kspent*y(j(met,i))*xpt(i)^2.*phiL(xpt(i));

met=23; 
nd=j(met,i); 
ndp=j(met,i+1); 
fe=(y(ndp)-y(nd))/dx; 
fw=0; 
phie=phiD(xe); phiw=phiD(xw);
dydt(nd)=Dd/dx*(xe^2*phie*fe-xw^2*phiw*fw);
for count=1:2:19
dydt(nd)=dydt(nd)+kd(xpt(i)).*phiD(xpt(i)).*xpt(i)^2*y(j(count,i))*content(count);
end
dydt(nd)=dydt(nd)-kup(xpt(i)).*xpt(i).^2*phiD(xpt(i)).*y(nd)-kspent*y(nd)*xpt(i)^2.*phiD(xpt(i)); 
% Spent
metB=met+1;
nd=j(metB,i); 
ndp=j(metB,i+1); 
fe=(y(ndp)-y(nd))/dx; 
fw=0; 
phie=phiD(xe); phiw=phiD(xw);
dydt(nd)=Dd/dx*(xe^2*phie*fe-xw^2*phiw*fw);
for count=2:2:20
dydt(nd)=dydt(nd)+kd(xpt(i)).*phiD(xpt(i)).*xpt(i)^2*y(j(count,i))*content(count);
end
dydt(nd)=dydt(nd)-kup(xpt(i)).*xpt(i).^2*phiD(xpt(i)).*y(nd)+kspent*y(j(met,i))*xpt(i)^2.*phiD(xpt(i)); 

met=25;
nd=j(met,i);
dydt(nd)=kup(xpt(i)).*phiD(xpt(i)).*y(j(met-2,i))-kspent*y(nd);
% Spent
metB=met+1;
nd=j(metB,i);
dydt(nd)=kup(xpt(i)).*phiD(xpt(i)).*y(j(metB-2,i))+kspent*y(j(met,i));

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
fe=(y(ndp)-y(nd))/dx;
fw=(y(nd)-y(ndm))/dx;
dydt(nd)=Deff/dx*(xe.^2.*phie.*fe-xw.^2.*phiw.*fw)-kd(xpt(i))/0.1.*phiL(xpt(i)).*xpt(i).^2.*y(nd)*content(met)-kspent*y(nd).*xpt(i).^2.*phiL(xpt(i));
% Spent
metB=met+1;
nd=j(metB,i)'; 
ndp=j(metB,i+1)'; 
ndm=j(metB,i-1)'; 
fe=(y(ndp)-y(nd))/dx;
fw=(y(nd)-y(ndm))/dx;
dydt(nd)=Deff/dx*(xe.^2.*phie.*fe-xw.^2.*phiw.*fw)-kd(xpt(i))/0.1.*phiL(xpt(i)).*xpt(i).^2.*y(nd)*content(metB)+kspent*y(j(met,i)').*xpt(i).^2.*phiL(xpt(i));

for met=3:2:19
    nd=j(met,i)'; 
    ndp=j(met,i+1)';
    ndm=j(met,i-1)';      
    fe=(y(ndp)-y(nd))/dx;
    fw=(y(nd)-y(ndm))/dx;    
    dydt(nd)=Deff/dx*(xe.^2.*phie.*fe-xw.^2.*phiw.*fw)+kd(xpt(i))/0.1.*phiL(xpt(i)).*xpt(i).^2.*(y(j(met-2,i))*content(met-2)-y(nd)*content(met))-kspent*y(nd).*xpt(i).^2.*phiL(xpt(i));
    % Spent
    metB=met+1;
    nd=j(metB,i)'; 
    ndp=j(metB,i+1)';
    ndm=j(metB,i-1)';      
    fe=(y(ndp)-y(nd))/dx;
    fw=(y(nd)-y(ndm))/dx;    
    dydt(nd)=Deff/dx*(xe.^2.*phie.*fe-xw.^2.*phiw.*fw)+kd(xpt(i))/0.1.*phiL(xpt(i)).*xpt(i).^2.*(y(j(metB-2,i))*content(metB-2)-y(nd)*content(metB))+kspent*y(j(met,i)').*xpt(i).^2.*phiL(xpt(i));
end
met=21; 
nd=j(met,i)'; 
ndp=j(met,i+1)'; 
ndm=j(met,i-1)'; 
fe=(y(ndp)-y(nd))/dx;
fw=(y(nd)-y(ndm))/dx;
dydt(nd)=Deff/dx*(xe.^2.*phie.*fe-xw.^2.*phiw.*fw)+kd(xpt(i))/0.1.*phiL(xpt(i)).*xpt(i).^2.*(y(j(met-2,i))*content(met-2))-kspent*y(nd).*xpt(i).^2.*phiL(xpt(i));
% Spent
metB=met+1;
nd=j(metB,i)'; 
ndp=j(metB,i+1)'; 
ndm=j(metB,i-1)'; 
fe=(y(ndp)-y(nd))/dx;
fw=(y(nd)-y(ndm))/dx;
dydt(nd)=Deff/dx*(xe.^2.*phie.*fe-xw.^2.*phiw.*fw)+kd(xpt(i))/0.1.*phiL(xpt(i)).*xpt(i).^2.*(y(j(metB-2,i))*content(metB-2))+kspent*y(j(met,i)').*xpt(i).^2.*phiL(xpt(i));

met=23; 
nd=j(met,i)'; 
ndp=j(met,i+1)'; 
ndm=j(met,i-1)'; 
fe=(y(ndp)-y(nd))/dx; 
fw=(y(nd)-y(ndm))/dx; 
phie=phiD(xe); phiw=phiD(xw);
dydt(nd)=Dd/dx*(xe.^2.*phie.*fe-xw.^2.*phiw.*fw);
for count=1:2:19
dydt(nd)=dydt(nd)+kd(xpt(i)).*phiD(xpt(i)).*xpt(i).^2.*y(j(count,i))*content(count);
end
dydt(nd)=dydt(nd)-kup(xpt(i)).*xpt(i).^2*phiD(xpt(i)).*y(nd)-kspent*y(nd).*xpt(i).^2.*phiD(xpt(i)); 
% Spent
metB=met+1; 
nd=j(metB,i)'; 
ndp=j(metB,i+1)'; 
ndm=j(metB,i-1)'; 
fe=(y(ndp)-y(nd))/dx; 
fw=(y(nd)-y(ndm))/dx; 
phie=phiD(xe); phiw=phiD(xw);
dydt(nd)=Dd/dx*(xe.^2.*phie.*fe-xw.^2.*phiw.*fw);
for count=2:2:20
dydt(nd)=dydt(nd)+kd(xpt(i)).*phiD(xpt(i)).*xpt(i).^2.*y(j(count,i))*content(count);
end
dydt(nd)=dydt(nd)-kup(xpt(i)).*xpt(i).^2*phiD(xpt(i)).*y(nd)+kspent*y(j(met,i)').*xpt(i).^2.*phiD(xpt(i)); 

met=25;
nd=j(met,i);
dydt(nd)=kup(xpt(i)).*phiD(xpt(i)).*y(j(met-2,i))-kspent*y(nd);
% Spent
metB=met+1;
nd=j(metB,i);
dydt(nd)=kup(xpt(i)).*phiD(xpt(i)).*y(j(metB-2,i))+kspent*y(j(met,i));

% Right boundary - spheroid radius
i=np; 
% West node
xw = (xpt(i)+xpt(i-1))/2; 
phiw=phiL(xw);
%fw = (Ltot(i)-Ltot(i-1))/dx; 
% East node
xe = 1;
phie=phiL(xe);

cinfR=cinf*exp(-kspent*t);
cinfB=cinf-cinfR;
met=1;
nd=j(met,i); 
ndm=j(met,i-1); 
Le = (hd*cinfR+2*Deff*phie/dx*y(nd))/(2*Deff*phie/dx+hd);
fe = (Le-y(nd))/(dx/2);
fw = (y(nd)-y(ndm))/dx;
dydt(nd)=Deff/dx*(xe.^2.*phie.*fe-xw.^2.*phiw.*fw)-kd(xpt(i))/0.1.*phiL(xpt(i)).*xpt(i).^2*y(nd)-kspent*y(nd).*xpt(i).^2.*phiL(xpt(i));
% Spent
metB=met+1;
nd=j(metB,i); 
ndm=j(metB,i-1); 
Le = (hd*cinfB+2*Deff*phie/dx*y(nd))/(2*Deff*phie/dx+hd);
fe = (Le-y(nd))/(dx/2);
fw = (y(nd)-y(ndm))/dx;
dydt(nd)=Deff/dx*(xe.^2.*phie.*fe-xw.^2.*phiw.*fw)-kd(xpt(i))/0.1.*phiL(xpt(i)).*xpt(i).^2*y(nd)+kspent*y(j(met,i)).*xpt(i).^2.*phiL(xpt(i));
for met=3:2:19
    nd=j(met,i); 
    ndm=j(met,i-1);     
    Le = (hd*0+2*Deff*phie/dx*y(nd))/(2*Deff*phie/dx+hd);
    fe = (Le-y(nd))/(dx/2);
    fw = (y(nd)-y(ndm))/dx;    
    dydt(nd)=Deff/dx*(xe.^2*phie.*fe-xw.^2.*phiw.*fw)+kd(xpt(i))/0.1.*phiL(xpt(i)).*xpt(i).^2.*(y(j(met-2,i))*content(met-2)-y(nd)*content(met))-kspent*y(nd).*xpt(i).^2.*phiL(xpt(i));
    % Spent
    metB=met+1;
    nd=j(metB,i); 
    ndm=j(metB,i-1);     
    Le = (hd*0+2*Deff*phie/dx*y(nd))/(2*Deff*phie/dx+hd);
    fe = (Le-y(nd))/(dx/2);
    fw = (y(nd)-y(ndm))/dx;    
    dydt(nd)=Deff/dx*(xe.^2*phie.*fe-xw.^2.*phiw.*fw)+kd(xpt(i))/0.1.*phiL(xpt(i)).*xpt(i).^2.*(y(j(metB-2,i))*content(metB-2)-y(nd)*content(metB))+kspent*y(j(met,i)).*xpt(i).^2.*phiL(xpt(i));
end
met=21; 
nd=j(met,i); 
ndm=j(met,i-1); 
Le = (hd*0+2*Deff*phie/dx*y(nd))/(2*Deff*phie/dx+hd);
fe = (Le-y(nd))/(dx/2);
fw = (y(nd)-y(ndm))/dx;
dydt(nd)=Deff/dx*(xe.^2.*phie.*fe-xw.^2.*phiw.*fw)+kd(xpt(i))/0.1.*phiL(xpt(i)).*xpt(i).^2.*(y(j(met-2,i))*content(met-2))-kspent*y(nd).*xpt(i).^2.*phiL(xpt(i));
% Spent
metB=met+1; 
nd=j(metB,i); 
ndm=j(metB,i-1); 
Le = (hd*0+2*Deff*phie/dx*y(nd))/(2*Deff*phie/dx+hd);
fe = (Le-y(nd))/(dx/2);
fw = (y(nd)-y(ndm))/dx;
dydt(nd)=Deff/dx*(xe.^2.*phie.*fe-xw.^2.*phiw.*fw)+kd(xpt(i))/0.1.*phiL(xpt(i)).*xpt(i).^2.*(y(j(metB-2,i))*content(metB-2))+kspent*y(j(met,i)).*xpt(i).^2.*phiL(xpt(i));

met=23; 
nd=j(met,i); 
ndm=j(met,i-1); 

ce=(hd*0+2*Dd*phiD(xe)/dx*y(nd))/(2*Dd*phiD(xe)/dx+hd);
fe=(ce-y(nd))/(dx/2); 
fw=(y(nd)-y(ndm))/dx; 
phie=phiD(xe); phiw=phiD(xw);
dydt(nd)=Dd/dx*(xe.^2.*phie.*fe-xw.^2.*phiw.*fw);
for count=1:2:19
dydt(nd)=dydt(nd)+kd(xpt(i)).*phiD(xpt(i)).*xpt(i).^2.*(y(j(count,i))*content(count));
end
dydt(nd)=dydt(nd)-kup(xpt(i)).*xpt(i).^2*phiD(xpt(i)).*y(nd)-kspent*y(nd).*xpt(i).^2.*phiD(xpt(i)); 
% Spent
metB=met+1; 
nd=j(metB,i); 
ndm=j(metB,i-1); 

ce=(hd*0+2*Dd*phiD(xe)/dx*y(nd))/(2*Dd*phiD(xe)/dx+hd);
fe=(ce-y(nd))/(dx/2); 
fw=(y(nd)-y(ndm))/dx; 
phie=phiD(xe); phiw=phiD(xw);
dydt(nd)=Dd/dx*(xe.^2.*phie.*fe-xw.^2.*phiw.*fw);
for count=2:2:20
dydt(nd)=dydt(nd)+kd(xpt(i)).*phiD(xpt(i)).*xpt(i).^2.*(y(j(count,i))*content(count));
end
dydt(nd)=dydt(nd)-kup(xpt(i)).*xpt(i).^2*phiD(xpt(i)).*y(nd)+kspent*y(j(met,i)).*xpt(i).^2.*phiD(xpt(i)); 


met=25;
nd=j(met,i);
dydt(nd)=kup(xpt(i)).*phiD(xpt(i)).*y(j(met-2,i))-kspent*y(nd);
% Spent
metB=met+1;
nd=j(metB,i);
dydt(nd)=kup(xpt(i)).*phiD(xpt(i)).*y(j(metB-2,i))+kspent*y(j(met,i));

