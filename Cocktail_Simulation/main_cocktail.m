clear all
% Create liposomes distribution 
liposomes_uptake_clearance_simulation

cnt=1:-0.1:0; % content percentage
for i=1:length(cnt)
    content(2*(i-1)+(1:2))=cnt(i); % this accounts for spent/unspent species
end
Dox_retained=zeros(1,np);
DoxB_retained=zeros(1,np);
i=1:np; 

Drug_Lip = 1*0.06*0.50; %uM
Lipsol = Drug_Lip/0.002; %uM
Doxfree_intern_Lip=[];Dox_retained_Lip=[];
DoxBfree_intern_Lip=[];DoxB_retained_Lip=[];
for count=1:12
    Dox_retained=zeros(1,np); 
    DoxB_retained=zeros(1,np);
    for met=1:2:21
        Dox_retained=Dox_retained+Ytheorup(count,26*(i-1)+met)*content(met)*Drug_Lip*1000; %nM
    end
    
    for met=2:2:22
        DoxB_retained=DoxB_retained+Ytheorup(count,26*(i-1)+met)*content(met)*Drug_Lip*1000; %nM
    end
    

Doxfree_intern = (Ytheorup(count,26*(i-1)+23)+Ytheorup(count,26*(i-1)+25))*Drug_Lip*1000; %nM
DoxBfree_intern = (Ytheorup(count,26*(i-1)+24)+Ytheorup(count,26*(i-1)+26))*Drug_Lip*1000; %nM

Doxfree_intern_Lip(end+1,:)=Doxfree_intern;
Dox_retained_Lip(end+1,:)=Dox_retained;

DoxBfree_intern_Lip(end+1,:)=DoxBfree_intern;
DoxB_retained_Lip(end+1,:)=DoxB_retained;
end

for count=1:193
    Dox_retained=zeros(1,np); 
    for met=1:2:21
        Dox_retained=Dox_retained+Ytheorcl(count,26*(i-1)+met)*content(met)*Drug_Lip*1000; %nM
    end
    

    DoxB_retained=zeros(1,np); 
    for met=2:2:22
        DoxB_retained=DoxB_retained+Ytheorcl(count,26*(i-1)+met)*content(met)*Drug_Lip*1000; %nM
    end
    

Doxfree_intern = (Ytheorcl(count,26*(i-1)+23)+Ytheorcl(count,26*(i-1)+25))*Drug_Lip*1000; %nM
Doxfree_intern_Lip(end+1,:)=Doxfree_intern;
Dox_retained_Lip(end+1,:)=Dox_retained;

DoxBfree_intern = (Ytheorcl(count,26*(i-1)+24)+Ytheorcl(count,26*(i-1)+26))*Drug_Lip*1000; %nM
DoxBfree_intern_Lip(end+1,:)=DoxBfree_intern;
DoxB_retained_Lip(end+1,:)=DoxB_retained;
end


save liposomes_drug 'Dox_retained_Lip' 'Doxfree_intern_Lip' 'DoxB_retained_Lip' 'DoxBfree_intern_Lip'
clear



Drug_Ab = 0.06*0.50; %uM (if 1 mole of drug corresponds to 4 moles of Ab)
load antibodies_experimental.mat
Abbulk=Drug_Ab*4*1000; %nM
antibodies_uptake_clearance_simulation
j=1:np; 

Dox_Antibodies=[]; DoxB_Antibodies=[];
for count=1:48
y1=Ytheorup(count,6*(j-1)+1);
y1B=Ytheorup(count,6*(j-1)+2);
y2=Ytheorup(count,6*(j-1)+3)*Rt/Abbulk;
y2B=Ytheorup(count,6*(j-1)+4)*Rt/Abbulk;
y3=Ytheorup(count,6*(j-1)+5);
y3B=Ytheorup(count,6*(j-1)+6);

Dox_Ab = (y1+y2+y3)*Abbulk/4; %nM
DoxB_Ab=(y1B+y2B+y3B)*Abbulk/4; %nM

Dox_Antibodies(end+1,:)=Dox_Ab; 
DoxB_Antibodies(end+1,:)=DoxB_Ab; 
end

for count=1:193
y1=1*Ytheorcl(count,6*(j-1)+1);
y1B=1*Ytheorcl(count,6*(j-1)+2);
y2=1*Ytheorcl(count,6*(j-1)+3)*Rt/Abbulk;
y2B=1*Ytheorcl(count,6*(j-1)+4)*Rt/Abbulk;
y3=Ytheorcl(count,6*(j-1)+5);
y3B=Ytheorcl(count,6*(j-1)+6);

Dox_Ab = (y1+y2+y3)*Abbulk/4; %nM
DoxB_Ab = (y1B+y2B+y3B)*Abbulk/4; %nM
Dox_Antibodies(end+1,:)=Dox_Ab; 
DoxB_Antibodies(end+1,:)=DoxB_Ab; 
end


load liposomes_drug
save drug_distribution_Ab_100_Lip_0

Drug_critical=5;
% For chemo
Fr_Chemo=[];
Drug_Chemo=Dox_Antibodies;
Drug_Chemo(37:end,:)=Drug_Chemo(37:end,:)+Doxfree_intern_Lip;
sz=size(Drug_Chemo); 
for count=1:sz(1)
    [i,j]=find(Drug_Chemo(count,:)>=Drug_critical);
    if isempty(j)
        Fr_Chemo(end+1)=0;
    else
        Fr_Chemo(end+1)=100*((xpt(j(end))/xpt(end))^3-(xpt(j(1))/xpt(end))^3);
    end
end

% For radio
% For radio
Fr_Radio=[];
Drug_Radio=Dox_Antibodies;
DrugB_Radio=DoxB_Antibodies;
Drug_Radio(37:end,:)=Drug_Radio(37:end,:)+Doxfree_intern_Lip+Dox_retained_Lip;
DrugB_Radio(37:end,:)=DrugB_Radio(37:end,:)+DoxBfree_intern_Lip+DoxB_retained_Lip;


Ti=0:0.5:5*24;

Int_Radio=zeros(size(Drug_Radio));
for i=1:100
for j=2:241
Int_Radio(j,i)=trapz(Ti(1:j),Drug_Radio(1:j,i));
end
end

Int2_Radio=zeros(size(Int_Radio));
for i=1:100
for j=2:241
Int2_Radio(j,i)=trapz(Ti(1:j),Int_Radio(1:j,i));
end
end

save results_radio 'Int_Radio' 'Ti' 'Int2_Radio' 'xpt'


    


