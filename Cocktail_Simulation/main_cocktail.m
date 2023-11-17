clear all, close all,clc
% Create liposomes distribution 
liposomes_uptake_clearance_simulation

cnt=1:-0.1:0; % content percentage
for i=1:length(cnt)
    content(2*(i-1)+(1:2))=cnt(i); % this accounts for spent/unspent species
end
% unspent/active drug retained in liposomes
Drug_retained=zeros(1,np);
% spent/inactive drug retained in liposomes
DrugB_retained=zeros(1,np);
i=1:np; 

% Fraction of drug carried by liposomes
Lf=0.5;

Drug_Lip = 1*0.06*Lf; %uM
Lipsol = Drug_Lip/0.002; %uM
Drugfree_intern_Lip=[];Drug_retained_Lip=[];
DrugBfree_intern_Lip=[];DrugB_retained_Lip=[];
for count=1:12
    Drug_retained=zeros(1,np); 
    DrugB_retained=zeros(1,np);
    for met=1:2:21
        Drug_retained=Drug_retained+Ytheorup(count,26*(i-1)+met)*content(met)*Drug_Lip*1000; %nM
    end
    
    for met=2:2:22
        DrugB_retained=DrugB_retained+Ytheorup(count,26*(i-1)+met)*content(met)*Drug_Lip*1000; %nM
    end
    

Drugfree_intern = (Ytheorup(count,26*(i-1)+23)+Ytheorup(count,26*(i-1)+25))*Drug_Lip*1000; %nM
DrugBfree_intern = (Ytheorup(count,26*(i-1)+24)+Ytheorup(count,26*(i-1)+26))*Drug_Lip*1000; %nM

Drugfree_intern_Lip(end+1,:)=Drugfree_intern;
Drug_retained_Lip(end+1,:)=Drug_retained;

DrugBfree_intern_Lip(end+1,:)=DrugBfree_intern;
DrugB_retained_Lip(end+1,:)=DrugB_retained;
end

for count=1:17
    Drug_retained=zeros(1,np); 
    for met=1:2:21
        Drug_retained=Drug_retained+Ytheorcl(count,26*(i-1)+met)*content(met)*Drug_Lip*1000; %nM
    end
    

    DrugB_retained=zeros(1,np); 
    for met=2:2:22
        DrugB_retained=DrugB_retained+Ytheorcl(count,26*(i-1)+met)*content(met)*Drug_Lip*1000; %nM
    end
    

Drugfree_intern = (Ytheorcl(count,26*(i-1)+23)+Ytheorcl(count,26*(i-1)+25))*Drug_Lip*1000; %nM
Drugfree_intern_Lip(end+1,:)=Drugfree_intern;
Drug_retained_Lip(end+1,:)=Drug_retained;

DrugBfree_intern = (Ytheorcl(count,26*(i-1)+24)+Ytheorcl(count,26*(i-1)+26))*Drug_Lip*1000; %nM
DrugBfree_intern_Lip(end+1,:)=DrugBfree_intern;
DrugB_retained_Lip(end+1,:)=DrugB_retained;
end


save liposomes_drug 'Drug_retained_Lip' 'Drugfree_intern_Lip' 'DrugB_retained_Lip' 'DrugBfree_intern_Lip'
clearvars -except Lf xpt



Drug_Ab = 0.06*(1-Lf); %uM (if 1 mole of drug corresponds to 4 moles of Ab)

if Drug_Ab>0
    load antibodies_experimental.mat
    Abbulk=Drug_Ab*4*1000; %nM

    antibodies_uptake_clearance_simulation
    j=1:np; 

    Drug_Antibodies=[]; DrugB_Antibodies=[];
    for count=1:48
        y1=Ytheorup(count,6*(j-1)+1);
        y1B=Ytheorup(count,6*(j-1)+2);
        y2=Ytheorup(count,6*(j-1)+3)*Rt/Abbulk;
        y2B=Ytheorup(count,6*(j-1)+4)*Rt/Abbulk;
        y3=Ytheorup(count,6*(j-1)+5);
        y3B=Ytheorup(count,6*(j-1)+6);

        Drug_Ab = (y1+y2+y3)*Abbulk/4; %nM
        DrugB_Ab=(y1B+y2B+y3B)*Abbulk/4; %nM

        Drug_Antibodies(end+1,:)=Drug_Ab; 
        DrugB_Antibodies(end+1,:)=DrugB_Ab; 
    end

    for count=1:17
        y1=1*Ytheorcl(count,6*(j-1)+1);
        y1B=1*Ytheorcl(count,6*(j-1)+2);
        y2=1*Ytheorcl(count,6*(j-1)+3)*Rt/Abbulk;
        y2B=1*Ytheorcl(count,6*(j-1)+4)*Rt/Abbulk;
        y3=Ytheorcl(count,6*(j-1)+5);
        y3B=Ytheorcl(count,6*(j-1)+6);

        Drug_Ab = (y1+y2+y3)*Abbulk/4; %nM
        DrugB_Ab = (y1B+y2B+y3B)*Abbulk/4; %nM
        Drug_Antibodies(end+1,:)=Drug_Ab; 
        DrugB_Antibodies(end+1,:)=DrugB_Ab; 
    end
else
    Drug_Antibodies=zeros(65,100);
    DrugB_Antibodies=zeros(65,100);
end


load liposomes_drug
save drug_distribution_coktail_Ab_Lip



% For radio
Drug_Radio=Drug_Antibodies;
DrugB_Radio=DrugB_Antibodies;
Drug_Radio(37:end,:)=Drug_Radio(37:end,:)+Drugfree_intern_Lip+Drug_retained_Lip;
DrugB_Radio(37:end,:)=DrugB_Radio(37:end,:)+DrugBfree_intern_Lip+DrugB_retained_Lip;


Ti=0:0.5:32;

Int_Radio=zeros(size(Drug_Radio));
for i=1:100
for j=2:length(Ti)
Int_Radio(j,i)=trapz(Ti(1:j),Drug_Radio(1:j,i)); % integrated drug (active/unspent) concentration
end
end

save results_radio 'Int_Radio' 'Ti' 'xpt'

% Spatio-temporal evolution of total (active and inactive) drug concentration 
[X,T]=meshgrid(xpt,Ti);
% Generates Figure 1(B) left panel
drug_3d_Ab_50_Lip_50(X*200, T, Drug_Radio+DrugB_Radio, [0 200], [18 18], [0 0])

% Generate figure illustrating the time integrated drug concentration 
% Figure 1(B) right panel
integrated_drug_3d_Ab_50_Lip_50(X*200, T, Int_Radio, [0 200], [18 18], [0 0])

% Generate figure illustrating the survival fraction of cancer cells when
% using a 50/50 Ab/Liposome cocktail 
% this reproduces the orange line in Figure 1 (C) middle panel
survival_fraction(xpt*200,exp(-0.015*Int_Radio(65,:))*100,Lf)

% Compute average survival fraction over all tumor spheroid volume 
x=xpt*200; y=exp(-0.015*Int_Radio(65,:))*100; y=y';
SF_average=trapz(x,x.^2.*y)*3/x(end)^3;

fprintf(['Cocktail: ',num2str(Lf*60), ' nM carried by Liposomes \n'])
fprintf(['              ',num2str((1-Lf)*60), ' nM carried by Antibodies \n'])
fprintf(['Average survival fraction for 200 um spheroid = ', num2str(SF_average),'%% \n'])
fprintf('\n')



    


