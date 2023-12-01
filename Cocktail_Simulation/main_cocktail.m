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


% Radioactivity concentration 3.7 MBq/L
Ac_tot=3.7/(4*log(2)/(9.9*3600*24)*6.022)*10^(-11); %μΜ

% Fraction of drug carried by liposomes
Lf=0.5;

Drug_Lip = 1*Ac_tot*Lf; %uM
% 0.34 MBq in 1 μmol lipid & 100,000 lipids per liposome
Ac_ratioLip = 0.34/(4*log(2)/(9.9*3600*24)*6.022)*10^(-6);


Lipsol = Drug_Lip/Ac_ratioLip; %uM
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
clearvars -except Lf xpt Ac_tot



Drug_Ab = Ac_tot*(1-Lf); %uM 

% 1.87 MBq in 1 mg Cetuximab (152kDa the M.W. of Cetuximab)
Ac_ratioAb =(1.87/(4*log(2)/(9.9*3600*24)*6.022))*152*10^(-11);


if Drug_Ab>0
    load antibodies_experimental.mat
    Abbulk=Drug_Ab/Ac_ratioAb*1000; %nM

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

        Drug_Ab = (y1+y2+y3)*Abbulk*Ac_ratioAb; %nM
        DrugB_Ab=(y1B+y2B+y3B)*Abbulk*Ac_ratioAb; %nM

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

        Drug_Ab = (y1+y2+y3)*Abbulk*Ac_ratioAb; %nM
        DrugB_Ab = (y1B+y2B+y3B)*Abbulk*Ac_ratioAb; %nM
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
drug_3d_Ab_50_Lip_50(300*(max(max(X))-X), T, Drug_Radio+DrugB_Radio, [0 200], [18 18], [0 0])

% Generate figure illustrating the time integrated drug concentration 
% Figure 1(B) right panel
integrated_drug_3d_Ab_50_Lip_50(300*(max(max(X))-X), T, Int_Radio, [0 300], [18 18], [0 0])

% Generate figure illustrating the survival fraction of cancer cells when
% using a 50/50 Ab/Liposome cocktail 
% this reproduces the orange line in Figure 1 (C) middle panel
survival_fraction(300*(max(xpt)-xpt),exp(-400*Int_Radio(65,:))*100,Lf,Ac_tot*1000)

% Compute average survival fraction over all tumor spheroid volume 
x=xpt*300; y=exp(-400*Int_Radio(65,:))*100; y=y';
SF_average=trapz(x,x.^2.*y)*3/x(end)^3;

fprintf(['Cocktail: ',num2str(Lf*Ac_tot*1000), ' nM carried by Liposomes \n'])
fprintf(['              ',num2str((1-Lf)*Ac_tot*1000), ' nM carried by Antibodies \n'])
fprintf(['Average survival fraction for 300 um spheroid = ', num2str(SF_average),'%% \n'])
fprintf('\n')



    


