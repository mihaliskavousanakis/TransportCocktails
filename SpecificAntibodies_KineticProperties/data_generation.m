% Script generating normalized variables X,Y used for fitting of kinetic
% parameters of specific Ab (Trastuzumab)

% Load data
T=readtable('Specific_ab_Trastuzumab_Uptake and Clearance_Antibody.xlsx','Sheet','Sheet2');

% radius 
Ri=T.Var1; 

% Ab Concentration at the end of incubation time,  24hrs
Y24=T.Avg;
% Ab Concentration after 1 hr clearance time
Y1=T.Avg_2; 
% Ab Concentration after 2 hrs clearance time
Y2=T.Avg_3; 
% Ab Concentration after 4 hrs clearance time
Y4=T.Avg_4;

% Remove NaN entries
Y24=rmmissing(Y24);
Y1=rmmissing(Y1); 
Y2=rmmissing(Y2);
Y4=rmmissing(Y4);
Ri=rmmissing(Ri);


% find minimum data
l=min([length(Y1),length(Y2),length(Y4),length(Y24),length(Ri)]);

Y=[]; X=[];

Ri=Ri(1:l);
% 24 hrs uptake
Y(end+1:end+l)=fliplr(Y24(1:l));
X(end+1:end+l)=fliplr(Ri);
% 1 hr clearance
Y(end+1:end+l)=fliplr(Y1(1:l));
X(end+1:end+l)=fliplr(Ri);
% 2 hrs clearance
Y(end+1:end+l)=fliplr(Y2(1:l));
X(end+1:end+l)=fliplr(Ri);
% 4 hrs clearance
Y(end+1:end+l)=fliplr(Y4(1:l));
X(end+1:end+l)=fliplr(Ri);

% Normalize data
X=X/max(Ri); 
Y=Y/0.06; % incubation concentration of Ab is 0.06 uM
X=X'; Y=Y'; 

save data_fit 'X' 'Y'




