% Script generating normalized variables X,Y used for fitting of kinetic
% parameters of non-adhering/non-releasing liposomes

% Load data
T=readtable('Non-Adhering Liposomes_Uptake and Clearance_Final Concentration.xlsx');

% radius 
Ri=T.Var1; 

% Liposome Concentration at the end of incubation time,  6hrs
Y6=T.Average_3;
% Liposome Concentration after 0.5 hrs clearance time
Y05=T.Average_4; 
% Liposome Concentration after 1 hr clearance time
Y1=T.Average_5; 
% Liposome Concentration after 2 hrs clearance time
Y2=T.Average_6;

% Remove NaN entries
Y6=rmmissing(Y6);
Y05=rmmissing(Y05); 
Y1=rmmissing(Y1);
Y2=rmmissing(Y2);

% find minimum data
l=min([length(Y05),length(Y1),length(Y2),length(Y6)]);

Y=[]; X=[];

Ri=Ri(1:l);
% 6 hrs uptake
Y(end+1:end+l)=fliplr(Y6(1:l));
X(end+1:end+l)=fliplr(Ri);
% 0.5 hrs clearance
Y(end+1:end+l)=fliplr(Y05(1:l));
X(end+1:end+l)=fliplr(Ri);
% 1 hr clearance
Y(end+1:end+l)=fliplr(Y1(1:l));
X(end+1:end+l)=fliplr(Ri);
% 2 hrsclearance
Y(end+1:end+l)=fliplr(Y2(1:l));
X(end+1:end+l)=fliplr(Ri);

% Normalize data
X=X/max(Ri); 
Y=Y/0.5; % incubation concentration of liposomes is 0.5 mM
X=X'; Y=Y'; 

save data_fit 'X' 'Y'




