% Fitting of transport parameters for non-adhering/non-releasing liposomes 
% The vector of unknowns, x, includes (in dimensionless form): 
% x(1): The diffusion coefficient of liposomes in spheroids 
% x(2): The mass transfer coefficient during uptake, PL,up
% x(3): The mass transfer coefficient during clearance, PL,cl

% Fitting is performed against experimental measurements included in the
% the matlab file, data_fit.mat
% The generation of data_fit.mat is performed in file: data_generation.m

% Load the normalized experimental measurements
% The Y variable represens normalized liposome concentration in spheroids
% at time, t=6 hrs (end of incubation), t=6+0.5 hrs (0.5 hr clearance),
% t=6+1 hrs (1 hr clearance), t=6+2 hrs (2 hrs clearance)

% generate data
data_generation; clear all, close all
load('data_fit.mat');

% Initial guess for fitting parameters
load('nonspecliposomes_res')
x=mean(Xtest); 



opts = statset('nlinfit');
opts.MaxIter=1000; 
[xnew,R,J,CovB]=nlinfit(X,Y,@experimental_fitting_uptake_clearance,x,opts);
ci=nlparci(xnew,R,'covar',CovB); % This is to obtain the 95% confidence intervals for each fitted parameter

% xnew(1) is the dimensionless diffusion coefficient of liposomes in
% spheroids 
% If the radius of the spheroid is R (m), then: 
% DL = xnew(1)*R^2 / 3600 (m^2/s)

% xnew(2) is the dimensionless mass transfer coefficient during uptake
% PL,up = xnew(2)*R / 3600 (m/s)

% xnew(3) is the dimensionless mass transfer coefficient during clerance
% experiments
% PL,cl = xnew(2)*R / 3600 (m/s)

% Generate figures
generate_figures(xnew)
