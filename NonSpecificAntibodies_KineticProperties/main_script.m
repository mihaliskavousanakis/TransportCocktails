% Fitting of transport parameters for non-specific antibodies (Rituximab) 
% The vector of unknowns, x, includes (in dimensionless form): 
% x(1): The diffusion coefficient of antibodies in spheroids 
% x(2): The mass transfer coefficient during uptake, PAb,up
% x(3): The mass transfer coefficient during clearance, PAb,cl

% Fitting is performed against experimental measurements included in the
% the matlab file, data_fit.mat
% The generation of data_fit.mat is performed in file: data_generation.m

% Load the normalized experimental measurements
% The Y variable represens normalized antibodies concentration in spheroids
% at time, t=24 hrs (end of incubation), t=24+1 hr (1 hr clearance),
% t=24+2 hrs (2 hrs clearance), t=24+4 hrs (4 hrs clearance)

% generate data
data_generation; clear all, close all
load('data_fit.mat');

% Initial guess for fitting parameters
load('nonspecificAb_res.mat')
x=mean(Xtest); 



opts = statset('nlinfit');
opts.MaxIter=1000; 
[xnew,R,J,CovB]=nlinfit(X,Y,@experimental_fitting_uptake_clearance,x,opts);
ci=nlparci(xnew,R,'covar',CovB); % This is to obtain the 95% confidence intervals for each fitted parameter

% xnew(1) is the dimensionless diffusion coefficient of anitbodies in
% spheroids 
% If the radius of the spheroid is R (m), then: 
% DAb = xnew(1)*R^2 / 3600 (m^2/s)

% xnew(2) is the dimensionless mass transfer coefficient during uptake
% PAb,up = xnew(2)*R / 3600 (m/s)

% xnew(3) is the dimensionless mass transfer coefficient during clerance
% experiments
% PAb,cl = xnew(2)*R / 3600 (m/s)

% Generate figures
generate_figures(xnew)



