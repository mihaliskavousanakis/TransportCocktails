% Fitting of kinetic parameters for specific antibodies (Trastuzumab) 
% The vector of unknowns, x, includes: 
% x(1): The dissociation rate constant, koff
% x(2): The equilibrium dissociation constant between the antibody and its antigen, KD
% x(3): The internalization rate constant of the antibody, kint

% Fitting is performed against experimental measurements included in the
% the matlab file, data_fit.mat
% The generation of data_fit.mat is performed in file: data_generation.m

% Load the normalized experimental measurements
% The Y variable represens normalized antibodies concentration in spheroids
% at time, t=24 hrs (end of incubation), t=24+1 hr (1 hr clearance),
% t=24+2 hrs (2 hrs clearance), t=24+4 hrs (4 hrs clearance)
warning off
% generate data
data_generation; clear all, close all
load('data_fit.mat');

% Initial guess for fitting parameters
load('initial_guess.mat')

opts = statset('nlinfit');
opts.MaxIter=1000; 
[xnew,R,J,CovB]=nlinfit(X,Y,@experimental_fitting_uptake_clearance,x,opts);
ci=nlparci(xnew,R,'covar',CovB); % This is to obtain the 95% confidence intervals for each fitted parameter

% Generate figures
generate_figures(xnew)

