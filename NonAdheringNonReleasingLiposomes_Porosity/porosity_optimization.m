% Fitting of porosity functional for non-specific / non-releasing Liposomes
% The general form is: 
% phi (r) = x(1) * r^x(2) + x(3)
% with r the normalized radius of the spheroid (normalization is performed
% with the maximum radiues)
clear, clc
% Initial guess for phi functional 
x0=[0.4 5 0.6];

%% Start with the default options
options = optimset;
%% Modify options setting
options = optimset(options,'Display', 'off');
[x,fval,exitflag,output] = ...
fminsearch(@experimental_fitting,x0,options);

% Fitted porosity functional 
fprintf('Porosity functional for non-specific/non-releasing Liposomes \n')
fprintf(['phi=',num2str(x(1)),'*r^',num2str(x(2)),'+',num2str(x(3)),'\n'])
