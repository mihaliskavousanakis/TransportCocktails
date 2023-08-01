# TransportCocktails
The folder NonSpecificAntibodies_Porosity includes all files required to reproduce the fitting of the porosity profile adopted for the simulation 
of Antibodies in spheroids. 
The porosity profile can be obtained running the file: porosity_optimization.m and should produce 
phi=0.43948*r^3.1974+0.56052

The folder NonAdheringNonReleasingLiposomes_Porosity includes all files required to reproduce the fitting of the porosity profile adopted for the simulation 
of Liposomes in spheroids. 
The porosity profile can be obtained running the file: porosity_optimization.m and should produce 
phi=0.83209*r^5.1257+0.16791

The folder NonAdheringNonReleasingLiposomes_KineticProperties includes all files required to reproduce the fitting of kinetic parameters adopted for the simulation of non-releasing/non-adhering liposomes. 
Executing the script: main_script
a) generates the mat file (running data_generation) originating from the raw data measurements (Non-Adhering Liposomes_Uptake and Clearance_Final Concentration.xlsx)
b) executes nonlinear fitting (function nlnfit) to perform the fit of liposome effective diffusion coefficient, and the mass transfer coefficients during uptake and clearance experiments.
% The diffusion coefficient can be computed from: xnew(1)*R^2/3600 sec  (R the spheroid radius, 196.25 um)
% The uptake mass transfer coefficient is computed from: xnew(2)*R/3600 sec
% The clearance mass transfer coefficient is computed from: xnew(3)*R/3600 sec

Cocktal_Simulation folder includes all files required to reproduce cocktail carrier simulations 
