# TransportCocktails
The folder NonSpecificAntibodies_Porosity includes all files required to reproduce the fitting of the porosity profile adopted for the simulation 
of Antibodies in spheroids. 
The porosity profile can be obtained running the file: porosity_optimization.m and should produce 
phi=0.43948*r^3.1974+0.56052

The folder NonAdheringNonReleasingLiposomes_Porosity includes all files required to reproduce the fitting of the porosity profile adopted for the simulation 
of Liposomes in spheroids. 
The porosity profile can be obtained running the file: porosity_optimization.m and should produce 
phi=0.83209*r^5.1257+0.16791

The folder NonAdheringNonReleasingLiposomes_TransportProperties includes all files required to reproduce the fitting of transport parameters adopted for the simulation of non-releasing/non-adhering liposomes. 
Executing the script: main_script
a) generates the mat file (running data_generation) originating from the raw data measurements (Non-Adhering Liposomes_Uptake and Clearance_Final Concentration.xlsx)
b) executes nonlinear fitting (function nlinfit) to perform the fit of liposome effective diffusion coefficient, and the mass transfer coefficients during uptake and clearance experiments.
The diffusion coefficient can be computed from: xnew(1)*R^2/3600 sec  (R the spheroid radius, 196.25 um)
The uptake mass transfer coefficient is computed from: xnew(2)*R/3600 sec
The clearance mass transfer coefficient is computed from: xnew(3)*R/3600 sec

The folder NonSpecificAntibodies_TransportProperties includes all files required to reproduce the fitting of transport parameters adopted for the simulation of non-specific Antibodies (Rituximab). 
Executing the script: main_script
a) generates the mat file (running data_generation) originating from the raw data measurements (Non-specific_Ab_Rituximab_Uptake and Clearance_Antibody.xlsx)
b) executes nonlinear fitting (function nlinfit) to perform the fit of Ab effective diffusion coefficient, and the mass transfer coefficients during uptake and clearance experiments.
The diffusion coefficient can be computed from: xnew(1)*R^2/3600 sec  (R the spheroid radius, 205 um)
The uptake mass transfer coefficient is computed from: xnew(2)*R/3600 sec
The clearance mass transfer coefficient is computed from: xnew(3)*R/3600 sec

The folder SpecificAntibodies_KineticProperties includes all files required to reproduce the fitting of kinetic properties adopted for the simulation of specific Antibodies (Trastuzumab). 
Executing the script: main_script
a) generate the mat file (running data_generation) originating from the raw data measurements (Specific_ab_Trastuzumab_Uptake and Clearance_Antibody.xlsx)
b) executes nonlinear fitting (function nlinfit) to perform the fit of dissociation rate constant, koff (in hrs^(-1)), dissociation equilibrium constant, KD (in nM), and internalization rate constant, kint (in hrs^(-1)). 

Cocktal_Simulation folder includes all files required to reproduce cocktail carrier simulations 
