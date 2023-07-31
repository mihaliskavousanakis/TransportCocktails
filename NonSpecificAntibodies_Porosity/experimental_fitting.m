function f=experimental_fitting(x)

% Fitting of porosity functional for non-specific Antibodies
% The general form is: 
% phi (r) = x(1) * r^x(2) + x(3)
% with r the normalized radius of the spheroid (normalization is performed
% with the maximum radiues)



% Call the script which performs the simulation 
try
                      
     % We first load the normalized experimental data
     % Radius is normalized with its maximum value, and concentration with
     % the antibody incubation concentration, 0.06 mM
     % In particular we first load the raw experimental measurements: 
     load porosity_data.mat
     %
     % normalize the radius, R with its maximum value
     R=R/max(R); 
     % normalize the concentration with incubation concentration, 0.06
     C24 = C24/0.06; % measured concentration at 24 hrs uptake
     C48 = C48/0.06; % measured concentration at 48 hrs uptake
     C72 = C72/0.06; % measured concentration at 72 hrs uptake
     % The Ab concentration equilibrates after 24 hrs and we compute the
     % average concentration: 
     Cavg = (C24+C48+C72)/3;

     % From the data, we exclude the exterior region of the spheroid: Ab
     % are washed out due to the process of tumor fishing
     % We select the data up to the maximum measured concentration of Ab.
     [k,l]=max(Cavg);
     Rdata = R(l:end); Cdata=Cavg(l:end);

     
     phi=x(1)*Rdata.^x(2)+x(3); 
     f=(sum((phi-Cdata).^2))^(1/2);
     
     phi=x(1)*R.^x(2)+x(3); 
     % Set a penalty if the fitting produces unrealistic values (porosity>1
     % or negative porosity)
     if max(phi)>1 || min(phi)<0
         f=1000; 
     end
     
     
    
catch
    f=1000;
end
