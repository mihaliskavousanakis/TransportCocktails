function Ytheor=experimental_fitting_uptake_clearance(x,X)

% x is the parameter vector consisting of: x= [Pup, Pcl,  koff, KD, kint ], i.e.:
% the mass transfer coefficient of specific Ab during uptake
% the mass transfer coefficient of specific Ab during clearance
% the dissociation rate constant, koff
% the dissociation equilibrium constant, KD
% the internalization rate constant, kint.

%Ytheor corresponds to the vector with normalized Ab concentrations
%as obtained after 24hrs incubation, 1 hr, 2 hrs, and 4 hrs clearance
%experiments. 

% Call the script which performs the simulation 
try
     ndata=length(X)/4;
     Ri=X(1:ndata); 
     Ytheor=[];
     % number of nodes 
     np = 100; 
     % Discretize domain
     xpt=zeros(np,1); 
     dx=1/np; 
     xpt(1)=dx/2; 
     for i=2:np
        xpt(i)=xpt(i-1)+dx;
     end
              
     % Load parameters from non-specific Ab simulations 
     load nonspecific_data.mat 
     Deff=xnew(1)*(205/200)^2; 
     % 205 um is the radius of spheroids where fitting was performed for transport properties of non-specific Ab and 200 um is the radius of spheroids 
     % where the fitting for specific Ab is now performed
     
     % load values for concentration of Ab (in nM) in the solution (Abbulk)
     % and Rt the concentration of antibody receptors (in nM)
     load experimental_Abbulk_Rt
     
     hup=x(1); 
     hcl=x(2);
     koff=x(3); 
     KD=x(4); 
     kr=x(5); 
     

     % Uptake experiments
     [xpt,ytheor_up]=specantibody_uptake([Deff hup koff koff/KD Rt Abbulk kr], [0 1 2 4 6 12 16 24]);
     % Pick the last solution that corresponds to 24 hrs incubation and set
     % it as initial condition for the clearance experiments 
     y0=ytheor_up(end,:); y0=y0';
         
     % Clearance experiments
     tspan=[0, 1, 2, 4]; % 0 corresponds to 24 hrs incubation 
     [xpt,ytheor_cl]=antibodies_clearance([Deff hcl koff koff/KD Rt Abbulk kr],y0, tspan); 
     
     % Simulate 1/2 minute fishing experiment
     j=1:np;
     for i=1:4
         y0=ytheor_cl(i,:); y0=y0';
         tspan=[0,0.25/60,0.5/60];
         [xpt,ytheor]=antibodies_clearance([Deff 10000 koff koff/KD Rt Abbulk kr],y0, tspan); 
         y1=interp1(xpt,ytheor(end,3*(j-1)+1),Ri,'pchip','extrap');
         y2=interp1(xpt,ytheor(end,3*(j-1)+2),Ri,'pchip','extrap');
         y3=interp1(xpt,ytheor(end,3*(j-1)+3),Ri,'pchip','extrap');
         Ytheor(end+1:end+length(Ri))=y1+y2*Rt/Abbulk+y3; % this includes all forms of antibodies (interstitium, bound, and internalized forms)
     end
     
    Ytheor=Ytheor';
   % Check for non-realistic parameters and set 0 (as penalty for
    %prediction of concentration values)
    for i=1:length(x)
        if x(i)<0
            Ytheor=zeros(size(X)); 
            break
        end
    end
    if (x(1)/xnew(2)*(205/200)>100) || (x(1)/xnew(2)*(205/200)<0.002)
            Ytheor=zeros(size(X)); 
            return
    end
    if (x(2)/xnew(3)*(205/200)>100) || (x(2)/xnew(3)*(205/200)<0.002)
            Ytheor=zeros(size(X)); 
            return
        end
        load experimental_measurements
        if (x(3)/koff>100) || (x(3)/koff<0.002)
           Ytheor=zeros(size(X)); 
            return
        end
        if (x(4)/KD>100) || (x(4)/KD<0.002)
           Ytheor=zeros(size(X)); 
            return
        end
        if (x(5)/0.323>100) || (x(5)/0.323<0.002)
            Ytheor=zeros(size(X)); 
            return
        end
catch
    Ytheor=zeros(size(X)); 
end



    
     
     
    

     
     
     
 
