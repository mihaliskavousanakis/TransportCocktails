function Ytheor=experimental_fitting_uptake_clearance(x,X)

% x is the parameter vector consisting of: x= [ Deff, Pup, Pcl ], i.e.:
% the effective diffusivity of liposomes, Deff and Pup,Pcl the mass transfer
% coefficient during uptake and clearance experiments, respectively.

%Ytheor corresponds to the vector with normalized liposome concentrations
%as obtained after 6hrs incubation, 0.5 hr, 1 hr, and 2 hrs clearance
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
              
     
     % Uptake experiments
     [xpt,ytheor_up]=nonspecliposome_uptake([x(1) x(2)]);
     % Pick the last solution that corresponds to 6 hrs incubation
     i=4; 
     y0=ytheor_up(i,:); y0=y0';
     % Simulate 1/2 minute fishing
     tspan=[0,0.25/60,0.5/60];
     [xpt,ytheor]=liposomes_clearance([x(1) 1000],y0,tspan); 
     Ytheor(end+1:end+ndata)=interp1(xpt,ytheor(end,:),Ri,'pchip','extrap');
     
     % Clearance experiments
     tspan=[0, 0.5, 1, 2];  
     [xpt,ytheor_cl]=liposomes_clearance([x(1) x(3)],y0, tspan); 
     
     % Simulate 1/2 minute fishing experiment
     for i=1:3
         y0=ytheor_cl(i,:); y0=y0';
         tspan=[0,0.25/60,0.5/60];
         [xpt,ytheor]=liposomes_clearance([x(1) 1000],y0,tspan); 
         Ytheor(end+1:end+ndata)=interp1(xpt,ytheor(end,:),Ri,'pchip','extrap');
     end
     
    Ytheor=Ytheor';
    % Check for non-realistic parameters and set 0 (as penalty for
    % prediction of concentration values)
    for i=1:length(x)
        if x(i)<0
            Ytheor=zeros(size(X)); 
            break
        end
    end
    
catch
    Ytheor=zeros(size(X)); 
end
