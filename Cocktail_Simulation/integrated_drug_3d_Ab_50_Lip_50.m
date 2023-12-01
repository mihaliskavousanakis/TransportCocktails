function integrated_drug_3d_Ab_50_Lip_50(x1, y1, Int_Radio, X2, Y2, Z1)
%CREATEFIGURE(x1, y1, Int_Radio1, X2, Y2, Z1)
%  X1:  surface xdata
%  Y1:  surface ydata
%  INT_RADIO1:  surface zdata
%  X2:  vector of plot3 x data
%  Y2:  vector of plot3 y data
%  Z1:  vector of plot3 z data

%  Auto-generated by MATLAB on 08-Aug-2023 12:30:50

% Create figure
figure('InvertHardcopy','off','Color',[1 1 1],'Renderer','painters');
colormap(jet);

% Create axes
axes1 = axes;
hold(axes1,'on');

% Create surf
surf(x1,y1,Int_Radio,'FaceAlpha',0.8);

% Create plot3
plot3(X2,Y2,Z1,'LineWidth',2,'Color',[1 0 0]);

% Create zlabel
zlabel({'Time integrated drug','concentration (nM $\cdot$ hrs)'},...
    'FontSize',26.4,...
    'Interpreter','latex');

% Create ylabel
ylabel('time (hrs)','FontSize',26.4,'Interpreter','latex');

% Create xlabel
xlabel({'distance from ','spheroid surface ($\mu m$)'},'FontSize',26.4,...
    'Interpreter','latex');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[0 200]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[0 36]);
view(axes1,[141.454856512141 34.0688344713915]);
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'FontSize',24,'XDir','reverse','XTick',[0 50 100 150 200 250 300],'YTick',[0 8 16 24 32],...
    'YTickLabel',{'0','8','16','24','32'},'ZTick',[0 0.005 0.01]);
