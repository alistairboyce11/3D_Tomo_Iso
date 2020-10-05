clear

% Script to plot Vp isosurfaces in synthetic 3D tomographic model

load /Users/ab4810/Google_Drive/GITHUB_AB/3D_model_info/AFRP20_MOD.txt

dep=AFRP20_MOD(:,1);
lat=AFRP20_MOD(:,2);
lon=AFRP20_MOD(:,3);
dVp=AFRP20_MOD(:,4);
Vp=AFRP20_MOD(:,5);

fprintf('Done loading\n')
[LAT,LON,DEP] = meshgrid(unique(lat),unique(lon),unique(dep));

dVp_grid = griddata(lat,lon,dep,dVp,LAT,LON,DEP);
% Vp_grid  = griddata(lat,lon,dep,Vp,LAT,LON,DEP);

% % % %
% Attempt a 3D plot of the synthetic model

data = smooth3(dVp_grid,'gaussian', [15,15,15]);
vel_lim=-0.5;
fprintf('Done making grid\n')

f=figure();


%%%%%%%%%%%%%%%%%%%%%%%%%%% Do Figure 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%

% % This deals with where the contours get to the edge of the model domain i.e. surface and CMB
get_caps=isocaps(LAT,LON,DEP,data,vel_lim,'below');
col_1=get_caps.vertices(:,3);
patch(get_caps,'FaceVertexCData',col_1,'FaceColor','interp','EdgeColor','none','FaceAlpha',1);
colormap(flipud(jet));
hold on

% % This part deals wiht the outside edge of the plumes - i.e. the velocity contour
get_surf=isosurface(LAT,LON,DEP,data,vel_lim);
col_2=get_surf.vertices(:,3);
p1 = patch(get_surf,'FaceVertexCData',col_2,'FaceColor','interp','EdgeColor','none','FaceAlpha',1);
colormap(flipud(jet));
isonormals(data,p1);

view([135,30]); 
% view([135,-30]); 
% view([90,30]); 
% view([180,30]); 

axis vis3d tight

% Sort out X labels:
x_ax_labels=linspace(-40,40,5);
xlabel('Latitude')
set(gca, 'XDir','reverse')
xlim([-44,44])
xticks(x_ax_labels)
xticklabels(x_ax_labels)

% Sort out Y labels:
y_ax_labels=linspace(-20,60,5);
ylabel('Longitude')
ylim([-20,64])
yticks(y_ax_labels)
yticklabels(y_ax_labels)

% Sort out Z labels:
dep_ax_labels=linspace(500,2500,5);
zlim([50,2900])
zlabel('Depth (km)')
set(gca, 'ZDir','reverse')
zticks(dep_ax_labels)
zticklabels(dep_ax_labels)

cbh=colorbar;
set(cbh, 'YDir', 'reverse' );
cbh.Ticks = dep_ax_labels;
cbh.TickLabels = {dep_ax_labels};
set(cbh,'position',[.88 .2 .02 .6])

load coastlines

C=zeros(length(coastlon),2);
C(:,1)=(coastlat);
C(:,2)=(coastlon);
C(:,3)=50;

hold on
plot3(C(:,1),C(:,2),C(:,3),'-','Color', 'k','LineWidth',2)

set(gca,'FontSize',14)
set(gca,'box','on')
camlight left
daspect([0.03,0.03,1])

lighting flat

saveas(f,'../3D_model_info/model_3d_plot_1.png','png')

%%%%%%%%%%%%%%%%%%%%%%%%%%% Do Figure 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%
clf()
% % This deals with where the contours get to the edge of the model domain i.e. surface and CMB
get_caps=isocaps(LAT,LON,DEP,data,vel_lim,'below');
col_1=get_caps.vertices(:,3);
patch(get_caps,'FaceVertexCData',col_1,'FaceColor','interp','EdgeColor','none','FaceAlpha',1);
colormap(flipud(jet));
hold on

% % This part deals wiht the outside edge of the plumes - i.e. the velocity contour
get_surf=isosurface(LAT,LON,DEP,data,vel_lim);
col_2=get_surf.vertices(:,3);
p1 = patch(get_surf,'FaceVertexCData',col_2,'FaceColor','interp','EdgeColor','none','FaceAlpha',1);
colormap(flipud(jet));
isonormals(data,p1);

% view([135,30]); 
view([135,-30]); 
% view([90,30]); 
% view([180,30]); 

axis vis3d tight

% Sort out X labels:
x_ax_labels=linspace(-40,40,5);
xlabel('Latitude')
set(gca, 'XDir','reverse')
xlim([-44,44])
xticks(x_ax_labels)
xticklabels(x_ax_labels)

% Sort out Y labels:
y_ax_labels=linspace(-20,60,5);
ylabel('Longitude')
ylim([-20,64])
yticks(y_ax_labels)
yticklabels(y_ax_labels)

% Sort out Z labels:
dep_ax_labels=linspace(500,2500,5);
zlim([50,2900])
zlabel('Depth (km)')
set(gca, 'ZDir','reverse')
zticks(dep_ax_labels)
zticklabels(dep_ax_labels)

cbh=colorbar;
set(cbh, 'YDir', 'reverse' );
cbh.Ticks = dep_ax_labels;
cbh.TickLabels = {dep_ax_labels};
set(cbh,'position',[.88 .2 .02 .6])

load coastlines

C=zeros(length(coastlon),2);
C(:,1)=(coastlat);
C(:,2)=(coastlon);
C(:,3)=50;

hold on
plot3(C(:,1),C(:,2),C(:,3),'-','Color', 'k','LineWidth',2)

set(gca,'FontSize',14)
set(gca,'box','on')
camlight left
daspect([0.03,0.03,1])

lighting flat

saveas(f,'../3D_model_info/model_3d_plot_2.png','png')

%%%%%%%%%%%%%%%%%%%%%%%%%%% Do Figure 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%

clf()
% % This deals with where the contours get to the edge of the model domain i.e. surface and CMB
get_caps=isocaps(LAT,LON,DEP,data,vel_lim,'below');
col_1=get_caps.vertices(:,3);
patch(get_caps,'FaceVertexCData',col_1,'FaceColor','interp','EdgeColor','none','FaceAlpha',1);
colormap(flipud(jet));
hold on

% % This part deals wiht the outside edge of the plumes - i.e. the velocity contour
get_surf=isosurface(LAT,LON,DEP,data,vel_lim);
col_2=get_surf.vertices(:,3);
p1 = patch(get_surf,'FaceVertexCData',col_2,'FaceColor','interp','EdgeColor','none','FaceAlpha',1);
colormap(flipud(jet));
isonormals(data,p1);

% view([135,30]); 
% view([135,-30]); 
% view([90,30]); 
view([180,30]); 

axis vis3d tight

% Sort out X labels:
x_ax_labels=linspace(-40,40,5);
xlabel('Latitude')
set(gca, 'XDir','reverse')
xlim([-44,44])
xticks(x_ax_labels)
xticklabels(x_ax_labels)

% Sort out Y labels:
y_ax_labels=linspace(-20,60,5);
ylabel('Longitude')
ylim([-20,64])
yticks(y_ax_labels)
yticklabels(y_ax_labels)

% Sort out Z labels:
dep_ax_labels=linspace(500,2500,5);
zlim([50,2900])
zlabel('Depth (km)')
set(gca, 'ZDir','reverse')
zticks(dep_ax_labels)
zticklabels(dep_ax_labels)

cbh=colorbar;
set(cbh, 'YDir', 'reverse' );
cbh.Ticks = dep_ax_labels;
cbh.TickLabels = {dep_ax_labels};
set(cbh,'position',[.88 .2 .02 .6])

load coastlines

C=zeros(length(coastlon),2);
C(:,1)=(coastlat);
C(:,2)=(coastlon);
C(:,3)=50;

hold on
plot3(C(:,1),C(:,2),C(:,3),'-','Color', 'k','LineWidth',2)

set(gca,'FontSize',14)
set(gca,'box','on')
camlight left
daspect([0.03,0.03,1])

lighting flat

saveas(f,'../3D_model_info/model_3d_plot_3.png','png')

%%%%%%%%%%%%%%%%%%%%%%%%%%% Do Figure 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%
clf()
% % This deals with where the contours get to the edge of the model domain i.e. surface and CMB
get_caps=isocaps(LAT,LON,DEP,data,vel_lim,'below');
col_1=get_caps.vertices(:,3);
patch(get_caps,'FaceVertexCData',col_1,'FaceColor','interp','EdgeColor','none','FaceAlpha',1);
colormap(flipud(jet));
hold on

% % This part deals wiht the outside edge of the plumes - i.e. the velocity contour
get_surf=isosurface(LAT,LON,DEP,data,vel_lim);
col_2=get_surf.vertices(:,3);
p1 = patch(get_surf,'FaceVertexCData',col_2,'FaceColor','interp','EdgeColor','none','FaceAlpha',1);
colormap(flipud(jet));
isonormals(data,p1);

% view([135,30]); 
% view([135,-30]); 
view([90,30]); 
% view([180,30]); 

axis vis3d tight

% Sort out X labels:
x_ax_labels=linspace(-40,40,5);
xlabel('Latitude')
set(gca, 'XDir','reverse')
xlim([-44,44])
xticks(x_ax_labels)
xticklabels(x_ax_labels)

% Sort out Y labels:
y_ax_labels=linspace(-20,60,5);
ylabel('Longitude')
ylim([-20,64])
yticks(y_ax_labels)
yticklabels(y_ax_labels)

% Sort out Z labels:
dep_ax_labels=linspace(500,2500,5);
zlim([50,2900])
zlabel('Depth (km)')
set(gca, 'ZDir','reverse')
zticks(dep_ax_labels)
zticklabels(dep_ax_labels)

cbh=colorbar;
set(cbh, 'YDir', 'reverse' );
cbh.Ticks = dep_ax_labels;
cbh.TickLabels = {dep_ax_labels};
set(cbh,'position',[.88 .2 .02 .6])

load coastlines

C=zeros(length(coastlon),2);
C(:,1)=(coastlat);
C(:,2)=(coastlon);
C(:,3)=50;

hold on
plot3(C(:,1),C(:,2),C(:,3),'-','Color', 'k','LineWidth',2)
set(gca,'FontSize',14)
set(gca,'box','on')
camlight left
daspect([0.03,0.03,1])

lighting flat

saveas(f,'../3D_model_info/model_3d_plot_4.png','png')

