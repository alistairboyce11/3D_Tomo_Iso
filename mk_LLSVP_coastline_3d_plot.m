clear

f=figure();
set(gca,'color',[1 1 1])

% load dep_cluster_Vs_CottaarLekic2016.txt
% lon=dep_cluster_Vs_CottaarLekic2016(:,1);
% lat=dep_cluster_Vs_CottaarLekic2016(:,2);
% Vt=dep_cluster_Vs_CottaarLekic2016(:,3);

load dep_cluster_Vp_CottaarLekic2016.txt
lon=dep_cluster_Vp_CottaarLekic2016(:,1);
lat=dep_cluster_Vp_CottaarLekic2016(:,2);
Vt=dep_cluster_Vp_CottaarLekic2016(:,3);
colat=lat+90;
LLVP=zeros(length(colat),3);
LLVP(:,1)=lon;
LLVP(:,2)=colat;
LLVP(:,3)=Vt;
[Z, refvec] = geoloc2grid(colat, lon , Vt,1);

LON=unique(lon);
LON=LON(1:end-1);
COLAT=unique(colat);
COLAT=COLAT(1:end-1);

contour3(LON,COLAT,290-Z,[286,286],'--','Color', 'red','LineWidth',2)

load coastlines

coastcolat=coastlat+90;
C=zeros(length(coastlon),3);
C(:,1)=(coastlon);
C(:,2)=(coastcolat);
hold on

plot3(C(:,1),C(:,2),C(:,3),'-','Color', 'k','LineWidth',2)


view([70,25]); 
axis vis3d tight
xlabel('Longitude')
xlim([-20,60])
ylabel('Latitude')

yticks([50 70 90 110 130])
yticklabels({'-40','-20','0','20','40'})
ylim([50,130])

zlabel('Depth (km)')

zticks([0 100 200])
zticklabels({'0','1000','2000',})
% 
zlim([0,290])
set(gca, 'ZDir','reverse')

set(gca,'FontSize',14)
set(gca,'box','on')
camlight left
daspect([.15,.15,1])


lighting flat
% set(f, 'EdgeColor', 'none');
grid off;
% saveas(f,'../3D_model_info/LLVSP_coastline_3d_plot.eps','eps')
% exportgraphics(f,'../3D_model_info/LLVSP_coastline_3d_plot.pdf','ContentType','vector')

% Use the figure editor to copy the output as a vector and just paste into Illustrator...



