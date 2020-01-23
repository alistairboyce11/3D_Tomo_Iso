clear

% Script to plot Vp isosurfaces in synthetic 3D tomographic model


%%
%Control model output
SP=1;  %Plot SUPERPLUME
VP=1;  %Plot vertical PLUME

GP=0; % Mid mantle veritcal gap in structures. 300,500,700,1000
GSD=0; % GAP seed depth: depth of vertical gap centre in anomalies.

resF=1;
nlon=360*resF;nlat=180*resF;nr=290; % nr smooth at 2900, a little rough at 29. 700-1900

lon = linspace(0,360*pi/180,nlon);
colat = linspace(0,180*pi/180,nlat);
lat = linspace(-90*pi/180,90*pi/180,nlat);
rad = linspace(6371,6371-2900,nr);

[LON,LAT,R] = meshgrid(lon,lat,rad);
A=zeros(size(LON));
VAL=zeros(size(LON));

% PLOT superplume 1 from African LLSVP

if SP == 1
    % Define starting lon, colat, alpha, a, b (ellipse angle, semi major, semi minor axes)
    x0_i=(40)*pi/180;
    x0_f=(15)*pi/180;

    y0_i=(90+10)*pi/180; 
    y0_f=(90+-35)*pi/180; 

    alpha_i=-90*pi/180;
    alpha_f=-75*pi/180;

    % This is distorted because the grid is not square
    % Want it to increase in size with depth.
    a0_i=0.15;
    a0_f=0.3;

    b0_i=0.1;
    b0_f=0.2;

    % This gives straight line varition between alpha x and y.
    % Curvature of superplume downwards - raise x and y int to power.
    curv=1.5;
    lambda_alpha=linspace(0,1,length(rad));
    lambda_x=linspace(0,1,length(rad)).^(1/curv);
    lambda_y=linspace(0,1,length(rad)).^(1/curv);
    lambda_a=linspace(0,1,length(rad));
    lambda_b=linspace(0,1,length(rad));

    x=zeros(1,length(rad));
    y=zeros(1,length(rad));
    alpha=zeros(1,length(rad));
    a=zeros(1,length(rad));
    b=zeros(1,length(rad));

    for k=1:length(rad)
    %     calcualte angle of ellispe for given depth
    %     calculate where centre of ellipse should be in xy space.
        alpha0=      alpha_i +   lambda_alpha(1,k)*(alpha_f - alpha_i);
        x0=         x0_i +      lambda_x(1,k)*(x0_f - x0_i);
        y0=         y0_i +      lambda_y(1,k)*(y0_f - y0_i);
        a0=         a0_i +      lambda_a(1,k)*(a0_f - a0_i);
        b0=         b0_i +      lambda_b(1,k)*(b0_f - b0_i);

        x(k)=x0;
        y(k)=y0;
        alpha(k)=alpha0;
        a(k)=a0;
        b(k)=b0;

    %     Calculate the co-efficients of the ellipse
        C1= (cos(alpha0).^2 / a(k).^2)  +  (sin(alpha0).^2 / b(k).^2) ;
        C2= (sin(alpha0).^2 / a(k).^2)  +  (cos(alpha0).^2 / b(k).^2) ; 
        C3= ((1/ a(k).^2)  + (1 / b(k).^2))*sin(2*alpha0); 

    % Its i,j,k because matrices read (ROW,COLUMN,LAYER)
    % This is (COLAT,LON,LAYER)
        for j=1:length(lon)
            for i=1:length(colat)
                elp = C1*(lon(j) - x0).^2 + C2*(colat(i) - y0).^2 + C3*(lon(j) - x0)*(colat(i) - y0);
                if elp < 1
                    A(i,j,k)=-0.02;
                end
    %             VAL(i,j,k)=elp;
            end
        end
    end
    
clear x0_i y0_i a0_i b0_i alpha_i x0_f y0_f a0_f b0_f alpha_f lambda_a lambda_b lambda_x lambda_y lambda_alpha a b x y alpha C1 C2 C3 x0 y0 alpha0 a0 b0 elp  
end


% PLOT vertical plume 2 from CMB below north Africa.

if VP == 1
    % Define starting lon, colat, alpha, a, b (ellipse angle, semi major, semi minor axes)
    x0_i=(42)*pi/180;
    x0_f=(44)*pi/180;

    y0_i=(90+10)*pi/180; 
    y0_f=(90+13)*pi/180; 

    alpha_i=-90*pi/180;
    alpha_f=-90*pi/180;

    % This is distorted because the grid is not square
    % Want it to increase in size with depth.
%     a0_i=0.1;
%     a0_f=0.2;
    a0_i=0.05;
    a0_f=0.1;


%     b0_i=0.02;
%     b0_f=0.05;
    b0_i=0.04;
    b0_f=0.08;

    % This gives straight line varition between alpha x and y.
    % Curvature of superplume downwards - raise x and y int to power. 1 is
    % vertical
    curv=1;
    lambda_alpha=linspace(0,1,length(rad));
    lambda_x=linspace(0,1,length(rad)).^(1/curv);
    lambda_y=linspace(0,1,length(rad)).^(1/curv);
    lambda_a=linspace(0,1,length(rad));
    lambda_b=linspace(0,1,length(rad));

    x=zeros(1,length(rad));
    y=zeros(1,length(rad));
    alpha=zeros(1,length(rad));
    a=zeros(1,length(rad));
    b=zeros(1,length(rad));

    for k=1:length(rad)
    %     calcualte angle of ellispe for given depth
    %     calculate where centre of ellipse should be in xy space.
        alpha0=      alpha_i +   lambda_alpha(1,k)*(alpha_f - alpha_i);
        x0=         x0_i +      lambda_x(1,k)*(x0_f - x0_i);
        y0=         y0_i +      lambda_y(1,k)*(y0_f - y0_i);
        a0=         a0_i +      lambda_a(1,k)*(a0_f - a0_i);
        b0=         b0_i +      lambda_b(1,k)*(b0_f - b0_i);

        x(k)=x0;
        y(k)=y0;
        alpha(k)=alpha0;
        a(k)=a0;
        b(k)=b0;

    %     Calculate the co-efficients of the ellipse
        C1= (cos(alpha0).^2 / a(k).^2)  +  (sin(alpha0).^2 / b(k).^2) ;
        C2= (sin(alpha0).^2 / a(k).^2)  +  (cos(alpha0).^2 / b(k).^2) ; 
        C3= ((1/ a(k).^2)  + (1 / b(k).^2))*sin(2*alpha0); 

    % Its i,j,k because matrices read (ROW,COLUMN,LAYER)
    % This is (COLAT,LON,LAYER)
        for j=1:length(lon)
            for i=1:length(colat)
                elp = C1*(lon(j) - x0).^2 + C2*(colat(i) - y0).^2 + C3*(lon(j) - x0)*(colat(i) - y0);
                if elp < 1
                    A(i,j,k)=-0.02;
                end
    %             VAL(i,j,k)=elp;
            end
        end
    end
    clear x0_i y0_i a0_i b0_i alpha_i x0_f y0_f a0_f b0_f alpha_f lambda_a lambda_b lambda_x lambda_y lambda_alpha a b x y alpha C1 C2 C3 x0 y0 alpha0 a0 b0 elp
end

% % % % % 
% Attempt a 3D plot of the synthetic model

vel_lim=-2;
data = 100*smooth3(A,'gaussian');

f=figure();


% This deals wiht where the contours get to the edge of the model domain i.e. surface and CMB
get_caps=isocaps(data,vel_lim,'below');
col_1=get_caps.vertices(:,3);
patch(get_caps,'FaceVertexCData',col_1,'FaceColor','interp','EdgeColor','none','FaceAlpha',.8);
colormap('Autumn');
hold on

% This part deals wiht the outside edge of the plumes - i.e. the velocity contour
get_surf=isosurface(data,vel_lim);
col_2=get_surf.vertices(:,3);
p1 = patch(get_surf,'FaceVertexCData',col_2,'FaceColor','interp','EdgeColor','none','FaceAlpha',.8);
colormap('Autumn');
isonormals(data,p1);

view([45,30]); 
axis vis3d tight
xlabel('Longitude')
xlim([-20,60])
ylabel('Co-Latitude')
ylim([50,130])
zlabel('Depth (km)')

zticks([0 50 100 150 200 250])
zticklabels({'0','500','1000','1500','2000','2500'})

zlim([0,290])
set(gca, 'ZDir','reverse')

cbh=colorbar;
set(cbh, 'YDir', 'reverse' );
cbh.Ticks = [1 50 100 150 200 250];
cbh.TickLabels = {'0','500','1000','1500','2000','2500'};
set(cbh,'position',[.88 .2 .02 .6])

load coastlines

coastcolat=coastlat+90;
C=zeros(length(coastlon),3);
C(:,1)=(coastlon);
C(:,2)=(coastcolat);
hold on

plot3(C(:,1),C(:,2),C(:,3),'-','Color', 'k')

set(gca,'FontSize',14)
set(gca,'box','on')
camlight left
daspect([.2,.2,1])
lighting gouraud

saveas(f,'../3D_corrections/synth_3d_plot.png','png')






