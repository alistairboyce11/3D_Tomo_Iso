
% Make some tomography style data
% data = 8.*rand(80,80,29)-4;
data = (8.*rand(10,10,10))-4;
data = smooth3(data,'box',1);


patch(isocaps(data,-2,'below'),'FaceColor','interp','EdgeColor','none');
p1 = patch(isosurface(data,-2),'FaceColor','blue','EdgeColor','none');
isonormals(data,p1);
view(3); 
axis vis3d tight
camlight left
colormap('hot');
% daspect([.4,.4,1])
lighting gouraud