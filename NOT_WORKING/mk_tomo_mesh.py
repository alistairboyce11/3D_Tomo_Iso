

"""
Created on Fri Jan 17 2020

@author: ab2568
Script to plot Vp isosurfaces in 3D tomographic model


"""
from pathlib import Path
home = str(Path.home())
#
import sys
#import obspy
from obspy.taup import TauPyModel
import matplotlib.pyplot as plt
import os.path
import glob
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import math
from skimage import measure
from skimage.draw import ellipsoid
from geographiclib.geodesic import Geodesic as geo
from scipy.interpolate import interp1d
from mpl_toolkits.mplot3d.axes3d import Axes3D
from mpl_toolkits.mplot3d import proj3d

taupmodel = TauPyModel(model='ak135') # Could change to AK135 to setup in accordance with BBAFRP19

# from matplotlib import rc
# rc('font',size=28)
# rc('font',family='serif')
# rc('axes',labelsize=32)



sys.path.append(home+'/Google_Drive/GITHUB_AB/3D_corrections/PLOTTING')
import Africa_BBAFRP19


def init_model_BBAFRP19():
    # Read crustal and mantle models
    global mod
    mod=Africa_BBAFRP19.BBAFRP19_model()
    mod.read(home+'/Google_Drive/GITHUB_AB/3D_corrections/MODELS/BBAFRP19/', filenames = ['BBAFRP19_P_model.dat'], verbose=True)
    

def init_model_BBAFRP20():
    # Read crustal and mantle models
    global mod
    mod=Africa_BBAFRP20.BBAFRP20_model()
    mod.read(home+'/Google_Drive/GITHUB_AB/3D_corrections/MODELS/BBAFRP20/', filenames = ['BBAFRP20_P_model.dat'], verbose=True)
    


def plot_BBAFRP19_model_3D():
    
    init_model_BBAFRP19()

    model=getattr(mod,'dVp')

    dVp_mod=np.transpose(model,(2,1,0))

    print('model variable is type: '+str(type(dVp_mod)))
    print('model variable is shape: ' +str(np.shape(dVp_mod)))

    print('Lon interval = '+str(mod.dlon)+', Lat interval = '+str(mod.dlat))
    print('Lon min = '+str(mod.lon_min)+', Lon max = '+str(mod.lon_max))
    print('Lat min = '+str(mod.lat_min)+', Lat max = '+str(mod.lat_max)+'\n')

    min_depth=np.min(mod.depth)
    max_depth=np.max(mod.depth)
    d_depth=mod.depth[1]-mod.depth[0]
    print('Min depth = '+str(min_depth)+', Max depth = '+str(max_depth)+'\n')

    lon=40
    lat=10
    depthrange=np.arange(0,2800,10)
    dvs=mod.get_value(0,lon,lat,'dVs')
    dvp=mod.get_value(0,lon,lat,'dVp')
    # vsref=taupmodel.model.s_mod.v_mod.evaluate_below(0,'s')
    # vpref=taupmodel.model.s_mod.v_mod.evaluate_below(0,'p')
    print('At '+str(lat)+'N,'+str(lon)+'E, dVs = '+str(dvs)+'%, dVp = '+str(dvp)+'%\n')


    begin_depth_labels = 200 #km
    depth_label_interval = 200 # km

    depth_tick_interval=int(depth_label_interval/d_depth)
    depth_tick_start=np.where(mod.depth == begin_depth_labels)[0][0]

    depth_label_pos=np.arange(depth_tick_start,len(mod.depth),depth_tick_interval)
    depth_label_val=np.arange(begin_depth_labels,max_depth,depth_label_interval)

    lat_lon_label_interval=10
    lat_lon_tick_interval=int(lat_lon_label_interval/mod.dlon)
    lon_tick_start=np.where(mod.lon == math.ceil(mod.lon_min/lat_lon_label_interval)*lat_lon_label_interval)[0][0]
    lat_tick_start=np.where(mod.lat == math.ceil(mod.lat_min/lat_lon_label_interval)*lat_lon_label_interval)[0][0]

    lat_label_pos=np.arange(lat_tick_start,len(mod.lat),lat_lon_tick_interval)
    lat_label_val=np.arange(math.ceil(mod.lat_min/lat_lon_label_interval)*lat_lon_label_interval,mod.lat_max,lat_lon_label_interval)
    lon_label_pos=np.arange(lon_tick_start,len(mod.lon),lat_lon_tick_interval)
    lon_label_val=np.arange(math.ceil(mod.lon_min/lat_lon_label_interval)*lat_lon_label_interval,mod.lon_max,lat_lon_label_interval)



    print(min_depth)
    print(max_depth)
    


    # ellip_base = ellipsoid(6, 10, 16, levelset=True)
    # ellip_double = np.concatenate((ellip_base[:-1, ...],ellip_base[2:, ...]), axis=0)
    # print(type(ellip_double))
    # print(np.shape(ellip_double))
    # # Use marching cubes to obtain the surface mesh of these ellipsoids
    # verts, faces, normals, values = measure.marching_cubes_lewiner(ellip_double, 0)

    verts, faces, normals, values = measure.marching_cubes_lewiner(dVp_mod,-0.5,spacing=(1.0, 1.0, 1.0),gradient_direction='descent')

    


    fig = plt.figure(figsize=(10, 10))
    ax1 = fig.add_subplot(111, projection='3d',facecolor='white',frame_on=False)

    # Fancy indexing: `verts[faces]` to generate a collection of triangles
    # mesh = Poly3DCollection(verts[faces])
    # mesh.set_edgecolor('none')
    # ax1.add_collection3d(mesh)

    ax1.plot_trisurf(verts[:, 0], verts[:,1], faces, verts[:, 2], lw=.1, cmap="jet")
    

    ax1.set_xlabel("x-axis: Lon.")
    ax1.set_ylabel("y-axis: Lat")
    # ax1.set_zlabel("z-axis: Depth")

    ax1.set_xlim(mod.lon_min, mod.lon_max)  # a = 6 (times two for 2nd ellipsoid)
    ax1.set_ylim(mod.lat_min, mod.lat_max)  # b = 10
    ax1.set_zlim(0,len(mod.depth))  # c = 16

    ax1.set_xticks(lon_label_pos)
    ax1.set_xticklabels(lon_label_val)

    ax1.set_yticks(lat_label_pos)
    ax1.set_yticklabels(lat_label_val)

    # ax1.set_zticks(depth_label_pos)
    # ax1.set_zticklabels(depth_label_val)
    
    
    # ax3.set_xticks(x_label_pos)
    # ax3.set_xticklabels(x_label_val)
    
    # ax1.view_init(elev=30, azim=135)
    # ax1.grid(False)
    ax1.xaxis.pane.set_edgecolor('black')
    ax1.yaxis.pane.set_edgecolor('black')
    ax1.zaxis.pane.set_edgecolor('black')
    ax1.xaxis.pane.fill = False
    ax1.yaxis.pane.fill = False
    ax1.zaxis.pane.fill = False
    ax1.spines['left'].set_position(('axes', 1))
    ax1.spines['bottom'].set_position(('axes', 1))
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    

    """                                                                                                                                                    
    Scaling is done from here...                                                                                                                           
    """
    x_scale=1
    y_scale=1
    z_scale=2

    scale=np.diag([x_scale, y_scale, z_scale, 1.0])
    scale=scale*(1.0/scale.max())
    scale[3,3]=1.0

    def short_proj():
      return np.dot(Axes3D.get_proj(ax1), scale)

    ax1.get_proj=short_proj
    """                                                                                                                                                    
    to here                                                                                                                                                
    """
    
    
    plt.gca().invert_zaxis()
    plt.tight_layout()
    plt.show()


plot_BBAFRP19_model_3D()




#
# # Generate a level set about zero of two identical ellipsoids in 3D
# ellip_base = ellipsoid(6, 10, 16, levelset=True)
# ellip_double = np.concatenate((ellip_base[:-1, ...],
#                                ellip_base[2:, ...]), axis=0)

#
# # Display resulting triangular mesh using Matplotlib. This can also be done
# # with mayavi (see skimage.measure.marching_cubes_lewiner docstring).

