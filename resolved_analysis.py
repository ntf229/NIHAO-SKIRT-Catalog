# Make images from spatially resolved fits files

import numpy as np
import matplotlib.pyplot as plt
import fitsio
import matplotlib.patches as mpatches
import os

def makeImages(galaxy):
    os.system('mkdir -p compositeImages/'+galaxy+'/')
    for i in range(len(flux[:,0,0,0])): # loop over orientations
        # RGB = z, r, g:
        r_grid = flux[i,6,:,:] # sdss_z band (see band_names for indicies)
        g_grid = flux[i,4,:,:] # sdss_r band
        b_grid = flux[i,3,:,:] # sdss_g band
        fuv_grid = flux[i,0,:,:]
        w4_grid = flux[i,13,:,:]
        # set brightest pixels to value of 50th brightest pixel
        numVisable = len(r_grid[r_grid > np.amax(r_grid)*1e-5])
        numBrightest = int(numVisable*0.001)
        max_r = r_grid.flatten()[np.argsort(r_grid.flatten())][-numBrightest]
        max_g = g_grid.flatten()[np.argsort(g_grid.flatten())][-numBrightest]
        max_b = b_grid.flatten()[np.argsort(b_grid.flatten())][-numBrightest]
        max_fuv = fuv_grid.flatten()[np.argsort(fuv_grid.flatten())][-numBrightest]
        max_w4 = w4_grid.flatten()[np.argsort(w4_grid.flatten())][-numBrightest]
        r_grid[r_grid > max_r] = max_r
        g_grid[g_grid > max_g] = max_g
        b_grid[b_grid > max_b] = max_b
        fuv_grid[fuv_grid > max_fuv] = max_fuv
        w4_grid[w4_grid > max_w4] = max_w4
        tot_max = np.amax([np.amax(r_grid),np.amax(g_grid),np.amax(b_grid)])
        stretch_image = 0.005 # increase to make dim pixels dimmer
        stretch_image2 = 0.001
        image_r = np.arcsinh((r_grid/np.amax(r_grid))/stretch_image) 
        image_g = np.arcsinh((g_grid/np.amax(g_grid))/stretch_image) 
        image_b = np.arcsinh((b_grid/np.amax(b_grid))/stretch_image) 
        image_fuv = np.arcsinh((fuv_grid/np.amax(fuv_grid))/stretch_image2) 
        image_w4 = np.arcsinh((w4_grid/np.amax(w4_grid))/stretch_image2) 
        fig = plt.figure()
        sizes = np.shape(r_grid)
        fig.set_size_inches(1. * sizes[0] / sizes[1], 1, forward = False)
        ax = plt.Axes(fig, [0., 0., 1., 1.])
        ax.set_axis_off()
        fig.add_axes(ax)
        red = (image_r*0.70+image_w4*0.3)/np.amax(image_r*0.70+image_w4*0.3)
        green = image_g/np.amax(image_g)*0.75
        blue = (image_b*0.70+image_fuv*0.3)/np.amax(image_b*0.70+image_fuv*0.3)
        image = np.transpose(np.asarray([red,green,blue]))
        ax.imshow(image, interpolation='none')
        plt.savefig('compositeImages/'+galaxy+'/axisRatio'+str(axisRatios[i])+'.png', dpi=sizes[0])
        plt.close()

# Broadband photometric bandpass names
band_names = ['FUV', 'NUV', 'u', 'g', 'r', 'i', 'z', '2MASS_J', '2MASS_H', '2MASS_KS', 'W1',
              'W2', 'W3', 'W4', 'PACS70', 'PACS100', 'PACS160', 'SPIRE250', 'SPIRE350', 'SPIRE500']

# NIHAO galaxy names
names = ['g1.88e10','g1.89e10','g1.90e10','g2.34e10','g2.63e10','g2.64e10','g2.80e10','g2.83e10','g2.94e10',
            'g3.44e10','g3.67e10','g3.93e10','g4.27e10','g4.48e10','g4.86e10','g4.94e10','g4.99e10',
            'g5.05e10','g6.12e10','g6.37e10','g6.77e10','g6.91e10','g7.12e10','g8.89e10','g9.59e10','g1.05e11',
            'g1.08e11','g1.37e11','g1.52e11','g1.57e11','g1.59e11','g1.64e11','g2.04e11','g2.19e11','g2.39e11',
            'g2.41e11','g2.42e11','g2.54e11','g3.06e11','g3.21e11','g3.23e11','g3.49e11','g3.55e11','g3.59e11',
            'g3.71e11','g4.90e11','g5.02e11','g5.31e11','g5.36e11','g5.38e11','g5.46e11','g5.55e11','g6.96e11',
            'g7.08e11','g7.44e11','g7.55e11','g7.66e11','g8.06e11','g8.13e11','g8.26e11','g8.28e11','g1.12e12',
            'g1.77e12','g1.92e12','g2.79e12']

# Make images
for i in range(len(names)):
    galaxies = fitsio.read('ResolvedData/'+names[i]+'_nihao-resolved-photometry.fits', ext='GALAXIES')
    summary = fitsio.read('ResolvedData/'+names[i]+'_nihao-resolved-photometry.fits', ext='SUMMARY')
    bands = summary['bands'][0]
    axisRatios = summary['axis_ratio']
    flux = summary['flux']
    makeImages(names[i])








