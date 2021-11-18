#this code reads in .mib files and outputs several visualisations of the data
import pyxem as pxm
import numpy as np
import hyperspy.api as hs
import os
import gc
import matplotlib.pyplot as plt 
import matplotlib as p
p.use('Qt5Agg')

from PyQt5 import QtCore, QtWidgets

#from pyxem.signals.diffraction_vectors import DiffractionVectors
system = 'mgstearate'
data_origin_folder = 'F:\\em18488_milled_pharma_sed\\raw\\' + system
destination_folder = "E:\\hl585\\duncan_data_processed\\em18488_milled_pharma_sed\\"+ system + "\\"
filenames = []
for root, dirs, files in os.walk(data_origin_folder):
    for file in files:
        if file.endswith(".mib"):
             filenames.append(os.path.join(root, file))

filenames = np.sort(filenames)

print(filenames)

for fname in filenames:
    print(fname)
    fname_label = fname[44:-5]
    print(fname_label)
    # Load the data
    dp = pxm.load_mib(fname, reshape=True)
    dp.set_signal_type('electron_diffraction')
    print(dp.metadata)
    #dp =pxm.signals.ElectronDiffraction2D(hs.load(fname))
    dp.compute()
    print('finished dp.compute')
    #dp.plot(cmap='inferno', vmax=50)
  
    
    # Find peaks and form diffracting pixels map
    pks = dp.find_peaks(method='difference_of_gaussian',
                        min_sigma=1.,
                        max_sigma=5.,
                        sigma_ratio=2.5,
                        threshold=0.04,
                        overlap=0.5,
                        show_progressbar=False,
                        interactive= False)
    recip_cal = 0.08
    from pyxem.signals.diffraction_vectors import DiffractionVectors
    pks = DiffractionVectors.from_peaks(pks, center=(126,126), calibration = recip_cal)


    xim = pks.get_diffracting_pixels_map(binary=False)
    print('got diffracting pixels map')
    xim.change_dtype('float32')
    #xim.save(fname[:-5] + '_xim.tif')
    xim.save(destination_folder + fname_label+ '_xim.tif')
    xim.save(destination_folder + fname_label+ '_xim.png')
    print('diffraction map saved')

    #now introduce VDFs. Take two rings (one inner, one outer,), both excluding the central beam
    roi1 = hs.roi.CircleROI(cx=126,cy=126,r=50,r_inner=15) #<------- The VDF bit
    vdf1= dp.get_integrated_intensity(roi=roi1)
    roi2 = hs.roi.CircleROI(cx=126,cy=126,r=85,r_inner=50)
    vdf2= dp.get_integrated_intensity(roi=roi2)


    plt.imshow(vdf1.data,vmax=np.max(vdf1.data),cmap='inferno') #<----- I'm using the hspy so I'm having to divide to 
    plt.savefig(destination_folder + fname_label+ 'vdf_15-50pix.png') # get rid of flyback. If using mib shouldn't need.
    plt.close()

    plt.imshow(vdf2.data,vmax=np.max(vdf1.data),cmap='inferno')
    plt.savefig(destination_folder + fname_label+ 'vdf_50-85pix.png')
    plt.close()

    # Calculate average diffraction
    dpeg = dp.mean((0,1))
    dpeg.change_dtype('float32')
    #dpeg.save(fname[:-5] + '_dpeg.tif')
    dpeg.save(destination_folder + fname_label+ '_dpeg.tif')
    dpeg.save(destination_folder + fname_label+ '_dpeg.png')
    print('diffraction pattern saved')

    # Calculate sum image
    mask = np.invert(dp.get_direct_beam_mask(25))
    dp = dp * mask
    dp = pxm.signals.ElectronDiffraction2D(dp)
    im = dp.sum((2,3))
    im = im.as_signal2D((0,1))
    im.change_dtype('float32')
    #im.save(fname[:-5] + '_im.tif')
    im.save(destination_folder + fname_label + '_im.tif')
    im.save(destination_folder + fname_label + '_im.png')
    print('got sum intensity image')

    gc.collect()
