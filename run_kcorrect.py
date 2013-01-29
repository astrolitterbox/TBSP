# -*- coding: utf-8 -*-
import numpy as np
import math as math
import matplotlib.pylab as plt
import kcorrect 
from kcorrect.sdss import SDSSFilterList, SDSSPhotoZ, SDSS, SDSSKCorrect
print kcorrect.__file__
from kcorrect.acs import ACSKCorrect, ACSFilterList



from kcorrect.utils.cosmology import ztodm
from kcorrect.utils.spec import lambda_to_centers
from scipy import integrate
from kcorrect.projectiontable import load_vmatrix

import string as string
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import pyfits as pyfits
import utils
import colorpy
from colorpy import ciexyz
from colorpy import colormodels

#import sami_db
import db as db
#import plot_survey as plot
#from sami_new import getDataArray as fullDataArray
#from sami_new import plotScatter as OldSchoolPlot
import matplotlib.font_manager 
import csv

imgDir = 'img/'
dbDir = '../../dev/db/'
dbFile = 'CALIFA.sqlite'
dataFile = 'data/Califa.csv'
morphDataFile = 'morph.dat'
observedList = 'list_observed.txt'
db_dataFile = 'db_data.txt'

#constants
c = 299792.458 
pi = 3.14159265
#     Cosmological parameters
H0 = 70.0 #km/s/Mpc, 1 Mpc= 3.08568025*1e19 km
tH = (3.08568025*1e19/H0)/(3600 * 365 * 24*10e9) #Hubble time in Gyr
dH = c/H0 #in Mpc
Wl = 0.728
Wm = 1 - Wl
Wh = c / H0 
Wk = 1 - Wm - Wl 
tH = (3.08568025*1e19/H0)/(3600 * 365 * 24*10e9) #Hubble time in Gyr



#utilities
def E_z(z): #time derivative of log(a(t)), used for integration
	return 1/(np.sqrt((Wm*(1+z)**3 + Wk*(1+z)**2 + Wl)))

def comovingDistance(z):
	Dc = dH * 1000* integrate.quad(E_z, 0, z)[0] #in kpc
	return Dc

def angular2physical(reff, z): #return physical effective diameter of the galaxy in kpc
    return (math.radians(2*reff/3600) *(comovingDistance(z)) / (1 + z)**2)

def getAbsMag(z, mag, ext):
	print z
	d = comovingDistance(z)
        dm = ztodm(z, (Wm, Wl, H0/100))
	absmag = mag - dm - ext
	print dm, 'dm', absmag, 'absmag'
	return absmag
	
def convert(data):
     tempDATA = []
     for i in data:
         tempDATA.append([float(j) for j in i])
     return np.asarray(tempDATA)

def findOverlap(a, b):
   a_multiset = collections.Counter(list(a))
   b_multiset = collections.Counter(list(b))
   return list((a_multiset & b_multiset).elements())

def read_frame(filename):
    image = pyfits.open(filename)
    #image.info()
    data = image['SCI'].data
    return data

def read_header(fname):
    image = pyfits.open(fname)
    primHead = image[0].header
    head = image['SCI'].header
   
    #print head, 'head', primHead
    #targname = head['TARGNAME']
    exptime = primHead['EXPTIME']
    photflam = head['PHOTFLAM']
    zpt = head['PHOTZPT']
    photplam = head['PHOTPLAM']
    photmode = head['PHOTMODE']
    return {'exptime':exptime, 'photflam':photflam, 'zpt':zpt, 'photplam': photplam, 'photmode': photmode}

def getFrameFilter(photmode):
   if string.find(photmode, 'F606') <> -1:
     return 'F606'
   elif string.find(photmode, 'F814') <> -1:
     return 'F814'
   elif string.find(photmode, 'F475') <> -1:
     return 'F475'
   else:
     print 'filter not found!'
 
 
def getFilterFilename(filterName):
   if filterName == 'F606':
     return 'acs_f606w'
   elif filterName == 'F814':
     return 'acs_f814w'
   elif filterName == 'F475':
     return 'acs_f475w'
   else:
     print 'filter not found!'
     exit()

def countsToMag(data,  frame_params): #as in http://www.stsci.edu/instruments/wfpc2/Wfpc2_hand4/ch8_calibration9.html
    mag = -2.5*np.log10(frame_params['photflam']*data/frame_params['exptime']) + frame_params['zpt']
    print mag.shape
    return mag

class FrameData():
   data = np.empty((10, 10)) 
   
  	
def main():
  filenames = ["data/j8cw52080_drz.fits", "data/j8cw52041_drz.fits", "data/j8cw51011_drz.fits"]
  inputRedshift = 0.022
  for filename in filenames:
  	frameData = FrameData()
	data = read_frame(filename)
  	frame_params = read_header(filename)
  	print frame_params
  	print data.shape
  	#HST filters and their UBVRI (Cousin) counterparts:
  	#U: F336W
  	#B: F439W
  	#V: F555W
  	#R: F675W
  	#I: F814W
  	frame_filter = getFrameFilter(frame_params['photmode'])
  	#due to sky subtraction, I pad non-positive data points w/ 0.001
  	data[np.where(data <= 0)] = 0.001
  	data = data[2000:2010, 2000:2010]
        frameData.data = data
  	magnitudes_array = countsToMag(data, frame_params)
  magnitudes_array = magnitudes_array.flatten()
  print magnitudes_array.shape
  magnitudes = np.hstack(magnitudes_array) 
  print magnitudes.shape
  uncertainties = np.ones((magnitudes.shape))*0.02
  redshift = np.ones((magnitudes.shape))*inputRedshift
  kc = ACSKCorrect(redshift, magnitudes, uncertainties, 0.02, cosmo=(Wm, Wl, H0/100))
  # redshift, stmag, stmag_sigma, extinction,
  lamb, templates = load_vmatrix()
  lamb = lambda_to_centers(lamb)
  print lamb.shape, templates.shape
  useful_wavel = np.where((np.round(lamb, 0) > 3000) & (np.round(lamb,0) < 9000))[0]
  print lamb.shape, 'lambda'
  lamb = lamb[useful_wavel]
  print lamb.shape, 'lambda'
  kcorr = kc.kcorrect()
  coeffs = kc.coeffs#[:, 1:4]
  print coeffs.shape
  coeffs = coeffs
  spec = np.dot(coeffs, templates)
  spec = spec[:, useful_wavel]
  #spec: y pixels, x (10 k) lambda datapoints
  #lamb = lamb[useful_wavel]
  lamb = lamb/10 #convert to Angstroms
  spec = spec*10**17
  print spec.shape, lamb.shape
  #xyz = np.empty((spec.shape[0], 1)) #test
  r = np.empty((spec.shape[0], 1))
  g = r.copy()
  b = r.copy()
  for i in range(0, spec.shape[0]):
  #rgb_vals = ciexyz.xyz_from_spectrum(np.transpose(np.array((lamb, spec[0, :]))))
      print np.array((lamb, spec[i, :])).shape, 'hstack'
      rgb_vals = ciexyz.xyz_from_spectrum(np.transpose(np.array((lamb, spec[i, :]))))

      print rgb_vals, rgb_vals.shape
      r[i], g[i], b[i] = colormodels.irgb_from_xyz(rgb_vals)
    #g = colormodels.irgb_from_xyz(rgb_vals)[1]
    #b = colormodels.irgb_from_xyz(rgb_vals)[2]
    
  
  rgb = np.dstack([r, g, b])
  
  print rgb.shape, 'after dstack'
  rgb = np.reshape(rgb, (10, 10, 3))
  fig = plt.figure()
  plt.imshow(rgb)
  plt.savefig('rgb')
  exit()
  fig = plt.figure()
  plt.plot(lamb[2000:4000], (10**17) * spec[2000:4000])
  plt.savefig('spec.png')
   
  exit()

  exit()

  data = np.empty((939, 16))

  califa_id = db.dbUtils.getFromDB('califa_id', dbDir+'CALIFA.sqlite', 'gc')
  
  u = db.dbUtils.getFromDB('petroMag_u', dbDir+'CALIFA.sqlite', 'mothersample')
  g = db.dbUtils.getFromDB('petroMag_g', dbDir+'CALIFA.sqlite', 'mothersample')
  r = db.dbUtils.getFromDB('petroMag_r', dbDir+'CALIFA.sqlite', 'mothersample')
  i = db.dbUtils.getFromDB('petroMag_i', dbDir+'CALIFA.sqlite', 'mothersample')
  z = db.dbUtils.getFromDB('petroMag_z', dbDir+'CALIFA.sqlite', 'mothersample')


  ext_u = db.dbUtils.getFromDB('extinction_u', dbDir+'CALIFA.sqlite', 'extinction')
  ext_g = db.dbUtils.getFromDB('extinction_g', dbDir+'CALIFA.sqlite', 'extinction')
  ext_r = db.dbUtils.getFromDB('extinction_r', dbDir+'CALIFA.sqlite', 'extinction')
  ext_i = db.dbUtils.getFromDB('extinction_i', dbDir+'CALIFA.sqlite', 'extinction')
  ext_z = db.dbUtils.getFromDB('extinction_z', dbDir+'CALIFA.sqlite', 'extinction')
  
  err_u = db.dbUtils.getFromDB('petroMagErr_u', dbDir+'CALIFA.sqlite', 'extinction')
  err_g = db.dbUtils.getFromDB('petroMagErr_g', dbDir+'CALIFA.sqlite', 'extinction')
  err_r = db.dbUtils.getFromDB('petroMagErr_r', dbDir+'CALIFA.sqlite', 'extinction')
  err_i = db.dbUtils.getFromDB('petroMagErr_i', dbDir+'CALIFA.sqlite', 'extinction')
  err_z = db.dbUtils.getFromDB('petroMagErr_z', dbDir+'CALIFA.sqlite', 'extinction')  
  
  redshift = db.dbUtils.getFromDB('z', dbDir+'CALIFA.sqlite', 'ned_z')  
  data[:, 0] = u[:]
  data[:, 1] = g[:]
  data[:, 2] = r[:]
  data[:, 3] = i[:]
  data[:, 4] = z[:]
  
  data[:, 5] = ext_u[:]
  data[:, 6] = ext_g[:]
  data[:, 7] = ext_r[:]
  data[:, 8] = ext_i[:]
  data[:, 9] = ext_z[:]
  
  data[:, 10] = err_u[:]
  data[:, 11] = err_g[:]
  data[:, 12] = err_r[:]
  data[:, 13] = err_i[:]
  data[:, 14] = err_z[:]
  
  data[:, 15] = redshift[:]
  
  maggies = data[:, 0:5]

  extinction = data[:, 5:10]

  maggies_err = data[:, 10:15] 
  print maggies.shape, extinction.shape, maggies_err.shape
  lamb, templates = load_vmatrix()
  lamb = lambda_to_centers(lamb)
  print lamb.shape, templates.shape

  #outputArray = np.empty((939, 9))
  kc =  SDSSKCorrect(redshift, maggies, maggies_err, extinction, cosmo=(Wm, Wl, H0/100))
  kcorr = kc.kcorrect()

  #absmag = getAbsMag(redshift, maggies[:, 2], extinction[:, 2])#kc.absmag() 
  #outputArray[:,0] = califa_id[:]
  #print kcorr[:, 2][:].shape
  
  #outputArray[:, 1:6] = kc.absmag()  
  coeffs = kc.coeffs#[:, 1:4]
  #print coeffs.shape
  #tmremain = np.array([[0.601525, 0.941511, 0.607033, 0.523732, 0.763937]])
  #ones = np.ones((1, len(redshift)))
  #prod = np.dot(tmremain.T, ones).T 
  print coeffs.shape
  coeffs = coeffs[-1]
  print coeffs
  spec = np.dot(coeffs, templates)
  np.savetxt('spec.txt', spec)
  np.savetxt('lambda.txt', lamb)
  print spec.shape, lamb.shape
  fig = plt.figure()

  plt.plot(lamb[2000:4000], (10**17) * spec[2000:4000])
  plt.savefig('spec.png')
  exit()

  modelMasses = coeffs*prod
  #print modelMasses.shape
  mass = np.sum(modelMasses, axis=1)
  for i in range (0, (len(data))):
    distmod = KC.utils.cosmology.ztodm(redshift[i])
    exp = 10 ** (0.4 * distmod)
    outputArray[i, 6] = mass[i] * exp
    #outputArray[i, 7] = getAbsMag(redshift[i], maggies[i, 2], extinction[i, 2])
    
    outputArray[i, 8] = distmod
  outputArray[:, 7] = kcorr[:, 2]  
  np.savetxt("kcorrect_sdss.csv", outputArray, fmt = '%i, %10.3f, %10.3f, %10.3f, %10.3e, %10.3f, %10.3e, %10.3e, %10.3e')  

  exit()
  data = np.empty((939, 16))

  califa_id = utils.convert(db.dbUtils.getFromDB('califa_id', dbDir+'CALIFA.sqlite', 'gc'))
  
  u = utils.convert(db.dbUtils.getFromDB('el_mag', dbDir+'CALIFA.sqlite', 'u_tot'))
  g = utils.convert(db.dbUtils.getFromDB('el_mag', dbDir+'CALIFA.sqlite', 'g_tot'))
  r = utils.convert(db.dbUtils.getFromDB('el_mag', dbDir+'CALIFA.sqlite', 'r_tot'))
  i = utils.convert(db.dbUtils.getFromDB('el_mag', dbDir+'CALIFA.sqlite', 'i_tot'))
  z = utils.convert(db.dbUtils.getFromDB('el_mag', dbDir+'CALIFA.sqlite', 'z_tot'))


  ext_u = utils.convert(db.dbUtils.getFromDB('extinction_u', dbDir+'CALIFA.sqlite', 'extinction'))
  ext_g = utils.convert(db.dbUtils.getFromDB('extinction_g', dbDir+'CALIFA.sqlite', 'extinction'))
  ext_r = utils.convert(db.dbUtils.getFromDB('extinction_r', dbDir+'CALIFA.sqlite', 'extinction'))
  ext_i = utils.convert(db.dbUtils.getFromDB('extinction_i', dbDir+'CALIFA.sqlite', 'extinction'))
  ext_z = utils.convert(db.dbUtils.getFromDB('extinction_z', dbDir+'CALIFA.sqlite', 'extinction'))
  
  err_u = utils.convert(db.dbUtils.getFromDB('petroMagErr_u', dbDir+'CALIFA.sqlite', 'extinction'))
  err_g = utils.convert(db.dbUtils.getFromDB('petroMagErr_g', dbDir+'CALIFA.sqlite', 'extinction'))
  err_r = utils.convert(db.dbUtils.getFromDB('petroMagErr_r', dbDir+'CALIFA.sqlite', 'extinction'))
  err_i = utils.convert(db.dbUtils.getFromDB('petroMagErr_i', dbDir+'CALIFA.sqlite', 'extinction'))
  err_z = utils.convert(db.dbUtils.getFromDB('petroMagErr_z', dbDir+'CALIFA.sqlite', 'extinction'))  
  
  redshift = utils.convert(db.dbUtils.getFromDB('z', dbDir+'CALIFA.sqlite', 'ned_z'))  
  
  data[:, 0] = u[:, 0]
  data[:, 1] = g[:, 0]
  data[:, 2] = r[:, 0]
  data[:, 3] = i[:, 0]
  data[:, 4] = z[:, 0]
  
  data[:, 5] = ext_u[:, 0]
  data[:, 6] = ext_g[:, 0]
  data[:, 7] = ext_r[:, 0]
  data[:, 8] = ext_i[:, 0]
  data[:, 9] = ext_z[:, 0]
  
  data[:, 10] = err_u[:, 0]
  data[:, 11] = err_g[:, 0]
  data[:, 12] = err_r[:, 0]
  data[:, 13] = err_i[:, 0]
  data[:, 14] = err_z[:, 0]
  
  data[:, 15] = redshift[:, 0]
  
  maggies = data[:, 0:5]

  extinction = data[:, 5:10]

  maggies_err = data[:, 10:15] 
  print maggies.shape, extinction.shape, maggies_err.shape
  
  outputArray = np.empty((939, 9))
  kc =  SDSSKCorrect(redshift, maggies, maggies_err, extinction, cosmo=(Wm, Wl, H0/100))
  kcorr = kc.kcorrect()

  #absmag = getAbsMag(redshift, maggies[:, 2], extinction[:, 2])#kc.absmag() 
  outputArray[:,0] = califa_id[:, 0]
  #print kcorr[:, 2][:].shape
  
  outputArray[:, 1:6] = kc.absmag()  
  coeffs = kc.coeffs#[:, 1:4]
  tmremain = np.array([[0.601525, 0.941511, 0.607033, 0.523732, 0.763937]])
  ones = np.ones((1, len(redshift)))
  prod = np.dot(tmremain.T, ones).T 
  modelMasses = coeffs*prod
  #print modelMasses.shape
  mass = np.sum(modelMasses, axis=1)
  for i in range (0, (len(data))):
    distmod = KC.utils.cosmology.ztodm(redshift[i])
    exp = 10 ** (0.4 * distmod)
    outputArray[i, 6] = mass[i] * exp
    #outputArray[i, 7] = getAbsMag(redshift[i], maggies[i, 2], extinction[i, 2])
    
    outputArray[i, 8] = distmod
  outputArray[:, 7] = kcorr[:, 2]  
  np.savetxt("absmag.csv", outputArray, fmt = '%i, %10.3f, %10.3f, %10.3f, %10.3e, %10.3f, %10.3e, %10.3e, %10.3e')  
  
if __name__ == '__main__':
    main()




