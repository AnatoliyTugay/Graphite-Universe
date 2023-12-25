# 1.8.2020. Generation of simulated image.
#
from matplotlib import pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy import units as u
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as ndi
from astropy.visualization import LogStretch
from astropy.modeling import models
import math
import photutils
import time
import statmorph

log_stretch = LogStretch(a=10000.0)
def normalize(image):
    m, M = np.min(image), np.max(image)
    return (image-m) / (M-m)

size = 20  # on each side from the center
sigma_psf = 2.0
ny, nx = 100, 100
y, x = np.mgrid[0:ny, 0:nx]
size = 20  
print('Size is number of pixels on each side from the center')
sigma_psf = 2.0
y, x = np.mgrid[-size:size+1, -size:size+1]
psf = np.exp(-(x**2 + y**2)/(2.0*sigma_psf**2))
psf /= np.sum(psf)

ny, nx = 100, 100
sersic_model = models.Sersic2D(amplitude=10, r_eff=reff, n=1.5, x_0=0.5*nx, y_0=0.5*ny, ellip=0.0, theta=0.5)
image = sersic_model(x, y)
size = 50  # on each side from the center
image2 = ndi.convolve(image, psf)
np.random.seed(j)
snp = 3.0
image2 += (1.0 / snp) * np.random.standard_normal(size=(ny, nx))
gain = 1000.0
threshold = photutils.detect_threshold(image2, 1.5)
npixels = 5  # minimum number of connected pixels
segm = photutils.detect_sources(image2, threshold, npixels)
label = np.argmax(segm.areas) + 1
segmap = segm.data == label
segmap_float = ndi.uniform_filter(np.float64(segmap), size=10)
segmap1 = segmap_float > 0.5

f2=open("flux.txt","a+")
f1=open("255names.txt","r")
f=open("result.txt","a+")
f.write(" ------------- \n")
fake_segmap=np.zeros((101, 101))
round=np.zeros((101, 101))
for i in range(0,101):
 for j in range(0,101):
  fake_segmap[i][j]=False
  round[i][j]=2
  r2=(i-50)*(i-50)+(j-50)*(j-50)
  if r2<225: 
   fake_segmap[i][j]=True
   round[i][j]=math.exp(-r2/50)



#for x in range(1, 600): s=f1.readline()
for x in range(1, 255): # 4459 total
 print('Galaxy number ',x)
 s=f1.readline()
 s1="f"+s[0:9]+".fits"
 if(float(s)>0):
  image_file = get_pkg_data_filename(s1)
  fits.info(image_file)
  image = fits.getdata(image_file, ext=0)
  threshold = photutils.detect_threshold(image, nsigma=2)
  npixels = 5  # minimum number of connected pixels
  segm = photutils.detect_sources(image, threshold, npixels=5)
  ns=2
  while segm==None:
   ns=ns*0.9
   threshold = photutils.detect_threshold(image, nsigma=ns)
   segm = photutils.detect_sources(image, threshold, npixels=5)
   
  label = np.argmax(segm.areas) + 1  
  print('Segmentation areas - ',label)
  segmap = segm.data == label
  segmap_float = ndi.uniform_filter(np.float64(segmap), size=10)
  segmap = segmap_float > 0.5 
  #plt.imshow(segmap, origin='lower', cmap='gray')
  #plt.savefig('s'+s[0:9]+'.png')
  if  (label==1): segmap=fake_segmap
  flux=0
  back=0
  area=0
  for i in range(45,55): 
   for j in range(45,55):
    if (segmap[i][j]==True): area=area+1
    flux=flux+image[i][j]
    back=back+image[i][j+20]


  image1=image
  print('Flux = ',flux)
  if (area<10): segmap=fake_segmap
  if (flux<40): image=image+round
  source_morphs = statmorph.source_morphology(image, segmap, gain=100.0, psf=psf)
  r1=100000
  j=0
  for i in range(0,len(source_morphs)):
   morph = source_morphs[i]
   r=(morph.xc_centroid-50)*(morph.xc_centroid-50)+(morph.xc_centroid-50)*(morph.xc_centroid-50)
   if r<r1:
    r1=r
    j=i


  morph = source_morphs[j]
  fake=False
  a=1
  if(morph.sersic_n==1 or morph.sersic_n==2): 
   fake=True
   segmap=fake_segmap

  while(fake):
   a=a*2
   image=image+a*round
   print('Central pixel value ',image[50][50])
   source_morphs = statmorph.source_morphology(image, segmap, gain=100.0, psf=psf)    
   morph = source_morphs[0]  
   fake=False
   if(morph.sersic_n==1 or morph.sersic_n==2): fake=True
   if(a>100): image=image2+image1/a    

  s2=s[0:9]
  f.write(s2)
  f2.write(s2)
  f2.write(" "+str(flux)+" "+str(back)+" "+str(area)+"\n")
  f.write(" "+str(morph.xc_centroid).rjust(25,' '))
  f.write(" "+str(morph.yc_centroid).rjust(25,' '))
  f.write(" "+str(morph.ellipticity_centroid).rjust(25,' '))
  f.write(" "+str(morph.elongation_centroid).rjust(25,' '))
  f.write(" "+str(morph.orientation_centroid).rjust(25,' '))
  f.write(" "+str(morph.xc_asymmetry).rjust(25,' '))
  f.write(" "+str(morph.yc_asymmetry).rjust(25,' '))
  f.write(" "+str(morph.ellipticity_asymmetry).rjust(25,' '))
  f.write(" "+str(morph.elongation_asymmetry).rjust(25,' '))
  f.write(" "+str(morph.orientation_asymmetry).rjust(25,' '))
  f.write(" "+str(morph.rpetro_circ).rjust(25,' '))
  f.write(" "+str(morph.rpetro_ellip).rjust(25,' '))
  f.write(" "+str(morph.rhalf_circ).rjust(25,' '))
  f.write(" "+str(morph.rhalf_ellip).rjust(25,' '))
  f.write(" "+str(morph.r20).rjust(25,' '))
  f.write(" "+str(morph.r80).rjust(25,' '))
  f.write(" "+str(morph.gini).rjust(25,' '))
  f.write(" "+str(morph.m20).rjust(25,' '))
  f.write(" "+str(morph.gini_m20_bulge).rjust(25,' '))
  f.write(" "+str(morph.gini_m20_merger).rjust(25,' '))
  f.write(" "+str(morph.sn_per_pixel).rjust(25,' '))
  f.write(" "+str(morph.concentration).rjust(25,' '))
  f.write(" "+str(morph.asymmetry).rjust(25,' '))
  f.write(" "+str(morph.smoothness).rjust(25,' '))
  f.write(" "+str(morph.sersic_amplitude).rjust(25,' '))
  f.write(" "+str(morph.sersic_rhalf).rjust(25,' '))
  f.write(" "+str(morph.sersic_n).rjust(25,' '))
  f.write(" "+str(morph.sersic_xc).rjust(25,' '))
  f.write(" "+str(morph.sersic_yc).rjust(25,' '))
  f.write(" "+str(morph.sersic_ellip).rjust(25,' '))
  f.write(" "+str(morph.sersic_theta).rjust(25,' '))
  f.write(" "+str(morph.sky_mean).rjust(15,' '))
  f.write(" "+str(morph.sky_median).rjust(15,' '))
  f.write(" "+str(morph.sky_sigma).rjust(15,' '))
  f.write(" "+str(morph.flag).rjust(5,' '))
  f.write(" "+str(morph.flag_sersic).rjust(5,' '))
  f.write("\n")
  print('xc_centroid =', morph.xc_centroid,' Sersic index = ',morph.sersic_n)

 
f1.close()
f.close()
f2.close()
