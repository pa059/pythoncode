import numpy as np

refmed= 1.0 # refreactive index of sorrounding medium

re1= 1.59  # refractive index core real
im1= 0.66  # refractive index core imaginary
re2= 1.409 # refractive index coat real
im2= .1747 # refractive index coat imaginary

radcor= 0.171 # radius of core
radcot= 6.265 # radius of coat
wavel= 3 # wavelength

refrel1=complex(re1,im1)/refmed
refrel2=complex(re2,im2)/refmed

x= (2.0 * np.pi * radcor * refmed)/wavel
y= (2.0 * np.pi * radcot * refmed)/wavel

print(refrel1)