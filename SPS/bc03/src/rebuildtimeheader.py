#!/usr/bin/env python

# Rebuilds headers of BC/CB model fits file
import os
import sys
import warnings
from   astropy.table import Table
from   astropy.io import fits
from   astropy.utils.exceptions import AstropyWarning

# Read fits table with BC models
with fits.open(str(sys.argv[1])) as hdul:
    warnings.simplefilter('ignore', category=AstropyWarning)
   #f = Table.read(hdul,hdu=1)      # hdu = 1 => SED's (luminosity vs. wavelength)
   #p = Table.read(hdul,hdu=2)      # hdu = 2 => galaxy physical properties
   #m = Table.read(hdul,hdu=3)      # hdu = 3 => photometric magnitude in different bands
   #d = Table.read(hdul,hdu=4)      # hdu = 4 => line spectral indices
    t = Table.read(hdul,hdu=5)      # hdu = 5 => time scale for spectral evolution (221 steps)
    t = t['age-yr']                 # time scale without pre MS evolution (agrees with t for BC03 and C&B models)

print("# Build fits compatible header for BINTABLE1 (model age)")
print()
print("sechdr = hdulist[1].header")

# Build fits compatible header for BINTABLE1 (model age)	
sechdr = hdul[1].header	
for i in range(len(t)+1):
    if i==0:
        s = 'Wavelength'
        u = 'A  (Wavelength)'
    else:
        c = t[i-1]
        if c < 10.E9:
            s = str('{:.5E}'.format(c))
        else:
            s = str('{:.4f}'.format(c*1.E-9)) + "E9"
        s = s.replace("E+0", "E")
        u = 'Lo/A  (SED at t = ' + s + ' yr)'
        s = s.replace(".", "p")
        s = 't' + s
    s = "sechdr['TTYPE" + str(i+1) + "']  = '" + s + "'"		# label for column i+1
    u = "sechdr['TUNIT" + str(i+1) + "']  = '" + u + "'"		# units for column i+1
    print(s)
    print(u)

print()
print("# Write fits file")
print("hdulist.writeto('newdata.fits')")
