import numpy as np
import pandas as pd
import bagpipes as bp
import math

filter_list_splus = ['filters/uJAVA.dat',
                     'filters/F0378.dat',
                     'filters/F0395.dat',
                     'filters/F0410.dat',
                     'filters/F0430.dat',
                     'filters/gSDSS.dat',
                     'filters/F0515.dat',
                     'filters/rSDSS.dat',
                     'filters/F0660.dat',
                     'filters/iSDSS.dat',
                     'filters/F0861.dat',
                     'filters/zSDSS.dat']


mags = ["U_PETRO_c", "F378_PETRO_c", "F395_PETRO_c", "F410_PETRO_c", "F430_PETRO_c", "G_PETRO_c", "F515_PETRO_c", "R_PETRO_c",
        "F660_PETRO_c", "I_PETRO_c", "F861_PETRO_c", "Z_PETRO_c"]
mags_err = ["e_U_PETRO", "e_F378_PETRO", "e_F395_PETRO", "e_F410_PETRO", "e_F430_PETRO", "e_G_PETRO",
            "e_F515_PETRO","e_R_PETRO", "e_F660_PETRO","e_I_PETRO", "e_F861_PETRO", "e_Z_PETRO"]

def load_Splus(ID):
    # Find the correct row for the object we want.
    row = int(ID) - 1

    cat = pd.read_csv("6_BCD.csv", delimiter=',', header=0)

    #Get our objects' ID
    objectInfoList = pd.read_csv("6_BCD.csv", delimiter=',', header=0).filter(regex=("ID|RA|DEC|zml"))
    print(f'This objects\' information is \nID: {objectInfoList["ID_"][row]}\nRA: {objectInfoList["RA"][row]}\nDEC: {objectInfoList["DEC"][row]}\nRedshift: {objectInfoList["zml"][row]}')

    # Extract the object we want from the catalogue.
    fluxes = []
    fluxerrs = []
    galaxy_param = cat[(cat['ID_'] == objectInfoList["ID_"][row])]

    for k in range(0, len(mags)):
        m = galaxy_param[mags[k]]
        if (math.isnan(m.values)) | (m.values == np.inf) | (m.values == 99.) | (m.values == -99.):
            f = np.array([99.])
            delta_f = np.array([99.])
        else:
            f = 10**(9.56) * 10**(-m/2.5)  # flux in mJy
            delta_m = galaxy_param[mags_err[k]]
            delta_f = f * (1/2.5) * np.log(10) * delta_m

        fluxes.append(f)  #mJy
        fluxerrs.append(delta_f)  #mJy

    # Turn these into a 2D array.
    photometry = np.c_[fluxes, fluxerrs]

    # blow up the errors associated with any missing fluxes.
    for i in range(len(photometry)):
        if (photometry[i, 0] == 0.) or (photometry[i, 1] <= 0):
            photometry[i,:] = [0., 9.9*10**99.]
            

    for i in range(len(photometry)):
        max_snr = 10.

        if photometry[i, 0]/photometry[i, 1] > max_snr:
            photometry[i, 1] = photometry[i, 0]/max_snr

    return photometry

galaxy = bp.galaxy("1", load_Splus, spectrum_exists=False, filt_list=filter_list_splus)
fig = galaxy.plot()

delayed = {}
delayed["age"] = (0., 15.)
delayed["tau"] = (0., 15.)
delayed["massformed"] = (1., 15.)
delayed["metallicity"] = (0., 2.5)

dust = {}
dust["type"] = "Calzetti"
dust["Av"] = (0., 2.)
dust["eta"] = 2.

nebular = {}
nebular["logU"] = -3.

fit_instructions = {}
fit_instructions["redshift"] = (0.,1)
fit_instructions["exponential"] = delayed
fit_instructions["dust"] = dust
fit_instructions["nebular"] = nebular

fit = bp.fit(galaxy, fit_instructions)

fit.fit(verbose=False)

fig = fit.plot_spectrum_posterior(save=False, show=True)
fig = fit.plot_sfh_posterior(save=False, show=True)