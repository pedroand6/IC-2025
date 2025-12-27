#!/usr/bin/env python

# Rebuilds headers of BC/CB model fits file

from astropy.io import fits
import sys

hdulist = fits.open(str(sys.argv[1]))
prihdr = hdulist[0].header
prihdr['stilts0']  = ''
prihdr['stilts']   = 'Fits file created by STILTS v3.1-1, including header up to NTABLE'
prihdr['modlhd0']  = ''
prihdr['modlhdr']  = '---------------- CB2016 model header follows below ----------------'
prihdr['biblio00']  = ''
prihdr['biblio01']  = 'BIBLIOGRAPHY:'
prihdr['biblio02']  = ''
prihdr['biblio03']  = 'EVOLUTIONARY TRACKS:'
prihdr['biblio04']  = '             PARSEC: Bressan, A., et al. 2012, MNRAS, 427'
prihdr['biblio05']  = '                     Chen, Y., et al. 2015, MNRAS, 452, 1068'
prihdr['biblio06']  = '             TP-AGB: Marigo, P. et al. 2013, MNRAS, 434, 488'
prihdr['biblio07']  = ''
prihdr['biblio08']  = 'STELLAR SPECTRA:'
prihdr['biblio09']  = '      Visible range: MILES stellar library'
prihdr['biblio10']  = '                     Sanchez-Blazquez et al. 2006, MNRAS, 371, 703'
prihdr['biblio11']  = '                     Falcon-Barroso, J., et al. 2011, A&A, 532, 95'
prihdr['biblio12']  = '                     Prugniel, Ph. 2011 A&A, 531, 165'
prihdr['biblio13']  = '  Supplemented with:'
prihdr['biblio14']  = '             Stelib: redward of Miles reddest point (7350 A)'
prihdr['biblio15']  = '                     Le Borgne, J.-F. et al. 2003, A&A, 402, 433'
prihdr['biblio16']  = '          BaSeL 3.1: redward of Stelib reddest point (8750 A)'
prihdr['biblio17']  = '                     Westera P. et al. 2002, A&A, 381, 524'
prihdr['biblio18']  = '            O-stars: Tlusty models'
prihdr['biblio19']  = '                     Lanz, T. & Hubeny, I. 2003, ApJS, 146, 417'
prihdr['biblio20']  = '            B-stars: Tlusty models'
prihdr['biblio21']  = '                     Lanz, T. & Hubeny, I. 2007, ApJS, 169, 83'
prihdr['biblio22']  = '   Massive MS stars: WM-Basic models'
prihdr['biblio23']  = '                     Leitherer, C., et al. 2010, ApJS, 189, 309'
prihdr['biblio24']  = '           WR stars: PoWR models'
prihdr['biblio25']  = '                     Hamann, W-R & Gafener, G. 2004, A&A, 427, 697'
prihdr['biblio26']  = '            A-stars: Martins et al. models'
prihdr['biblio27']  = '                     Martins, L.P. et al. 2005, MNRAS, 358, 49'
prihdr['biblio28']  = '   Lower mass stars: UV-Blue models'
prihdr['biblio29']  = '                     Rodriguez-Merino, LH, et al. 2005,ApJ,626,411'
prihdr['biblio30']  = '               CSPN: Rauch models'
prihdr['biblio31']  = '                     Rauch, T. 2002, RevMexAstronAstrof CS, 12, 150'
prihdr['biblio32']  = '       TP-AGB stars: Aringer models for C stars'
prihdr['biblio33']  = '                     Aringer et al. 2009, A&A, 503, 913'
prihdr['biblio34']  = '                     IRTF library'
prihdr['biblio35']  = '                     Rayner et al. 2009, ApJS, 185, 289'
prihdr['biblio36']  = '                     Dusty code'
prihdr['biblio37']  = '                     Ivezic Z. & Elitzur M., 1997, MNRAS, 287, 799'
prihdr['biblio38']  = '                     Gonzalez-L. R. et al. 2010, MNRAS, 403, 1213'
prihdr['biblio39']  = ''
prihdr['biblio40']  = 'LINE INDICES:'
prihdr['biblio41']  = '       Lick indices: Worthey, G., et al. 1994, ApJS, 94, 687'
prihdr['biblio42']  = '                     Trager, S.C., et al. 1998, ApJS, 116, 1'
prihdr['biblio43']  = '      Other indices: Worthey G & Ottaviani DL 1997, ApJS, 111, 377'
prihdr['biblio44']  = '                     Brodie, J., & Hanes, D.A. 1986, ApJ, 300, 258'
prihdr['biblio45']  = '                     Diaz, A. & Terlevich**2 1989, MNRAS, 239, 325'
prihdr['biblio46']  = '                     Gorgas J, Cardiel N et al. 1999, A&AS, 139, 29'
prihdr['biblio47']  = '                     Marcillac, D., et al. 2006, A&A, 458, 369'
prihdr['biblio48']  = '         UV indices: Fanelli, M. N., et al. 1987, ApJ, 321, 768'
prihdr['biblio49']  = '                     Fanelli, M. N., et al. 1990, ApJ, 364, 272'
prihdr['biblio50']  = '                     Fanelli, M. N., et al. 1992, ApJS, 82, 197'
prihdr['biblio51']  = '                     Chavez, M. et al. 2009, ApJ, 700, 694'
prihdr['biblio52']  = '                     Maraston, C. et al. 2009, A&A, 493, 425'
prihdr['biblio53']  = ''
prihdr['biblio54']  = 'Fuel Consumption Theorem: Renzini & Buzzoni 1986, in'
prihdr['biblio55']  = '  Spectral Evolution of Galaxies, eds. C. Chiosi & A. Renzini, p195'
prihdr['biblio56']  = ''
prihdr['biblio57']  = '-------------------------------------------------------------------'
prihdr['modlhd1']  =   ''
prihdr['modlid']   = 'CB2016 - Simple Stellar Population Synthesis Model'
prihdr['modlref']  = 'Charlot & Bruzual (2016, MNRAS, in preparation)'
prihdr['modlfile'] = str(sys.argv[1])

# BC2003 models, Padova 1994 tracks + VWPAGB
if str(sys.argv[2]) == '-22':
	prihdr['modltrks'] = ('0.7696, 0.2303, 0.0001   ','X, Y, Z, Padova 1994 + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-32':
	prihdr['modltrks'] = ('0.7686, 0.2310, 0.0004   ','X, Y, Z, Padova 1994 + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-42':
	prihdr['modltrks'] = ('0.7560, 0.2400, 0.0040   ','X, Y, Z, Padova 1994 + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-52':
	prihdr['modltrks'] = ('0.7420, 0.2500, 0.0080   ','X, Y, Z, Padova 1994 + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-62':
	prihdr['modltrks'] = ('0.7000, 0.2800, 0.0200   ','X, Y, Z, Padova 1994 + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-72':
	prihdr['modltrks'] = ('0.5980, 0.3520, 0.0500   ','X, Y, Z, Padova 1994 + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-82':
	prihdr['modltrks'] = ('0.4250, 0.4750, 0.1000   ','X, Y, Z, Padova 1994 + TP-AGB + VWPAGB')

# BC2003 models, Padova 2000 tracks + VWPAGB
elif str(sys.argv[2]) == '-122':
	prihdr['modltrks'] = ('0.7696, 0.2300, 0.0004   ','X, Y, Z, Padova 2000 + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-132':
	prihdr['modltrks'] = ('0.7690, 0.2300, 0.0010   ','X, Y, Z, Padova 2000 + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-142':
	prihdr['modltrks'] = ('0.7560, 0.2400, 0.0040   ','X, Y, Z, Padova 2000 + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-152':
	prihdr['modltrks'] = ('0.7420, 0.2500, 0.0080   ','X, Y, Z, Padova 2000 + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-162':
	prihdr['modltrks'] = ('0.7080, 0.2730, 0.0190   ','X, Y, Z, Padova 2000 + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-172':
	prihdr['modltrks'] = ('0.6700, 0.3000, 0.0300   ','X, Y, Z, Padova 2000 + TP-AGB + VWPAGB')

# CB2007 models, Padova 1994 tracks + VWPAGB + enhanced TP-AGB
elif str(sys.argv[2]) == '-92':
	prihdr['modltrks'] = ('0.7696, 0.2303, 0.0001   ','X, Y, Z, Padova 1994 + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-93':
	prihdr['modltrks'] = ('0.7686, 0.2310, 0.0004   ','X, Y, Z, Padova 1994 + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-94':
	prihdr['modltrks'] = ('0.7560, 0.2400, 0.0040   ','X, Y, Z, Padova 1994 + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-95':
	prihdr['modltrks'] = ('0.7420, 0.2500, 0.0080   ','X, Y, Z, Padova 1994 + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-96':
	prihdr['modltrks'] = ('0.7000, 0.2800, 0.0200   ','X, Y, Z, Padova 1994 + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-97':
	prihdr['modltrks'] = ('0.5980, 0.3520, 0.0500   ','X, Y, Z, Padova 1994 + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-98':
	prihdr['modltrks'] = ('0.4250, 0.4750, 0.1000   ','X, Y, Z, Padova 1994 + TP-AGB + VWPAGB')

# BC2019 models, PARSEC tracks + VWPAGB
elif str(sys.argv[2]) == '-400':
	prihdr['modltrks'] = ('0.7700, 0.2300, 0.0000   ','X, Y, Z, PARSEC + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-401':
	prihdr['modltrks'] = ('0.7509, 0.2490, 0.0001   ','X, Y, Z, PARSEC + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-402':
	prihdr['modltrks'] = ('0.7508, 0.2490, 0.0002   ','X, Y, Z, PARSEC + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-403':
	prihdr['modltrks'] = ('0.7505, 0.2490, 0.0005   ','X, Y, Z, PARSEC + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-404':
	prihdr['modltrks'] = ('0.7490, 0.2500, 0.0010   ','X, Y, Z, PARSEC + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-405':
	prihdr['modltrks'] = ('0.7460, 0.2520, 0.0020   ','X, Y, Z, PARSEC + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-406':
	prihdr['modltrks'] = ('0.7400, 0.2560, 0.0040   ','X, Y, Z, PARSEC + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-407':
	prihdr['modltrks'] = ('0.7350, 0.2590, 0.0060   ','X, Y, Z, PARSEC + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-408':
	prihdr['modltrks'] = ('0.7290, 0.2630, 0.0080   ','X, Y, Z, PARSEC + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-409':
	prihdr['modltrks'] = ('0.7230, 0.2670, 0.0100   ','X, Y, Z, PARSEC + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-410':
	prihdr['modltrks'] = ('0.7130, 0.2730, 0.0140   ','X, Y, Z, PARSEC + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-411':
	prihdr['modltrks'] = ('0.7040, 0.2790, 0.0170   ','X, Y, Z, PARSEC + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-412':
	prihdr['modltrks'] = ('0.6960, 0.2840, 0.0200   ','X, Y, Z, PARSEC + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-413':
	prihdr['modltrks'] = ('0.6680, 0.3020, 0.0300   ','X, Y, Z, PARSEC + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-414':
	prihdr['modltrks'] = ('0.6390, 0.3210, 0.0400   ','X, Y, Z, PARSEC + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-416':
	prihdr['modltrks'] = ('0.5840, 0.3560, 0.0600   ','X, Y, Z, PARSEC + TP-AGB + VWPAGB')

# CB2019 models, PARSEC tracks + MBPAGB
elif str(sys.argv[2]) == '-500':
	prihdr['modltrks'] = ('0.7700, 0.2300, 0.0000   ','X, Y, Z, PARSEC + TP-AGB + MBPAGB')
elif str(sys.argv[2]) == '-501':
	prihdr['modltrks'] = ('0.7509, 0.2490, 0.0001   ','X, Y, Z, PARSEC + TP-AGB + MBPAGB')
elif str(sys.argv[2]) == '-502':
	prihdr['modltrks'] = ('0.7508, 0.2490, 0.0002   ','X, Y, Z, PARSEC + TP-AGB + MBPAGB')
elif str(sys.argv[2]) == '-503':
	prihdr['modltrks'] = ('0.7505, 0.2490, 0.0005   ','X, Y, Z, PARSEC + TP-AGB + MBPAGB')
elif str(sys.argv[2]) == '-504':
	prihdr['modltrks'] = ('0.7490, 0.2500, 0.0010   ','X, Y, Z, PARSEC + TP-AGB + MBPAGB')
elif str(sys.argv[2]) == '-505':
	prihdr['modltrks'] = ('0.7460, 0.2520, 0.0020   ','X, Y, Z, PARSEC + TP-AGB + MBPAGB')
elif str(sys.argv[2]) == '-506':
	prihdr['modltrks'] = ('0.7400, 0.2560, 0.0040   ','X, Y, Z, PARSEC + TP-AGB + MBPAGB')
elif str(sys.argv[2]) == '-507':
	prihdr['modltrks'] = ('0.7350, 0.2590, 0.0060   ','X, Y, Z, PARSEC + TP-AGB + MBPAGB')
elif str(sys.argv[2]) == '-508':
	prihdr['modltrks'] = ('0.7290, 0.2630, 0.0080   ','X, Y, Z, PARSEC + TP-AGB + MBPAGB')
elif str(sys.argv[2]) == '-509':
	prihdr['modltrks'] = ('0.7230, 0.2670, 0.0100   ','X, Y, Z, PARSEC + TP-AGB + MBPAGB')
elif str(sys.argv[2]) == '-510':
	prihdr['modltrks'] = ('0.7130, 0.2730, 0.0140   ','X, Y, Z, PARSEC + TP-AGB + MBPAGB')
elif str(sys.argv[2]) == '-511':
	prihdr['modltrks'] = ('0.7040, 0.2790, 0.0170   ','X, Y, Z, PARSEC + TP-AGB + MBPAGB')
elif str(sys.argv[2]) == '-512':
	prihdr['modltrks'] = ('0.6960, 0.2840, 0.0200   ','X, Y, Z, PARSEC + TP-AGB + MBPAGB')
elif str(sys.argv[2]) == '-513':
	prihdr['modltrks'] = ('0.6680, 0.3020, 0.0300   ','X, Y, Z, PARSEC + TP-AGB + MBPAGB')
elif str(sys.argv[2]) == '-514':
	prihdr['modltrks'] = ('0.6390, 0.3210, 0.0400   ','X, Y, Z, PARSEC + TP-AGB + MBPAGB')
elif str(sys.argv[2]) == '-516':
	prihdr['modltrks'] = ('0.5840, 0.3560, 0.0600   ','X, Y, Z, PARSEC + TP-AGB + MBPAGB')

# BC2020 models, SISSA 2020 tracks + VWPAGB
elif str(sys.argv[2]) == '-600':
	prihdr['modltrks'] = ('0.7700, 0.2300, 0.0000   ','X, Y, Z, SISSA 2020 + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-601':
	prihdr['modltrks'] = ('0.7509, 0.2490, 0.0001   ','X, Y, Z, SISSA 2020 + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-602':
	prihdr['modltrks'] = ('0.7508, 0.2490, 0.0002   ','X, Y, Z, SISSA 2020 + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-603':
	prihdr['modltrks'] = ('0.7505, 0.2490, 0.0005   ','X, Y, Z, SISSA 2020 + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-604':
	prihdr['modltrks'] = ('0.7490, 0.2500, 0.0010   ','X, Y, Z, SISSA 2020 + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-605':
	prihdr['modltrks'] = ('0.7460, 0.2520, 0.0020   ','X, Y, Z, SISSA 2020 + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-606':
	prihdr['modltrks'] = ('0.7400, 0.2560, 0.0040   ','X, Y, Z, SISSA 2020 + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-607':
	prihdr['modltrks'] = ('0.7350, 0.2590, 0.0060   ','X, Y, Z, SISSA 2020 + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-608':
	prihdr['modltrks'] = ('0.7290, 0.2630, 0.0080   ','X, Y, Z, SISSA 2020 + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-609':
	prihdr['modltrks'] = ('0.7230, 0.2670, 0.0100   ','X, Y, Z, SISSA 2020 + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-610':
	prihdr['modltrks'] = ('0.7130, 0.2730, 0.0140   ','X, Y, Z, SISSA 2020 + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-611':
	prihdr['modltrks'] = ('0.7040, 0.2790, 0.0170   ','X, Y, Z, SISSA 2020 + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-612':
	prihdr['modltrks'] = ('0.6960, 0.2840, 0.0200   ','X, Y, Z, SISSA 2020 + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-613':
	prihdr['modltrks'] = ('0.6680, 0.3020, 0.0300   ','X, Y, Z, SISSA 2020 + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-614':
	prihdr['modltrks'] = ('0.6390, 0.3210, 0.0400   ','X, Y, Z, SISSA 2020 + TP-AGB + VWPAGB')
elif str(sys.argv[2]) == '-616':
	prihdr['modltrks'] = ('0.5840, 0.3560, 0.0600   ','X, Y, Z, SISSA 2020 + TP-AGB + VWPAGB')

# CB2020 models, SISSA 2020 tracks + MBPAGB
elif str(sys.argv[2]) == '-700':
	prihdr['modltrks'] = ('0.7700, 0.2300, 0.0000   ','X, Y, Z, SISSA 2020 + TP-AGB + MBPAGB')
elif str(sys.argv[2]) == '-701':
	prihdr['modltrks'] = ('0.7509, 0.2490, 0.0001   ','X, Y, Z, SISSA 2020 + TP-AGB + MBPAGB')
elif str(sys.argv[2]) == '-702':
	prihdr['modltrks'] = ('0.7508, 0.2490, 0.0002   ','X, Y, Z, SISSA 2020 + TP-AGB + MBPAGB')
elif str(sys.argv[2]) == '-703':
	prihdr['modltrks'] = ('0.7505, 0.2490, 0.0005   ','X, Y, Z, SISSA 2020 + TP-AGB + MBPAGB')
elif str(sys.argv[2]) == '-704':
	prihdr['modltrks'] = ('0.7490, 0.2500, 0.0010   ','X, Y, Z, SISSA 2020 + TP-AGB + MBPAGB')
elif str(sys.argv[2]) == '-705':
	prihdr['modltrks'] = ('0.7460, 0.2520, 0.0020   ','X, Y, Z, SISSA 2020 + TP-AGB + MBPAGB')
elif str(sys.argv[2]) == '-706':
	prihdr['modltrks'] = ('0.7400, 0.2560, 0.0040   ','X, Y, Z, SISSA 2020 + TP-AGB + MBPAGB')
elif str(sys.argv[2]) == '-707':
	prihdr['modltrks'] = ('0.7350, 0.2590, 0.0060   ','X, Y, Z, SISSA 2020 + TP-AGB + MBPAGB')
elif str(sys.argv[2]) == '-708':
	prihdr['modltrks'] = ('0.7290, 0.2630, 0.0080   ','X, Y, Z, SISSA 2020 + TP-AGB + MBPAGB')
elif str(sys.argv[2]) == '-709':
	prihdr['modltrks'] = ('0.7230, 0.2670, 0.0100   ','X, Y, Z, SISSA 2020 + TP-AGB + MBPAGB')
elif str(sys.argv[2]) == '-710':
	prihdr['modltrks'] = ('0.7130, 0.2730, 0.0140   ','X, Y, Z, SISSA 2020 + TP-AGB + MBPAGB')
elif str(sys.argv[2]) == '-711':
	prihdr['modltrks'] = ('0.7040, 0.2790, 0.0170   ','X, Y, Z, SISSA 2020 + TP-AGB + MBPAGB')
elif str(sys.argv[2]) == '-712':
	prihdr['modltrks'] = ('0.6960, 0.2840, 0.0200   ','X, Y, Z, SISSA 2020 + TP-AGB + MBPAGB')
elif str(sys.argv[2]) == '-713':
	prihdr['modltrks'] = ('0.6680, 0.3020, 0.0300   ','X, Y, Z, SISSA 2020 + TP-AGB + MBPAGB')
elif str(sys.argv[2]) == '-714':
	prihdr['modltrks'] = ('0.6390, 0.3210, 0.0400   ','X, Y, Z, SISSA 2020 + TP-AGB + MBPAGB')
elif str(sys.argv[2]) == '-716':
	prihdr['modltrks'] = ('0.5840, 0.3560, 0.0600   ','X, Y, Z, SISSA 2020 + TP-AGB + MBPAGB')

prihdr['modlimf']  = ('0.10, 100                ','[Mo], ML, MU, Chabrier (2003) IMF')
prihdr['modltstp'] = ('1E4, 15E9, 221           ','[yr], Initial, Final age, Nsteps')
prihdr['modlwavl'] = ('   5.6, 3.6E+8, var,13216','[A], Initial, Final wavl, Step, Npoints')
prihdr['modlstl1'] = ('   5.6,  911.0, 0.9, 1007','Tlusty,Martins,UVBlue,WMBasic,PoWR,Rauch')
prihdr['modlstl2'] = (' 911.5, 3540.5, 0.5, 5259','Tlusty,Martins,UVBlue,WMBasic,PoWR,Rauch')
prihdr['modlstl3'] = ('3541.4, 7349.3, 0.9, 4232','Miles')
prihdr['modlstl4'] = ('7351.0, 8750.0, 1.0, 1400','Stelib')
prihdr['modlstl5'] = ('8770.0, 3.6E+8, var, 1318','BaSeL3.1, Aringer+, IRTF, Dusty TP-AGB')
prihdr['modlbnt0']  =   ''
prihdr['modlbnt1'] = ('BINTABLE 1                    ','SED of SSP vs. age in Lo/A')
prihdr['modlbt1a'] = ('  SED(t)                      ','SED for 221 values of age t in Lo/A')
prihdr['modlbt1b'] = ('  SED(t=0) = SED(t=10 Myr)    ','First time step, n=1')
prihdr['modlbt1c'] = ('  SED(t=10 Myr)               ','Second time step, n=2')
prihdr['modlbt1d'] = ('  SED(t=15 Gyr)               ','Last time step, n=221')
prihdr['modlbt20']  =   ''
prihdr['modlbnt2'] = ('BINTABLE 2                    ','PHYSICAL PROPERTIES OF SSP vs. AGE')
prihdr['modlbt2a'] = ('  logNLy, logNHeI, logNHeII   ','H, HeI, and HeII ionizing photons')
prihdr['modlbt2b'] = ('  B912                        ','Amplitude of Lyman break')
prihdr['modlbt2c'] = ('  B4000VN, B4000SDSS, B4000   ','Amplitude of 4000 A break')
prihdr['modlbt2d'] = ('  SNR_yr                      ','Super nova rate /yr')
prihdr['modlbt2e'] = ('  PNBR_yr                     ','Planetary Nebula birth rate /yr')
prihdr['modlbt2f'] = ('  NBH                         ','Number of Black Holes formed')
prihdr['modlbt2g'] = ('  NNS                         ','Number of Neutron Stars formed')
prihdr['modlbt2h'] = ('  NWD                         ','Number of White Dwarfs formed')
prihdr['modlbt2i'] = ('  Lx(Lo)                      ','X ray luminosity (Lo)')
prihdr['modlbt2j'] = ('  Mstars                      ','Mass in living stars in Mo')
prihdr['modlbt2k'] = ('  Mremnants                   ','Mass in stellar remnants in Mo')
prihdr['modlbt2l'] = ('  Mretgas                     ','Mass of gas returned to ISM in Mo')
prihdr['modlbt2m'] = ('  Mgalaxy                     ','Mass of galaxy Mo')
prihdr['modlbt2n'] = ('  SFR_yr                      ','Star Formation Rate in Mo/yr')
prihdr['modlbt2o'] = ('  Mtot                        ','Mstars + Mremnants in Mo')
prihdr['modlbt2p'] = ('  Mtot_Lb, Mtot_Lv, Mtot_Lk   ','Mtot/light ratio in B,V,K bands')
prihdr['modlbt2q'] = ('  Mliv_Lb, Mliv_Lv, Mliv_Lk   ','Mstars/light ratio in B,V,K bands')
prihdr['modlbt2r'] = ('  bt_yr                       ','Evolutionary flux,Renzini & Buzzoni')
prihdr['modlbt2s'] = ('  Bt_yr_Lo                    ','b(t)/Bolometric luminosity.')
prihdr['modlbt2t'] = ('  TurnoffMass                 ','Turnoff Mass in Mo')
prihdr['modlbt2u'] = ('  BPMS_BMS                    ','Ratio of PostMS/MS bolometric flux')
prihdr['modlbt30']  =   ''
prihdr['modlbnt3'] = ('BINTABLE 3                    ','MAGNITUDES vs. AGE')
prihdr['modlbt3a'] = ('  VEGA mags:                  ','Vega magnitude system')
prihdr['modlbt3b'] = ('     Mbol:                    ','Bolometric magnitude')
prihdr['modlbt3c'] = ('     UBVRIJKL:                ','Johnson filters')
prihdr['modlbt3d'] = ('     RI:                      ','Cousins filters')
prihdr['modlbt3e'] = ('     JHK:                     ','Palomar filters')
prihdr['modlbt3f'] = ('     Kprime:                  ','Cowie Kprime filter')
prihdr['modlbt3g'] = ('     JHKs:                    ','2Mass filters')
prihdr['modlbt3h'] = ('     3.6,4.5,5.7,7.9 micron:  ','IRAC filters')
prihdr['modlbt3i'] = ('     12,25,60,100 micron:     ','IRAS filters')
prihdr['modlbt3j'] = ('     24,70,160 micron:        ','MIPS filters')
prihdr['modlbt3k'] = ('     F220w,F250w,F330w,F410w, ','HST ACS WFC wide filters')
prihdr['modlbt3l'] = ('     F435w,F475w,F555w,F606w, ',' "   "   "   "      "   ')
prihdr['modlbt3m'] = ('     F625w,F775w,F814w:       ',' "   "   "   "      "   ')
prihdr['modlbt3n'] = ('     f225w,f275w,f336w,f438w, ','HST UVIS1 filters')
prihdr['modlbt3o'] = ('     f547m,f555w,f606w,f625w, ',' "    "      "   ')
prihdr['modlbt3p'] = ('     f656n,f657n,f658n,f814w: ',' "    "      "   ')
prihdr['modlbt3q'] = ('  AB mags:                    ','AB magnitude system')
prihdr['modlbt3r'] = ('     ugriz:                   ','SDSS filters')
prihdr['modlbt3s'] = ('     ugriyzKs:                ','CFHT MegaCam filters')
prihdr['modlbt3t'] = ('     FUV,NUV:                 ','GALEX filters')
prihdr['modlbt40']  =   ''
prihdr['modlbnt4'] = ('BINTABLE 4                    ','LINE STRENGTH INDICES vs. AGE')
prihdr['modlbt4a'] = ('  CN1,CN2,Ca4227,G4300,Fe4383,','Lick Indices')
prihdr['modlbt4b'] = ('  Ca4455,Fe4531,Fe4668,Hbeta, ','  "     "')
prihdr['modlbt4c'] = ('  Fe5015,Mg1,Mg2,Mgb,Fe5270,  ','  "     "')
prihdr['modlbt4d'] = ('  Fe5335,Fe5406,Fe5709,Fe5782,','  "     "')
prihdr['modlbt4e'] = ('  NaD,TiO1,TiO2               ','  "     "')
prihdr['modlbt4f'] = ('  HdelA, HgamA, HdelF, HgamF  ','Worthey & Ottaviani (1997)')
prihdr['modlbt4g'] = ('  CaII8498, CaII8542, CaII8662','Diaz, Terlevich, &')
prihdr['modlbt4h'] = ('  MgI8807                     ','      Terlevich (1989)')
prihdr['modlbt4i'] = ('  H8_3889, H9_3835, H10_3798  ','Marcillac et al. (2006)')
prihdr['modlbt4j'] = ('  BHHK                        ','Brodie & Hanes HK index')
prihdr['modlbt4k'] = ('  D4000                       ','4000 A break index (Gorgas&Cardiel)')
prihdr['modlbt4l'] = ('  B4000VN                     ','Very narrow 4000 A break')
prihdr['modlbt4m'] = ('  BL1302,SiIV,BL1425,Fe1453,  ','UV, Fanelli et al. (1987,1990,1992)')
prihdr['modlbt4n'] = ('  CIV1548a,CIV1548c,CIV1548e, ',' "     "      "         "')
prihdr['modlbt4o'] = ('  BL1617,BL1664,BL1719,BL1853,',' "     "      "         "')
prihdr['modlbt4p'] = ('  FeII2402,BL2538,FeII2609,   ',' "     "      "         "')
prihdr['modlbt4q'] = ('  MgII,MgI,Mgwide,FeI,BL3096  ',' "     "      "         "')
prihdr['modlbt4r'] = ('  CIVa                        ','CIV  1540 A stellar absorption')
prihdr['modlbt4s'] = ('  HeIIe                       ','HeII 1640 A stellar emission')

# Build fits compatible header for BINTABLE2 (physical properties)
sechdr = hdulist[2].header
sechdr['TTYPE1']  = 'logage'             # label for column 1
sechdr['TTYPE2']  = 'logNLy'             # label for column 2
sechdr['TTYPE3']  = 'logNHeI'            # label for column 3
sechdr['TTYPE4']  = 'logNHeII'           # label for column 4
sechdr['TTYPE5']  = 'logNHeII-logNLy'    # label for column 5
sechdr['TTYPE6']  = 'B912'               # label for column 6
sechdr['TTYPE7']  = 'B4000VN'            # label for column 7
sechdr['TTYPE8']  = 'B4000SDSS'          # label for column 8
sechdr['TTYPE9']  = 'B4000'              # label for column 9

sechdr['TTYPE10'] = 'SNR'                # label for column 10
sechdr['TTYPE11'] = 'PISNR'              # label for column 11
sechdr['TTYPE12'] = 'RegSNR'             # label for column 12
sechdr['TTYPE13'] = 'TypeIaSNR'          # label for column 13
sechdr['TTYPE14'] = 'FailedSNR'          # label for column 14
sechdr['TTYPE15'] = 'PNBR'               # label for column 15
sechdr['TTYPE16'] = 'NBH'                # label for column 16
sechdr['TTYPE17'] = 'NNS'                # label for column 17
sechdr['TTYPE18'] = 'NWD'                # label for column 18
sechdr['TTYPE19'] = 'Lx'                 # label for column 19

sechdr['TTYPE20'] = 'Mstars'             # label for column 20
sechdr['TTYPE21'] = 'Mremnants'          # label for column 21
sechdr['TTYPE22'] = 'Mretgas'            # label for column 22
sechdr['TTYPE23'] = 'Mgalaxy'            # label for column 23
sechdr['TTYPE24'] = 'SFR'                # label for column 24
sechdr['TTYPE25'] = 'Mtot'               # label for column 25
sechdr['TTYPE26'] = 'Mtot_Lb'            # label for column 26
sechdr['TTYPE27'] = 'Mtot_Lv'            # label for column 27
sechdr['TTYPE28'] = 'Mtot_Lk'            # label for column 28
sechdr['TTYPE29'] = 'Mliv_Lb'            # label for column 29

sechdr['TTYPE30'] = 'Mliv_Lv'            # label for column 30
sechdr['TTYPE31'] = 'Mliv_Lk'            # label for column 31
sechdr['TTYPE32'] = 'EvolutionaryFlux'   # label for column 32
sechdr['TTYPE33'] = 'SpecificEvolFlux'   # label for column 33
sechdr['TTYPE34'] = 'TurnoffMass'        # label for column 34
sechdr['TTYPE35'] = 'BPMS_BMS'           # label for column 35
sechdr['TTYPE36'] = 'TotalMassLossRate'  # label for column 36
sechdr['TTYPE37'] = 'DustProductionRate' # label for column 37

sechdr['TUNIT1']  = 'yr            (log age of SSP)'
sechdr['TUNIT2']  = 'photons/Mo    (log number of HI ionizing photons)'
sechdr['TUNIT3']  = 'photons/Mo    (log number of HeI ionizing photons)'
sechdr['TUNIT4']  = 'photons/Mo    (log number of HeII ionizing photons)'
sechdr['TUNIT5']  = '..........    (log ratio number of ionizing HeII/HI)'
sechdr['TUNIT6']  = '..........    (Lyman break amplitude)'
sechdr['TUNIT7']  = '..........    (4000 A break, very narrow definition)'
sechdr['TUNIT8']  = '..........    (4000 A break, SDSS definition)'
sechdr['TUNIT9']  = '..........    (4000 A break, original definition)'

sechdr['TUNIT10'] = 'SN/yr         (Total Supernova rate)'
sechdr['TUNIT11'] = 'SN/yr         (PI Supernova rate)'
sechdr['TUNIT12'] = 'SN/yr         (Regular Supernova rate)'
sechdr['TUNIT13'] = 'SN/yr         (Type Ia Supernova rate)'
sechdr['TUNIT14'] = 'SN/yr         (Failed Supernova rate)'
sechdr['TUNIT15'] = 'PN/yr         (Planetary Nebula birth rate)'
sechdr['TUNIT16'] = 'number/Mo     (Number of Black holes)'
sechdr['TUNIT17'] = 'number/Mo     (Number of Neutron stars)'
sechdr['TUNIT18'] = 'number/Mo     (Number of White dwarfs)'
sechdr['TUNIT19'] = 'Lo            (X ray luminosity)'

sechdr['TUNIT20'] = 'Mo            (Mliv = Mass in living stars)'
sechdr['TUNIT21'] = 'Mo            (Mrem = Mass in stellar remnants)'
sechdr['TUNIT22'] = 'Mo            (Mgas = Mass of processed gas)'
sechdr['TUNIT23'] = 'Mo            (Mass of galaxy)'
sechdr['TUNIT24'] = 'Mo/yr         (Star formation rate)'
sechdr['TUNIT25'] = 'Mo            (Mtot = Mliv + Mrem)'
sechdr['TUNIT26'] = 'Mo/Lo         (Mtot/Light ratio B band)'
sechdr['TUNIT27'] = 'Mo/Lo         (Mtot/Light ratio V band)'
sechdr['TUNIT28'] = 'Mo/Lo         (Mtot/Light ratio K band)'
sechdr['TUNIT29'] = 'Mo/Lo         (Mliv/Light ratio B band)'

sechdr['TUNIT30'] = 'Mo/Lo         (Mliv/Light ratio V band)'
sechdr['TUNIT31'] = 'Mo/Lo         (Mliv/Light ratio K band)'
sechdr['TUNIT32'] = 'Nstars/yr     (Evolutionary flux)'
sechdr['TUNIT33'] = 'Nstars/yr/Lo  (Specific Evolutionary flux)'
sechdr['TUNIT34'] = 'Mo            (MS Turnoff mass)'
sechdr['TUNIT35'] = '..........    (PostMS/MS bolometric flux)'
sechdr['TUNIT36'] = 'M/Mo/yr       (Total Mass Loss Rate)'
sechdr['TUNIT37'] = 'M/Mo/yr       (Total Dust Production Rate)'

# Build fits compatible header for BINTABLE3 (photometry)
sechdr = hdulist[3].header
sechdr['TTYPE1']  = 'logage'
sechdr['TUNIT1']  = 'yr       (log age of SED)'        # label for column 1
sechdr['TUNIT2']  = 'Vega mag (Bolometric magnitude)'  # label for column 2
sechdr['TUNIT3']  = 'AB mag   (FUV_GALEX_AB)'          # label for column 3
sechdr['TUNIT4']  = 'AB mag   (ACSWFC_F220w_AB)'       # label for column 4
sechdr['TUNIT5']  = 'AB mag   (NUV_GALEX_AB)'          # label for column 5
sechdr['TUNIT6']  = 'AB mag   (WFC3_F225W_AB)'         # label for column 6
sechdr['TUNIT7']  = 'AB mag   (UVIS1_f225w_AB)'        # label for column 7
sechdr['TUNIT8']  = 'AB mag   (UVIS1_f275w_AB)'        # label for column 8
sechdr['TUNIT9']  = 'AB mag   (ACSWFC_F250w_AB)'       # label for column 9
sechdr['TUNIT10'] = 'AB mag   (ACSWFC_F330w_AB)'       # label for column 10
sechdr['TUNIT11'] = 'AB mag   (UVIS1_f336w_AB)'        # label for column 11
sechdr['TUNIT12'] = 'AB mag   (WFC3_F336W_AB)'         # label for column 12
sechdr['TUNIT13'] = 'AB mag   (u_SDSS_AB)'             # label for column 13
sechdr['TUNIT14'] = 'Vega mag (U_Johnson)'             # label for column 14
sechdr['TUNIT15'] = 'AB mag   (u3_CFHT_MC_AB)'         # label for column 15
sechdr['TUNIT16'] = 'AB mag   (u1_CFHT_MC_AB)'         # label for column 16
sechdr['TUNIT17'] = 'AB mag   (ACSWFC_F410w_AB)'       # label for column 17
sechdr['TUNIT18'] = 'AB mag   (WFC3_FR388N_AB)'        # label for column 18
sechdr['TUNIT19'] = 'AB mag   (ACSWFC_F435w_AB)'       # label for column 19
sechdr['TUNIT20'] = 'AB mag   (WFC3_F438W_AB)'         # label for column 20
sechdr['TUNIT21'] = 'AB mag   (UVIS1_f438w_AB)'        # label for column 21
sechdr['TUNIT22'] = 'Vega mag (B3_Johnson)'            # label for column 22
sechdr['TUNIT23'] = 'Vega mag (B2_Johnson)'            # label for column 23
sechdr['TUNIT24'] = 'AB mag   (g_SDSS_AB)'             # label for column 24
sechdr['TUNIT25'] = 'AB mag   (ACSWFC_F475w_AB)'       # label for column 25
sechdr['TUNIT26'] = 'AB mag   (g3_CFHT_MC_AB)'         # label for column 26
sechdr['TUNIT27'] = 'AB mag   (g1_CFHT_MC_AB)'         # label for column 27
sechdr['TUNIT28'] = 'AB mag   (UVIS1_f555w_AB)'        # label for column 28
sechdr['TUNIT29'] = 'AB mag   (WFC3_F555W_AB)'         # label for column 29
sechdr['TUNIT30'] = 'AB mag   (ACSWFC_F555w_AB)'       # label for column 30
sechdr['TUNIT31'] = 'AB mag   (UVIS1_f547m_AB)'        # label for column 31
sechdr['TUNIT32'] = 'Vega mag (V_Johnson)'             # label for column 32
sechdr['TUNIT33'] = 'AB mag   (ACSWFC_F606w_AB)'       # label for column 33
sechdr['TUNIT34'] = 'AB mag   (UVIS1_f606w_AB)'        # label for column 34
sechdr['TUNIT35'] = 'AB mag   (r_SDSS_AB)'             # label for column 35
sechdr['TUNIT36'] = 'AB mag   (r1_CFHT_MC_AB)'         # label for column 36
sechdr['TUNIT37'] = 'AB mag   (UVIS1_f625w_AB)'        # label for column 37
sechdr['TUNIT38'] = 'AB mag   (ACSWFC_F625w_AB)'       # label for column 38
sechdr['TUNIT39'] = 'AB mag   (r3_CFHT_MC_AB)'         # label for column 39
sechdr['TUNIT40'] = 'AB mag   (UVIS1_f656n_AB)'        # label for column 40
sechdr['TUNIT41'] = 'AB mag   (UVIS1_f657n_AB)'        # label for column 41
sechdr['TUNIT42'] = 'AB mag   (UVIS1_f658n_AB)'        # label for column 42
sechdr['TUNIT43'] = 'Vega mag (R_Cousins)'             # label for column 43
sechdr['TUNIT44'] = 'Vega mag (R_Johnson)'             # label for column 44
sechdr['TUNIT45'] = 'AB mag   (i_SDSS_AB)'             # label for column 45
sechdr['TUNIT46'] = 'AB mag   (i2_CFHT_MC_AB)'         # label for column 46
sechdr['TUNIT47'] = 'AB mag   (i3_CFHT_MC_AB)'         # label for column 47
sechdr['TUNIT48'] = 'AB mag   (ACSWFC_F775w_AB)'       # label for column 48
sechdr['TUNIT49'] = 'Vega mag (I_Cousins)'             # label for column 49
sechdr['TUNIT50'] = 'AB mag   (UVIS1_f814w_AB)'        # label for column 50
sechdr['TUNIT51'] = 'AB mag   (ACSWFC_F814w_AB)'       # label for column 51
sechdr['TUNIT52'] = 'AB mag   (WFC3_F814W_AB)'         # label for column 52
sechdr['TUNIT53'] = 'Vega mag (I_Johnson)'             # label for column 53
sechdr['TUNIT54'] = 'AB mag   (z1_CFHT_MC_AB)'         # label for column 54
sechdr['TUNIT55'] = 'AB mag   (z_SDSS_AB)'             # label for column 55
sechdr['TUNIT56'] = 'AB mag   (z3_CFHT_MC_AB)'         # label for column 56
sechdr['TUNIT57'] = 'AB mag   (WFC3_F110W_AB)'         # label for column 57
sechdr['TUNIT58'] = 'Vega mag (J_2Mass)'               # label for column 58
sechdr['TUNIT59'] = 'Vega mag (J_Johnson)'             # label for column 59
sechdr['TUNIT60'] = 'AB mag   (WFC3_F125W_AB)'         # label for column 60
sechdr['TUNIT61'] = 'Vega mag (J_Palomar)'             # label for column 61
sechdr['TUNIT62'] = 'AB mag   (WFC3_F160W_AB)'         # label for column 62
sechdr['TUNIT63'] = 'Vega mag (H_Palomar)'             # label for column 63
sechdr['TUNIT64'] = 'AB mag   (H_2MASS_AB)'            # label for column 64
sechdr['TUNIT65'] = 'Vega mag (H_2Massr)'              # label for column 65
sechdr['TUNIT66'] = 'Vega mag (Kprime_Cowie)'          # label for column 66
sechdr['TUNIT67'] = 'AB mag   (Ks_CFHT_WC_AB)'         # label for column 67
sechdr['TUNIT68'] = 'Vega mag (Ks_2Mass)'              # label for column 68
sechdr['TUNIT69'] = 'Vega mag (K_Johnson)'             # label for column 69
sechdr['TUNIT70'] = 'Vega mag (K_Palomar)'             # label for column 70
sechdr['TUNIT71'] = 'Vega mag (L_Johnson)'             # label for column 71
sechdr['TUNIT72'] = 'Vega mag (I3p6_IRAC)'             # label for column 72
sechdr['TUNIT73'] = 'Vega mag (I4p5_IRAC)'             # label for column 73
sechdr['TUNIT74'] = 'Vega mag (I5p7_IRAC)'             # label for column 74
sechdr['TUNIT75'] = 'Vega mag (I7p9_IRAC)'             # label for column 75
sechdr['TUNIT76'] = 'Vega mag (I12_IRAS)'              # label for column 76
sechdr['TUNIT77'] = 'Vega mag (M24_MIPS)'              # label for column 77
sechdr['TUNIT78'] = 'Vega mag (I25_IRAS)'              # label for column 78
sechdr['TUNIT79'] = 'Vega mag (I60_IRAS)'              # label for column 79
sechdr['TUNIT80'] = 'Vega mag (M70_MIPS)'              # label for column 80
sechdr['TUNIT81'] = 'Vega mag (I100_IRAS)'             # label for column 81
sechdr['TUNIT82'] = 'Vega mag (M160_MIPS)'             # label for column 82

# Build fits compatible header for BINTABLE4 (line indices)
sechdr = hdulist[4].header
sechdr['TTYPE1']  = 'logage'             # label for column 1
sechdr['TTYPE2']  = 'CN1'                # label for column 2
sechdr['TTYPE3']  = 'CN2'                # label for column 3
sechdr['TTYPE4']  = 'Ca4227'             # label for column 4
sechdr['TTYPE5']  = 'G4300'              # label for column 5
sechdr['TTYPE6']  = 'Fe4383'             # label for column 6
sechdr['TTYPE7']  = 'Ca4455'             # label for column 7
sechdr['TTYPE8']  = 'Fe4531'             # label for column 8
sechdr['TTYPE9']  = 'Fe4668'             # label for column 9
sechdr['TTYPE10'] = 'Hbeta'              # label for column 10
sechdr['TTYPE11'] = 'Fe5015'             # label for column 11
sechdr['TTYPE12'] = 'Mg1'                # label for column 12
sechdr['TTYPE13'] = 'Mg2'                # label for column 13
sechdr['TTYPE14'] = 'Mgb'                # label for column 14
sechdr['TTYPE15'] = 'Fe5270'             # label for column 15
sechdr['TTYPE16'] = 'Fe5335'             # label for column 16
sechdr['TTYPE17'] = 'Fe5406'             # label for column 17
sechdr['TTYPE18'] = 'Fe5709'             # label for column 18
sechdr['TTYPE19'] = 'Fe5782'             # label for column 19
sechdr['TTYPE20'] = 'NaD'                # label for column 20
sechdr['TTYPE21'] = 'TiO1'               # label for column 21
sechdr['TTYPE22'] = 'TiO2'               # label for column 22
sechdr['TTYPE23'] = 'HdeltaA'            # label for column 23
sechdr['TTYPE24'] = 'HgammaA'            # label for column 24
sechdr['TTYPE25'] = 'HdeltaF'            # label for column 25
sechdr['TTYPE26'] = 'HgammaF'            # label for column 26
sechdr['TTYPE27'] = 'D4000'              # label for column 27
sechdr['TTYPE28'] = 'B4000VN'            # label for column 28
sechdr['TTYPE29'] = 'CaII8498'           # label for column 29
sechdr['TTYPE30'] = 'CaII8542'           # label for column 30
sechdr['TTYPE31'] = 'CaII8662'           # label for column 31
sechdr['TTYPE32'] = 'MgI8807'            # label for column 32
sechdr['TTYPE33'] = 'H83889'             # label for column 33
sechdr['TTYPE34'] = 'H93835'             # label for column 34
sechdr['TTYPE35'] = 'H103798'            # label for column 35
sechdr['TTYPE36'] = 'BHHK'               # label for column 36
sechdr['TTYPE37'] = 'BL1302'             # label for column 37
sechdr['TTYPE38'] = 'SiIV'               # label for column 38
sechdr['TTYPE39'] = 'BL1425'             # label for column 39
sechdr['TTYPE40'] = 'Fe1453'             # label for column 40
sechdr['TTYPE41'] = 'CIV1548a'           # label for column 41
sechdr['TTYPE42'] = 'CIV1548c'           # label for column 42
sechdr['TTYPE43'] = 'CIV1548e'           # label for column 43
sechdr['TTYPE44'] = 'BL1617'             # label for column 44
sechdr['TTYPE45'] = 'BL1664'             # label for column 45
sechdr['TTYPE46'] = 'BL1719'             # label for column 46
sechdr['TTYPE47'] = 'BL1853'             # label for column 47
sechdr['TTYPE48'] = 'FeII2402'           # label for column 48
sechdr['TTYPE49'] = 'BL2538'             # label for column 49
sechdr['TTYPE50'] = 'FeII2609'           # label for column 50
sechdr['TTYPE51'] = 'MgII'               # label for column 51
sechdr['TTYPE52'] = 'MgI'                # label for column 52
sechdr['TTYPE53'] = 'Mgwide'             # label for column 53
sechdr['TTYPE54'] = 'FeI'                # label for column 54
sechdr['TTYPE55'] = 'BL3096'             # label for column 55
sechdr['TTYPE56'] = 'CIVabs'             # label for column 56
sechdr['TTYPE57'] = 'HeIIems'            # label for column 57
sechdr['TUNIT1']  = 'yr   (log age of SED)'
sechdr['TUNIT2']  = 'mag  (CN1 Lick index)'
sechdr['TUNIT3']  = 'mag  (CN2 Lick index)'
sechdr['TUNIT4']  = 'A    (Ca4227 Lick index)'
sechdr['TUNIT5']  = 'A    (G4300 Lick index)'
sechdr['TUNIT6']  = 'A    (Fe4383 Lick index)'
sechdr['TUNIT7']  = 'A    (Ca4455 Lick index)'
sechdr['TUNIT8']  = 'A    (Fe4531 Lick index)'
sechdr['TUNIT9']  = 'A    (Fe4668 Lick index)'
sechdr['TUNIT10'] = 'A    (Hbeta Lick index)'
sechdr['TUNIT11'] = 'A    (Fe5015 Lick index)'
sechdr['TUNIT12'] = 'mag  (Mg1 Lick index)'
sechdr['TUNIT13'] = 'mag  (Mg2 Lick index)'
sechdr['TUNIT14'] = 'A    (Mgb Lick index)'
sechdr['TUNIT15'] = 'A    (Fe5270 Lick index)'
sechdr['TUNIT16'] = 'A    (Fe5335 Lick index)'
sechdr['TUNIT17'] = 'A    (Fe5406 Lick index)'
sechdr['TUNIT18'] = 'A    (Fe5709 Lick index)'
sechdr['TUNIT19'] = 'A    (Fe5782 Lick index)'
sechdr['TUNIT20'] = 'A    (NaD Lick index)'
sechdr['TUNIT21'] = 'mag  (TiO1 Lick index)'
sechdr['TUNIT22'] = 'mag  (TiO2 Lick index)'
sechdr['TUNIT23'] = 'A    (Worthey & Ottaviani HdeltaA index)'
sechdr['TUNIT24'] = 'A    (Worthey & Ottaviani HgammaA index)'
sechdr['TUNIT25'] = 'A    (Worthey & Ottaviani HdeltaF index)'
sechdr['TUNIT26'] = 'A    (Worthey & Ottaviani HgammaF index)'
sechdr['TUNIT27'] = '     (Gorgas & Cardiel 4000 A break index)'
sechdr['TUNIT28'] = '     (Very narrow 4000 A break index)'
sechdr['TUNIT29'] = 'A    (Diaz & Terlevich**2 CaII8498 index)'
sechdr['TUNIT30'] = 'A    (Diaz & Terlevich**2 CaII8542 index)'
sechdr['TUNIT31'] = 'A    (Diaz & Terlevich**2 CaII8662 index)'
sechdr['TUNIT32'] = 'A    (Diaz & Terlevich**2 MgI8807 index)'
sechdr['TUNIT33'] = 'A    (Marcillac et al. H83889 index)'
sechdr['TUNIT34'] = 'A    (Marcillac et al. H93835 index)'
sechdr['TUNIT35'] = 'A    (Marcillac et al. H103798 index)'
sechdr['TUNIT36'] = 'A    (Brodie & Hanes HK index)'
sechdr['TUNIT37'] = 'A    (Fanelli et al. BL1302 index)'
sechdr['TUNIT38'] = 'A    (Fanelli et al. SiIV index)'
sechdr['TUNIT39'] = 'A    (Fanelli et al. BL1425 index)'
sechdr['TUNIT40'] = 'A    (Fanelli et al. Fe1453 index)'
sechdr['TUNIT41'] = 'A    (Fanelli et al. CIV1548a index)'
sechdr['TUNIT42'] = 'A    (Fanelli et al. CIV1548c index)'
sechdr['TUNIT43'] = 'A    (Fanelli et al. CIV1548e index)'
sechdr['TUNIT44'] = 'A    (Fanelli et al. BL1617 index)'
sechdr['TUNIT45'] = 'A    (Fanelli et al. BL1664 index)'
sechdr['TUNIT46'] = 'A    (Fanelli et al. BL1719 index)'
sechdr['TUNIT47'] = 'A    (Fanelli et al. BL1853 index)'
sechdr['TUNIT48'] = 'A    (Fanelli et al. FeII2402 index)'
sechdr['TUNIT49'] = 'A    (Fanelli et al. BL2538 index)'
sechdr['TUNIT50'] = 'A    (Fanelli et al. FeII2609 index)'
sechdr['TUNIT51'] = 'A    (Fanelli et al. MgII index)'
sechdr['TUNIT52'] = 'A    (Fanelli et al. MgI index)'
sechdr['TUNIT53'] = 'A    (Fanelli et al. Mgwide index)'
sechdr['TUNIT54'] = 'A    (Fanelli et al. FeI index)'
sechdr['TUNIT55'] = 'A    (Fanelli et al. BL3096 index)'
sechdr['TUNIT56'] = 'A    (CIV 1540 A stellar absorption line index)'
sechdr['TUNIT57'] = 'A    (HeII 1640 A stellar emission line index)'
