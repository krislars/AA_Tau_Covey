def getMag(band,wavelength,flux):
    
    import numpy as np
    import scipy.interpolate as interp
    import scipy.integrate as integrate
    
    
    # Convert filter files (except 2MASS) from nm to microns
    # Optical filter curves from http://www.aip.de/en/research/facilities/stella/instruments/data
    # SDSS centers/widths/F0 from http://www.astronomy.ohio-state.edu/~martini/usefuldata.html
    # NOTE: ugriz filters are on the AB magnitude system, UBVRIJHKs are on the Vega system
    
    if band=='Ks':
        bandwav,bandpass=np.loadtxt('utils/filters/Ks_2MASS.txt',unpack=True)
        center=2.159        # micron
        F0=4.283E-10        # W m^-2 micron^-1
        dlambda=0.262

    elif band=='H':
        bandwav,bandpass=np.loadtxt('utils/filters/H_2MASS.txt',unpack=True)
        center=1.662
        F0=1.133E-9
        dlambda=0.251

    elif band=='J':
        bandwav,bandpass=np.loadtxt('utils/filters/J_2MASS.txt',unpack=True)
        center=1.235
        F0=3.129E-9
        dlambda=0.162

    elif band=='U':
        bandwav_A,bandpass_100=np.loadtxt('utils/filters/Bessel_U.txt',unpack=True)
        bandwav=bandwav_A/1.0E3
        bandpass=bandpass_100/100.0
        center=0.365
        F0=4.18E-8
        dlambda=0.066

    elif band=='B':
        bandwav_A,bandpass_100=np.loadtxt('utils/filters/Bessel_B.txt',unpack=True)
        bandwav=bandwav_A/1.0E3
        bandpass=bandpass_100/100.0
        center=0.445
        F0=6.32E-8
        dlambda=0.094

    elif band=='V':
        bandwav_A,bandpass_100=np.loadtxt('utils/filters/Bessel_V.txt',unpack=True)
        bandwav=bandwav_A/1.0E3
        bandpass=bandpass_100/100.0
        center=0.551
        F0=3.63E-8
        dlambda=0.088

    elif band=='R':
        bandwav_A,bandpass_100=np.loadtxt('utils/filters/Bessel_R.txt',unpack=True)
        bandwav=bandwav_A/1.0E3
        bandpass=bandpass_100/100.0
        center=0.658
        F0=2.18E-8
        dlambda=0.138

    elif band=='I':
        bandwav_A,bandpass_100=np.loadtxt('utils/filters/Bessel_I.txt',unpack=True)
        bandwav=bandwav_A/1.0E3
        bandpass=bandpass_100/100.0
        center=0.806
        F0=1.13E-8
        dlambda=0.149
                                      
        

    filterband=np.zeros(wavelength.size)
        
    bandinterp=interp.interp1d(bandwav,bandpass)
            # 1D function between x=bandwav and y=bandpass, y=f(x)
        
    if bandwav[0]<bandwav[-1]: # increasing wavelength
        inband=np.logical_and(wavelength>bandwav[0],wavelength<bandwav[-1] )
    else: # decreasing wavelength
        inband=np.logical_and(wavelength<bandwav[0],wavelength>bandwav[-1] )
             #Array of Trues whereever wavelength is in range of bandwav... basically just an index of bandwav in wavelength array
        
    filterband[inband]=bandinterp(wavelength[inband])
            # Now filterband is the same shape as wavelength and flux, zeros and bandpass
    
    
    mag=2.5*np.log10(F0*dlambda/integrate.trapz(flux*filterband, x=wavelength))


    return mag