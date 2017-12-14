def getMag(band,wavelength,flux):
    
    import numpy as np
    import scipy.interpolate as interp
    import scipy.integrate as integrate
    
    from astropy.io import fits
    from astropy import units as u
    
    optical=['U','B','V','R','I']
    infrared=['J','H','Ks']
    
    if band in optical:
        x,y=np.loadtxt('utils/BessellMurphy2012/BessellMurphy_'
                       +band+'.txt',unpack=True)
        bandwav=x*u.AA
        bandpass=y
    elif band in infrared:
        x,y=np.loadtxt('utils/2MASS/'+band+'_2MASS.txt',unpack=True)
        bandpass=y/x # see explanatory supplement
        bandwav=x*u.micron
    
    hdu = fits.open('utils/BessellMurphy2012/alpha_lyr_stis_008.fits')
    rawwavelength=hdu[1].data['WAVELENGTH']
    rawflux=hdu[1].data['FLUX']
    vegaflux=rawflux * u.erg / (u.cm * u.cm * u.s * u.AA)
    vegawavelength=rawwavelength * u.AA
    # note: interpolation doesn't work in mixed units, so must convert    
    vegawavelength=vegawavelength.to(wavelength.unit)
    
    
    # function to compute band at arbitrary wavelength
    # note: interpolation doesn't work in mixed units, so must convert
    bandinterp=interp.interp1d(bandwav.to(wavelength.unit),bandpass)
    
    # has the size of wavelengths, TRUE where wavelength in the filter bandpass
    inband=np.logical_and(wavelength>np.min(bandwav),wavelength<np.max(bandwav))
    vegainband=np.logical_and(vegawavelength>np.min(bandwav),
                              vegawavelength<np.max(bandwav))
        
    sourcesum=integrate.trapz(flux[inband]*wavelength[inband]*
                              bandinterp(wavelength[inband]),x=wavelength[inband])
    vegasum=integrate.trapz(vegaflux[vegainband]*vegawavelength[vegainband]*
                              bandinterp(vegawavelength[vegainband]),
                              x=vegawavelength[vegainband])
    
    return -2.5*np.log10(np.abs(sourcesum)/vegasum)