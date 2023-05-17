# NIHAO-SKIRT-Catalog

The following data product has free model parameters of $f_{\rm dust}=0.1$ and $\tau_{\rm clear}=2.5$ Myrs.

### Link to [data product](https://drive.google.com/drive/folders/1eDouinsNBrEwPaQjvM7gVHa8kFe1rNGs?usp=share_link)

The two analysis files in this repository show how to load the and interpret the data. A brief summary is given below. 

Spatially integrated data: nihao-integrated-seds.fits 
  * ext='GALAXIES': one entry per galaxy (65 total) 
    * ['name'] = Name of the NIHAO galaxy 
    * ['stellar_mass'] = Stellar mass in $M_{\odot}$ 
    * ['sfr'] = Star formation rate in $M_{\odot}\ /\ \rm year$ averaged over 100 Myrs
    * ['sfh'] = Star formation history in $M_{\odot}$ as a function of age (['ages']) 
    * ['ceh'] = Chemical evolution hisotry in $M_{\odot}$ as a function of age (['ages']) and metallicity (['metals']) 
    * ['ages'] = Age bins in years corresponding to ['sfh'] and ['ceh'] 
    * ['metals'] = Metallicity bins corresponding to ['ceh'] 
    * ['dust_mass'] = Dust mass in $M_{\odot}$ 
    * ['axis_ratios'] = Axis ratios (b/a) 
    * ['size'] = Side length of modelled the cube in pcs 
  * ext='SUMMARY': one per orientation (650 total) 
    * ['name'] = Name of the NIHAO galaxy 
    * ['stellar_mass'] = Stellar mass in $M_{\odot}$ 
    * ['sfr'] = Star formation rate in $M_{\odot}\ /\ \rm year$ averaged over 100 Myrs
    * ['dust_mass'] = Dust mass in $M_{\odot}$ 
    * ['axis_ratio'] = Axis ratios (b/a) 
    * ['Av'] = Attenuation in the v-band in magnitudes 
    * ['attenuated_energy'] = Energy attenuated by dust in $10^{-23}\ erg\ s^{-1}\ cm^{-2}$ 
    * ['emitted_energy'] = Energy emitted by dust in $10^{-23}\ erg\ s^{-1}\ cm^{-2}$ 
    * ['bands'] = Name of broadband photometric bandpass  
    * ['flux'] = Flux of broadband photometric bandpass ['bands'] in Jansky 
    * ['flux_nodust'] = Dust-free flux of broadband photometric bandpass ['bands'] in Jansky 
    * ['size'] = Side length of modelled the cube in pcs 
  * ext='WAVE': Wavelengths of 'SPEC' and 'SPECNODUST' in Å (single array) 
  * ext='SPEC': Spectra in Jansky with wavelengths given by 'WAVE' (one array per orientation; 650 total) 
  * ext='SPECNODUST': Dust-free spectra in Jansky with wavelengths given by 'WAVE' (one array per orientation; 650 total) 
  * ext='ATTENUATION_WAVE': Wavelengths of 'ATTENUATION_MAGS' in Å (single array) 
  * ext='ATTENUATION_MAGS': Attenuation curves in magnitudes with wavelengths given by 'ATTENUATION_WAVE' 
    
Spatially resolved data: ResolvedData/[name]_nihao-resolved-photometry.fits (65 files)
  * ext='GALAXIES': one entry per file
    * ['name'] = Name of the NIHAO galaxy 
    * ['stellar_mass'] = Stellar mass in $M_{\odot}$ 
    * ['sfr'] = Star formation rate in $M_{\odot}\ /\ \rm year$ averaged over 100 Myrs
    * ['sfh'] = Star formation history in $M_{\odot}$ as a function of age (['ages']) 
    * ['ceh'] = Chemical evolution hisotry in $M_{\odot}$ as a function of age (['ages']) and metallicity (['metals']) 
    * ['ages'] = Age bins in years corresponding to ['sfh'] and ['ceh'] 
    * ['metals'] = Metallicity bins corresponding to ['ceh'] 
    * ['dust_mass'] = Dust mass in $M_{\odot}$ 
    * ['axis_ratios'] = Axis ratios (b/a) 
    * ['size'] = Side length of modelled the cube in pcs 
  * ext = 'SUMMARY': one per orientation (10 for each file) 
    * ['name'] = Name of the NIHAO galaxy 
    * ['stellar_mass'] = Stellar mass in $M_{\odot}$ 
    * ['sfr'] = Star formation rate in $M_{\odot}\ /\ \rm year$ averaged over 100 Myrs
    * ['dust_mass'] = Dust mass in $M_{\odot}$ 
    * ['axis_ratio'] = Axis ratios (b/a) 
    * ['Av'] = Attenuation in the v-band in magnitudes 
    * ['bands'] = Name of broadband photometric bandpass  
    * ['flux'] = Flux of broadband photometric bandpass ['bands'] in Jansky 
    * ['flux_nodust'] = Dust-free flux of broadband photometric bandpass ['bands'] in Jansky 
    * ['size'] = Side length of modelled the cube in pcs 
    
    
