## Column Description

- sp\_cat: The name of the spectroscopic catalogue where a match was found, i.e. deep2-dr4, sdss-dr10, vipers, and vvds.

- sp\_cat\_id: The object ID in the corrsponding spectroscopic catalogue.

- true\_class: The spectroscopic label. 0 for a galaxy, 1 for a star.

- z\_spec: Spectroscopic redshift.

### CFHTLenS Column Description

For details, see
[CFTHTLenS CATALOGUE DATA DESCRIPTION](http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/community/CFHTLens/README_catalogs_release.txt)

- field: Field identifier composed of the region (first two characters), the relative position in R.A. with respect to the centre of this region (3rd&4th character; "m" indicates minus, i.e. lower R.A. than the centre and "p" indicates plus, i.e. larger R.A. than the centre), and the relative position in Dec with respect to the centre of this region (3rd&4th character; "m" indicates minus, i.e. lower Dec than the centre and "p" indicates plus, i.e. larger Dec than the centre), 

- MASK: If MASK>0 the object lies within a mask. However, objects with MASK<=1 can safely be used for most scientific purposes.

- T\_B: BPZ spectral type. 1=CWW-Ell, 2=CWW-Sbc, 3=CWW-Scd, 4=CWW-Im, 5=KIN-SB3, 6=KIN-SB2. Note that we use a recalibrated template set described in Capak et al. (2004) and that the templates are interpolated, hence fractional types occur.

- NBPZ\_FILT, NBPZ\_FLAGFILT, NBPZ\_NONDETFILT: The number of filters in which an object has 'reliable photometry' (NBPZ\_FILT, magnitude errors < 1 mag and object brighter than the limiting magnitude); number of filters in which an object has formal magnitude errors of 1 mag or larger (NBPZ\_FLAGFILT); number of filters in which an object is fainter than the formal limiting magnitude (NBPZ\_NONDETFILT). If an object would fall into FLAGFILT AND NONDETFILT it is listed under FLAGFILT! Magnitude errors refer to MAG\_ISO.

- BPZ\_FILT, BPZ\_FLAGFILT, BPZ\_NONDETFILT: These keys contain a 'binary' encoding which filters fulfil the respective conditions. Filter u* is assigned a '1', "g'" = '2', "r'" = '4', "i'" = '8' and "z'" = '16'.  The keys BPZ\_FILT etc.  represent the sums of the filters fulfilling the corresponding criteria.

- PZ\_full: This is the full P(z) to z=3.5.  There are 70 columns sampling P(z) at intervals of dz=0.05.  The ﬁrst bin is centered at z=0.025. Note these 70 columns do not always sum to 1.  There is a ﬁnal bin not included in the catalogues with z>3.5 that, in a small number of cases, has non-zero probability. In CFHTLenS analysis we set a hard prior that there is zero probability past z>3.5, which corresponds to normalising each P(z) to one. For future ﬂexibility however we do not impose this normalisation on the catalogue, leaving it to the user to apply.

- star\_flag: Stars and galaxies are separated using a combination of size, i/y-band magnitude and colour information.  For i<21, all objects with size smaller than the PSF are classified as stars.  For i>23, all objects are classified as galaxies.  In the range 21<i<23, a star is defined as size<PSF and chi2\_star<2.0*chi2\_gal, where the chi2's are the best fit chi2's from the galaxy and star libraries given by LePhare.  NOTE: star\_flag is optimized for galaxy studies, to keep an almost 100% complete galaxy sample with low (but not vanishing) stellar contamination.  CLASS\_STAR usually gives a cleaner star sample, but can lead to serious incompleteness in a galaxy sample.

- MAG\_LIM\_[ugriyz]: These are 1-sigma limiting magnitude measured in a circular aperture with a diameter of 2 * FWHM, where FWHM is the seeing in this band (see image header).

- FIELD\_POS: Column needed for compliance with the FITS\_LDAC standard.

- WEIGHT: The lensfit weight to be used in shear measurement for each galaxy as given by equation 8 of Miller et al.  (2012).  The weight is effectively an inverse-variance weight, including the variance of the intrinsic galaxy ellipticity distribution, and thus can be used to estimate the measurement uncertainty in the shear for an individual galaxy.

- FITCLASS: Object classification as returned by lensfit.  Classification values are 0: galaxy 1: star -1: bad fit: no useable data -2: bad fit: blended object -3: bad fit: miscellaneous reason -4: bad fit: chi-squared exceeds a critical value

- SCALELENGTH: lensfit galaxy model scalelength, marginalised over other parameters, as defined by Miller et al (2012).  To be used for the calculation of the multiplicative correction as in Miller et al (2012), equations 15-17.

- MODEL-FLUX: lensfit galaxy model total flux, in calibrated CCD data units

- SNRATIO: signal-to-noise ratio of the object, measured on a stack of the supplied exposures within a limiting isophote 2$\sigma$ above the noise.  To be used for the calculation of the multiplicative correction as in Miller et al (2012), equations 15-17.

- PSF-e1, PSF-e2: means of the PSF model ellipticity components measured on each exposure.  PSF ellipticities are derived from the PSF model at the location of each galaxy and are top-hat weighted with radius 8\,pixels

- PSF-Strehl-ratio: means of the PSF model "pseudo-Strehl ratio" values measured on each exposure.  The psuedo-Strehl ratio is defined as the fraction of light falling in the central pixel, and is inversely related to the PSF size.

- e1, e2: expectation values of galaxy ellipticity, from the marginalised galaxy ellipticity likelihood surface, to be used for shear measurement.  Ellipticity is defined as e=(a-b)/(a+b) where a,b are the semi-major and semi-minor axes.  The ellipticity coordinate system is aligned with the celestial sphere at the location of each galaxy.  Corrections should be applied to these values as described by Miller et al (2012) and Heymans et al (2012).  We strongly urge the user not to use these raw uncalibrated ellipticity values blindly in a lensing analysis.  First, any shear measurement must measure weighted averages using the lensfit weight.  Secondly an additive c2 correction must be applied to the e2 component.  This can be calculated from equation (19) of Heymans et al 2012 from the quatities scalelength and SNratio noting that equation (19) is given in physical units (arcsec) whereas scalelength is given in pixel units.  One MegaCam CCD pixel is 0.187 arcsec.  Thirdly a multiplicative shear calibration correction must be applied following equations 15-17 of Miller et al (2012).  Note that it is incorrect to apply this multiplicative correction on an object by object basis.  Instead this calibration correction must be applied as as an ensemble average (see section 4.1 of Heymans et al 2012 for a summary of the required calibration corrections).

- N-EXPOSURES-USED: the number of individual exposures used by lensfit for this galaxy.

- PSF-e<1,2>-exp<i>: the lensfit PSF model ellipticity components 1,2 (top-hat weighted as above) on each exposure i.  
