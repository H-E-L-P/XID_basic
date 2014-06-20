PRO simple_sim_map, srcs, map, naxis1, naxis2, fwhm, paxis1=paxis1, noise=noise, perfect_map=perfect_map, prf_out=prf_out, sig_map=sig_map, seed=seed, pix_noise=pix_noise, prf_norm=prf_norm

;+
; Name:
;  SIMPLE_SIM_MAP
;
;
; PURPOSE:
;  Takes an input list of sources and places them into a map with a given psf and noise
;
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE:
;
;
;
; INPUTS:
;      SRCS: Structure with 
;            SRCS.x: x-coordinate 
;            SRCS.y: y-coordinate
;            SRCS.f: flux
;      NAXIS1: size of image in x  (Default 
;      NAXIS2: size of image in y  (Default 
;      FWHM: FWHM of PSF in pixels (Default
;      PAXIS1: Size of smoothing PSF kernel in pixels (Default 
;              NOISE: Noise value (or image), sigma in units of 
;                     error in the final point-source flux estimate. 
;                     Various possible input forms
;      SEED: Random number seed (Default -1)
;
; OPTIONAL INPUTS:
;
;
;
; KEYWORD PARAMETERS:
;      PIX_NOISE: 0, noise per point source flux, 1: noise per pixel
;      PRF_NORM: 0, flux per pixel, 1: flux per beam
;
;
; OUTPUTS:
;      MAP: NOISY MAP in units of flux per pixel
;      PERFECT_MAP: No noise map
;                   SIG_MAP: Noise in the form of a map in units 
;                            of flux per pixel (N.B. not per point
;                            source as the noise is input!)
;      PRF:  Point Response function in units of per pixel
;
; OPTIONAL OUTPUTS:
;
;
;
; COMMON BLOCKS:
;
;
;
; SIDE EFFECTS:
;
;
;
; RESTRICTIONS:
;    Slightly clunky inputs due to history in list-driven photometry routine 
;    should be streamlined e.g. to take a fits header as input and any arbitrary psf
;    also should pull the random source generation stuff out
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;
;-

   
;----------------------------------------------------------------------
; parameter checking & defaults
;----------------------------------------------------------------------
; source parameters
IF NOT keyword_set(n_src) THEN n_src = 2
; map parameters
IF NOT keyword_set(naxis1) THEN naxis1 = 20
IF NOT keyword_set(naxis2) THEN naxis2 = naxis1
n_pix = long(naxis1)*long(naxis2)
; PRF parameters
IF NOT keyword_set(fwhm) THEN fwhm = 2
IF NOT keyword_set(paxis1) THEN paxis1 = (fix(fwhm)+1)*2+1
IF NOT keyword_set(paxis2) THEN paxis2 = paxis1

IF NOT keyword_set(seed) THEN seed = -1


x_src = srcs.x
y_src = srcs.y
f_src = srcs.f

;----------------------------------------------------------------------
; map
;----------------------------------------------------------------------



map = fltarr(naxis1, naxis2)


;----------------------------------------------------------------------
; putting sources in "zero-footprint" map
;----------------------------------------------------------------------

; really naive way of doing it, which doesn't allow for two sources in
; same pixel
; map[round(x_src), round(y_src)] = f_src

; better way of doing it

simple_map_xy, f_src, round(x_src), round(y_src), map


;----------------------------------------------------------------------
; Constructing PRF
;----------------------------------------------------------------------

paxis1 = paxis1
paxis2 = paxis2
fwhm = fwhm
prf=psf_gaussian(npixel=[paxis1,paxis2],fwhm=fwhm, /norm)  
CASE prf_norm OF
 0: prf = prf/total(prf)
 1: prf = prf/max(prf)
ENDCASE 
 
prf2_norm = total(prf^2)

;----------------------------------------------------------------------
; Convolution of PRF with zero-footprint map
;----------------------------------------------------------------------
;stop
map = convol(map, prf, /edge_zero)

;----------------------------------------------------------------------
; adding noise to map
;----------------------------------------------------------------------

perfect_map = map
IF NOT keyword_set(noise) THEN noise = 1
; assume noise is constant for now
sz = size(noise)
CASE  sz[0] OF  
   0: sig_map = replicate(noise, naxis1, naxis2)
   1: BEGIN
      IF sz[1] EQ 1 THEN sig_map = replicate(noise[0], naxis1, naxis2) ELSE message, 'Dont know how to use vector'
   END 
   2: BEGIN
      IF sz[1] NE naxis1 OR sz[2] NE naxis2 THEN message, 'Noise image wrong size'
      sig_map = noise
   END 
   ELSE: message, 'Dont know how to use noise array of '+sz[0]+'dimensions'
END

; scaling from flux error antifipacted for the point source in the optimally
; averaged data to error in each pixel.

CASE pix_noise OF 
   0: sig_map = sig_map*sqrt(prf2_norm)
   1: sig_map = sig_map
   ELSE: message, 'noise model not recognised'
ENDCASE 

; generating realisation of the noise

noisy_map = reform(randomn(seed, n_pix), naxis1, naxis2)*sig_map



; adding noise to map
noisy_map = noisy_map+map

map = noisy_map

prf_out = prf

END
