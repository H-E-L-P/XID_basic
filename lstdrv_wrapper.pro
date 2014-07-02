PRO lstdrv_wrapper, naxis1, naxis2, n_src, fwhm, paxis1, name, srcs=srcs, noise=noise, debug=debug

;----------------------------------------------------------------------
; Source structure
;----------------------------------------------------------------------
IF NOT keyword_set(srcs) THEN BEGIN 
   src = {x:0., y:0., f:0.0, ra:0.0d0, dec:0.0d0}
   srcs = replicate(src, n_src)
   
;----------------------------------------------------------------------
; synthetic sources generated with Poisson 
;----------------------------------------------------------------------
   seed = -1
   srcs.x = round(randomu(seed, n_src)*naxis1-0.5)
   srcs.y = round(randomu(seed, n_src)*naxis2-0.5)
   srcs.f = 10.^randomn(seed,  n_src)
ENDIF 


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

IF NOT keyword_set(name) THEN name = 'eg2'

;----------------------------------------------------------------------
; WCS header - not relevant for the exercise but just means we can
; plot things in astro software only used at setup and output
;----------------------------------------------------------------------

IF NOT keyword_set(cdelt) THEN cdelt = 4./3600. ; 4 arcsec pixels like MIPS
IF NOT keyword_set(ra) THEN ra = 53.1267
IF NOT keyword_set(dec) THEN dec =  -27.805
IF NOT keyword_set(roll) THEN roll = 0
icmkhdr, naxis1, naxis2, cdelt, cdelt,ra, dec, roll, hd
;pdh21: changed prf keyowrd to prf_out to match, added prf_norm as it required variable to be set, added pix_noise as it required variable to be set
simple_sim_map, srcs, noisy_map, naxis1, naxis2, fwhm, paxis1=paxis1, noise=noise, prf_out=prf, perfect=map, sig_map=sig_map, seed=seed,prf_norm=0,pix_noise=0

;----------------------------------------------------------------------
; stripped background
;----------------------------------------------------------------------
make_2d, findgen(naxis1), findgen(naxis2), x_pix, y_pix
back = sin(x_pix/4)
back2 = rot(back, 20)*10

back3 = rot(back, 60)*10
back5 = randomn(seed, naxis1)#replicate(1, naxis2)
back5 = rot(back5, 60)*10

;----------------------------------------------------------------------
; matrix mapping pixels <--> src positions
;----------------------------------------------------------------------
x_src = srcs.x
y_src = srcs.y
make_2d, findgen(naxis1), findgen(naxis2), x_pix, y_pix

Lstdrv_MATRIX, X_PIX, Y_PIX, X_SRC, Y_SRC, DX, DY; ,  DR

;----------------------------------------------------------------------
; getting PRF values at dx, dy offsets
;----------------------------------------------------------------------

; dx and dy are defined relative to the centre of the PRF image 
; we need them in IDL pixel coordinates
dx = temporary(dx)+paxis1/2
dy = temporary(dy)+paxis2/2
good = where(dx GE 0 AND dx LT paxis1 AND dy GE 0 AND dy LT paxis2, ngood, comp=bad, ncomp=nbad)
PRF_matrix = fltarr(n_pix, n_src)
IF ngood GT 0 THEN PRF_matrix[good] = prf[round(dx[good]), round(dy[good])]
IF nbad GT 0 THEN PRF_matrix[bad] = 0
; setting bad to zero to save memory
bad = 0



;----------------------------------------------------------------------
; building the map the other way
;----------------------------------------------------------------------
f_src = srcs.f
map2 = PRF_matrix#f_src
map2 = reform(map2, naxis1, naxis2)


;----------------------------------------------------------------------
; checks
;----------------------------------------------------------------------

; check maximum pixel is where we expect it to be

tmp = max(map, ij)
whereimage, size(map), ij, i, j
tmpmax = max(f_src, imax)
i2 = round(x_src[imax])
j2 = round(y_src[imax])
IF i-i2 NE 0 OR j-j2 NE 0 THEN message, 'Missmatch in peak', /inf




; check my making up an image using the matrix approach


tmp = max(map2, ij)
 whereimage, size(map2), ij, i, j   

good = where(map GT 0)
good2 = where(map2 GT 0)

; quick check
delt = map-map2
check = where(delt NE 0., nbad)

IF nbad GT 0 THEN BEGIN 
   dev = delt[check]/map[check]
   maxdev = max(dev)

   IF maxdev GT 1.e-3 THEN message, 'Maps do not agree'
ENDIF


IF keyword_set(debug) THEN BEGIN 
   window, 1
   icplot, map
   window, 2
   icplot, map
   window, 3
   icplot, map-map2
    print,map,format='(20f5.2)' 
    print,map2,format='(20f5.2)' 
    print, (map-map2)*1000,format='(20f5.2)'  
ENDIF


;----------------------------------------------------------------------
; writing out data in ascii for matrix based stuff
;----------------------------------------------------------------------


get_lun, unit
openw, unit, name+'.src'
printf, unit, '# '+string(n_src)+' elelement SOURCE FLUX vector'
printf, unit, f_src
free_lun, unit


get_lun, unit
openw, unit, name+'.data'
printf, unit, '# '+string(n_pix)+' elelement DATA vector'
data = reform(map, n_pix)
printf, unit, data
free_lun, unit

get_lun, unit
openw, unit, name+'.noisydata'
printf, unit, '# '+string(n_pix)+' elelement DATA+NOISE vector'
data = reform(noisy_map, n_pix)
printf, unit, data
free_lun, unit

get_lun, unit
openw, unit, name+'.sigma'
printf, unit, '# '+string(n_pix)+' elelement SIGMA vector'
data = reform(sig_map, n_pix)
printf, unit, data
free_lun, unit


get_lun, unit
openw, unit, name+'.matrix'
printf, unit, '# '+string(n_pix)+' X'+string(n_src)+' elelement N_DATA X N_SRC response matrix. Formatted as a sparse matrix i,j, M(i,j) with 0=<i<N_DATA, 0=<j<N_SRC and only if M(i,j) non zero'
; note resuing px and py variables to save memory
make_2d, lindgen(n_pix), lindgen(n_src), px, py
good = where(PRF_matrix NE 0, ngood)
FOR k=0L, ngood-1L DO printf, unit, px[good[k]], py[good[k]], PRF_matrix[good[k]]
free_lun, unit
;files now also saved as idl save files
save,map,noisy_map,n_pix,sig_map,filename=name+'.sav'

;----------------------------------------------------------------------
; writing out catalogue and map data in fits for astronomy routines
;----------------------------------------------------------------------
extast, hd, ast
xy2ad, srcs.x, srcs.y, ast, ra, dec
srcs.ra = ra
srcs.dec = dec
mwrfits, srcs, name+'_src.fits', /create

writefits, name+'_map.fits', map, hd
writefits, name+'_sigmap.fits', sig_map, hd
writefits, name+'_noisymap.fits', noisy_map, hd
writefits, name+'_noisymap_b2.fits', noisy_map+back2, hd
writefits, name+'_noisymap_b23.fits', noisy_map+back2+back3, hd
writefits, name+'_noisymap_b5.fits', noisy_map+back5, hd


End
