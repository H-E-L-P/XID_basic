; synthetic sources

x_src = [10, 12]
y_src = [10, 10]
f_src = [1, 2]
;x_src = [10]
;y_src = [10]
;f_src = [1]
n_src = n_elements(x_src)


  src = {x:0., y:0., f:0.0, ra:0.0d0, dec:0.0d0}
   srcs = replicate(src, n_src)
   srcs.x = x_src
   srcs.y = y_src
   srcs.f = f_src

; map
naxis1 = 20
naxis2 = 20
make_2d, findgen(naxis1), findgen(naxis2), x_pix, y_pix
map = fltarr(naxis1, naxis2)
n_pix = naxis1*naxis2

; sources in map

map[x_src, y_src] = f_src

; PRF
paxis1 = 5
paxis2 = 5
fwhm = 2
prf=psf_gaussian(npixel=[paxis1,paxis2],fwhm=fwhm)

; Convolution

map = convol(map, prf)

; check maximum pixel is where we expect it to be
tmp = max(map, ij)
whereimage, size(map), ij, i, j
print, i, j


; matrix of pixel <-->src positions

LSTDRV_MATRIX, X_PIX, Y_PIX, X_SRC, Y_SRC, DX, DY,  DR

; getting PRF values at dx, dy offsets

; dx and dy are defined relative to the centre of the PRF image 
; we need them in IDL pixel coordinates
px = dx+paxis1/2
py = dy+paxis2/2
good = where(px GE 0 AND px LT paxis1 AND py GE 0 AND py LT paxis2, ngood, comp=bad, ncomp=nbad)
PRF_matrix = fltarr(n_pix, n_src)
IF ngood GT 0 THEN PRF_matrix[good] = prf[px[good], py[good]]
IF nbad GT 0 THEN PRF_matrix[bad] = 0


; check my making up an image using the matrix approach

map2 = PRF_matrix#f_src
map2 = reform(map2, naxis1, naxis2)

tmp = max(map2, ij)
 whereimage, size(map2), ij, i, j   


good = where(map GT 0)
good2 = where(map2 GT 0)

; quick check
icplot, map-map2


;----------------------------------------------------------------------
; writing out data
;----------------------------------------------------------------------

name = 'eg1'

get_lun, unit
openw, unit, 'eg1.data'
printf, unit, '# '+string(n_pix)+' elelement DATA vector'
data = reform(map, n_pix)
printf, unit, data
free_lun, unit

get_lun, unit
openw, unit, 'eg1.matrix'
printf, unit, '# '+string(n_pix)+' X'+string(n_src)+' elelement N_DATA X N_SRC response matrix. Formatted as a sparse matrix i,j, M(i,j) with 0=<i<N_DATA, 0=<j<N_SRC and only if M(i,j) non zero'
make_2d, indgen(n_pix), indgen(n_src), i, j
good = where(PRF_matrix NE 0, ngood)
FOR k=0, ngood-1 DO printf, unit, i[good[k]], j[good[k]], PRF_matrix[good[k]]
free_lun, unit
