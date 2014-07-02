PRO simple_map_xy, data, x, y, map, nsamp, weight,  weight_map,  flags,  flags_map
;+
; NAME:
;    SIMPLE_MAP
;
;
; PURPOSE:
;    Takes a series of data samples with coordinates and maps them to
;    an image header
;
;
; CATEGORY:
;    
;
;
; CALLING SEQUENCE:
;    SIMPLE_MAP, data, x, y, map, [weight, weight_map]
;
;
; INPUTS:
;   data: data samples
;   x: coordinates of samples
;   y: coordinates of samples
;
;
;
; OPTIONAL INPUTS:
;  weight: weights for data sample
;  flags: flags (or other variable you want co-added)
;
; KEYWORD PARAMETERS:
;
;
;
; OUTPUTS:
;  map: map with geometry defined by header and data-samples co-added
;       at nearest pixel
;  nsamp: number of sampels contributing to a map pixel (output)
;
;
; OPTIONAL OUTPUTS:
;  weight_map: weight map on same system as map
;  flags_map: flags map co-added on same system as map
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
;
;
;
; PROCEDURE:
;  uses the histogram function for speed.
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;
;  Original version December 19th 2007.  Seb Oliver
;
;-


; extracting dimension of image from header

IF NOT keyword_set(naxis1) THEN BEGIN
   IF keyword_set(map) THEN BEGIN 
      sz =  size(map)
      naxis1 = sz[1]
      naxis2 = sz[2]
   ENDIF ELSE BEGIN
      message, 'naxis1 not defined'
   ENDELSE 
ENDIF 

IF NOT keyword_set(naxis2) THEN naxis2 = naxis1 


map = fltarr(naxis1, naxis2)
weight_map = fltarr(naxis1, naxis2)
flag_map = fltarr(naxis1, naxis2)

IF n_params() LE  6 THEN  BEGIN
   weight = replicate(1, n_elements(data))
   DO_weight = 0
ENDIF ELSE DO_weight = 1

IF n_params() LE 8 THEN  DO_flag = 0 ELSE DO_flag = 1

; weighting data samples

data_tmp = data*weight


; checking which samples within field

inmap = where(x GT -0.5 AND x LT naxis1-0.5 AND y GT -0.5 AND y LT naxis2-0.5, ngood)
IF ngood LE 0 then message, 'No good data'

; rounding to nearest pixel
ix = round(x[inmap])
iy = round(y[inmap])

; creating 1D index 
index = long(iy)*naxis1+long(ix)

; using histogram to map samples to image space

nsamp = histogram(index, reverse_indi=r, min=0, max=long(naxis1)*naxis2-1)

; only work with map pixels that have some data

good = where(nsamp GT 0, ngood)
IF ngood LE 0 THEN message, 'No good data'

; this loop appears to be necessary and is the shortest loop I can
; imagine, being only 1
 
FOR j=0L, n_elements(good)-1L DO BEGIN 
   i = good[j]
;   IF nsamp[i] EQ 10 THEN stop
; See IDL manual to see where next two lines comes from 
   IF R[i] NE R[i+1] THEN MAP[i] = total(data_tmp[inmap[ [R[R[I] : R[i+1]-1]]]] ) $
                                         ELSE print, 'This should not happen'
   IF keyword_set(DO_weight) THEN weight_MAP[i] = total(weight[inmap [R[R[I] : R[i+1]-1]]] )
   IF keyword_set(DO_flags) THEN flag_MAP[i] = total(flags[inmap [R[R[I] : R[i+1]-1]] ])
   
ENDFOR

;reformating arrays to 2D

map = reform(map, naxis1, naxis2)
nsamp = reform(nsamp, naxis1, naxis2)
IF keyword_set(DO_weight) THEN weight_map = reform(weight_MAP, naxis1, naxis2) 
IF keyword_set(DO_flags) THEN weight_map = reform(flags_MAP, naxis1, naxis2) 


END

