PRO LSTDRV_MATRIX, X_PIX, Y_PIX, X_SRC, Y_SRC, DX, DY,  DR
;+
; NAME:
;      LSTDRV_MATRIX
;
;
; PURPOSE:
;      Constructs the matrix 
;
;
; CATEGORY:
;     List driven photometry - using maximum likelihood solution
;
;
; CALLING SEQUENCE:
;
;
;
; INPUTS:
;     X_PIX: vector of x-values for "pixels"
;     Y_PIX: vector of y-values for "pixels"
;     X_SRC: vector of x-values for sources
;     Y_SRC: vector of y-values for sources
;
;
;
; OPTIONAL INPUTS:
;
;
;
; KEYWORD PARAMETERS:
;
;
;
; OUTPUTS:
;
;    DX: offset in x in matrix format for linear inversion
;    DY: offset in y  in matrix format for linear inversion
;
; OPTIONAL OUTPUTS:
;   DR:  
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
on_error, 2

n_src = n_elements(x_src)
n_pix = n_elements(x_pix)

one_src = replicate(1, n_src)
one_pix = replicate(1, n_pix)

; n.b. reforming in case positions or sources are sent down as matricies

dx = (reform(x_pix, n_pix)#one_src)-(one_pix#reform(x_src, n_src))
dy = (reform(y_pix, n_pix)#one_src)-(one_pix#reform(y_src, n_src))

IF n_params() GT 6 THEN dr = sqrt(dx^2+dy^2)


return
END
