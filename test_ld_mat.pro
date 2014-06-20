x_pix = findgen(5)
y_pix = findgen(5)

make_2d, x_pix, y_pix, x_pix, y_pix

x_pix = reform(x_pix, n_elements(x_pix))
y_pix = reform(y_pix, n_elements(y_pix))


x_src = [0, 0]
y_src = [3, 2]

dist_circle,dist,5,x_src[0],y_src[0]
dist1 = reform(dist, n_elements(dist))

dist_circle,dist,5,x_src[1],y_src[1]
dist2 = reform(dist, n_elements(dist))

dist = [[dist1],[dist2]]



LSTDRV_MATRIX, X_PIX, Y_PIX, X_SRC, Y_SRC, DX, DY, dr 


print, 'Should be zeros', minmax(dist-dr)
