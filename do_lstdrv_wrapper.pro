; synthetic sources

x_src = [10, 12]
y_src = [10, 10]
f_src = [1, 2]
n_src = n_elements(x_src)


  src = {x:0., y:0., f:0.0, ra:0.0d0, dec:0.0d0}
   srcs = replicate(src, n_src)
   srcs.x = x_src
   srcs.y = y_src
   srcs.f = f_src
lstdrv_wrapper,20,20,2,2,5,'eg2', srcs=srcs


lstdrv_wrapper,400,400,100,3,9,'eg3'                  


lstdrv_wrapper,400,400,600,3,9,'eg4'                  

naxis1 = 400
naxis2 = 400
 dist_circle,sig,[naxis1,naxis2]      
sig = sig/50.+1
lstdrv_wrapper,naxis1, naxis2,100,3,9,'eg5',noise=sig

naxis1 = 200
naxis2 = 200
 dist_circle,sig,[naxis1,naxis2]      
sig = (sig/50.+1)*5
Lstdrv_wrapper,naxis1, naxis2,600,3,9,'eg6',noise=sig
