data = mrdfits('eg3_src.fits', 1)
readcol, 'svd_out.g3', out
plot, data.f, out, psym=1

add_one_tag, data, 'fout', out, data2

mwrfits, data2, 'eg3_out.fits',/create


data = mrdfits('eg4_src.fits', 1)
readcol, 'svd_out.g4', out
plot, data.f, out, psym=1

add_one_tag, data, 'fout', out, data2

mwrfits, data2, 'eg4_out.fits',/create


