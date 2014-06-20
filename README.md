Code for prototyping the simultaneous list driven source profile photometric code.

---------------------------------------

##Overview
The full software will probably have three stages.  Stage 1: simulation, Stage 2. generation of linear matrices, Stage 3. solving matrix equations.  This prototype does the first two stages in one. 
The IDL code creates toy model source lists and maps and the required matrices required for linear fitting software. The fortran code does the linear algebra.

---------------------------------------

##Code

IDL Code (Seb Oliver).
* Highest level: do_lstdrv_wrapper.pro runs a few examples
	     : compare.pro creates binary fits file to compare input/output
* Mid level: lstdrv_wrapper.pro runs Stage 1 & 2 above
* Low level: lstdrv_matrix.pro calculates the pixel source offsets using matrix operations.
* Obsolete: test_ld_mat.pro lsrdrv_example.pro old testing routines
* tmp/ : all code including library routines required to run IDL code

Fortran code (Martin Kunz)
* ConjGrad.f: conjugate Gradient method for solving matrix equation.


#Data:
--------------------------------------------------
Various examples have been run

* eg1: two sources no noise
* eg2: two sources no noise
* eg3: a 400X400 image with 100 sources
* eg4: a 400X400 image with 600 sources
* eg5: a 400X400 image with 100 sources and noise increasing from the centre
* eg6: a 200X200 image with 600 sources and noise increasing from the centre

The data files are:
* '#.data'     : simulated detector signal       (Alt. Input to ConjGrad code)
	    : recorded as a vector
* '#.matrix'    : source-signal transfer matrix   (Input to ConjGrad code)
	    : recorded as a sparse matrix
* '#.noisydata' : simulated noisy data            (Input to ConjGrad code)
	    : recorded as a vector
* '#.sigma'     : sqrt(variance) vector           (Input to ConjGrad code)
	    : recorded as a vector
* '#.src'        : simulated input source flux     (for c.f. ConjGrad results)
	    : recorded as a vector

'DATA = MATRIX # SOURCE'

The fits files are
- '#_map.fits'        : simulated map with WCS
- '#_noisymap.fits'   : simulated noisy map
- '#_sigmap.fits'     : sqrt(variance) map
- '#_src.fits'        : sources in binary fits catalogue with RA, Dec
- '#_out.fits'	  : sources and fit fluxes in bin fits catalogue (produced by compare.pro)


Also
- svd_out.g3 outputs from early Matlab runs using SVD
- svd_out.g4 


Documentation
==============

* Readme: this file
* MaxLikeListDrivePhot.* : Latex document describing method









Seb Oliver December 2008.