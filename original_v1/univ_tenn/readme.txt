This code was written for a gentoo gnu/linux dual-socket dual-core amd64 workstation 
with 15 gigs of ram.  The code should run just fine on a machine with a different number 
of processors, but it will be sluggish at best if there is not 15 gigs of ram.

---

The following code is provided:

ctis.c:       performs the reconstruction

circulant.c:  writes columns of H to disk

complex.c:    complex number library (used by ctis.c, circulant.c)
complex.h

random.c:     random number library (used by ctis.c)
random.h

param:        #include in ctis.c and circulant.c; includes problem instance sizes

---

These are the steps for implementing our code:

1) install fftw (go to the fftw site and do the install)

   a) include --enable-threads in the flags to the configure script

2) populate the param file; it is clearly self documented

3) store the columns of your system matrix in the following file:  save-columns-Hmatrix.h\!

   use the following format

   col x
   y z 
   ..

   where x is a column number, y is a row number, and z is a system matrix value

   zero system matrix values need not be provided

   provide information for the first column of each new wavelength block only

4) compile circulant.c: 

   gcc -O3 -o circulant circulant.c -lfftw3

5) run circulant

6) decide the value for #defines to be used in step 8 (these are optional)

   FTHREADS     4 <--- number of threads for fftw (default)

   THREADS      4 <--- number of threads for other things (default)

   COMPUTE_PT   0 read precomputed P^T from disk (default)
                1 compute P^T and store to disk

   REMOVE_NOISE 0 do not remove noise
                1 remove noise (default)
                
   REPORT       0 do not report conjugate gradient error at each step
                1 report conjugate gradient error at each step (default)

   FOCAL_PLANE  "focal_plane" <--- g (default)

   PICTURE      0 do not produce pgm of reconstructed image
                1 produce pgm of reconstructed image (default)

   PICTURE_OUT  "reconstructed_image" <--- reconstructed image (default)

   SEED         37 <--- random number generator seed (default)

   EPSILON      (1e-8) <--- cg error theshold (defualt)

   MU           .01 <--- regularization parameter (default)

   LIM          2 <--- conjugate gradient iteration threshold (defualt)

   RESTART      4 <--- # of vose-horton steps (default)

7) compile complex.c:

   gcc -c complex.c

8) compile ctis.c:  

   gcc -O3 -o ctis ctis.c -D REMOVE_NOISE=1 -D REPORT=1 -D 'FOCAL_PLANE="focal_plane"' -D PICTURE=1 -D 'PICTURE_OUT="reconstructed_image"' -D SEED=37 -D EPSILON=(1e-8) -D MU=.01 -D LIM=2 -D RESTART=4 -D COMPUTE_PT=1 -D FTHREADS=4 -D THREADS=4 complex.o  -lfftw3_threads -lfftw3 -lm -lpthread  

   #defines are optional (default values will be used if not set here)

9) run ctis

