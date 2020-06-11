# ctis
ctis

## setup (skip if you have already done this)
 - install pthread: $ sudo apt-get install libpthread-stubs0-dev
 - extract the tar.gz: $ tar xvzf fftw-3.3.8.tar.gz
 - and cd into the directory $ cd fftw-3.3.8
 - $ ./configure --enable-threads
 - $ make
 - $ sudo make install

## Run

### Circulant System Matrix (skip if you don't want to re-generate A_0 or wisdom_327680)

For convenience we saved our output in *A_0*, but if you want to re-generate it you can do the following,

- given that the file "save-columns-Hmatrix.h!" is our system matrix for one wavelength.

$ gcc -o circulant circulant.c -lfftw3 -lm
$ ./circulant

- I changed the defines in *circulant.c* from
```
#define GX_DIM    2048
#define GY_DIM    2048     /*  GX_DIM*GY_DIM = # rows of H ; WVDIM*XDIM*YDIM = # columns of H */
#define WVDIM (int)    29
#define XDIM  (int)    81
#define YDIM  (int)    90
#define DX    (int)     1
#define DY    (int)     1
#define GX    (int)  2048
```

to
```
#define GX_DIM    512
#define GY_DIM    640     /*  GX_DIM*GY_DIM = # rows of H ; WVDIM*XDIM*YDIM = # columns of H */
#define WVDIM (int)    1
#define XDIM  (int)    132
#define YDIM  (int)    132
#define DX    (int)     1
#define DY    (int)     1
#define GX    (int)  512
```

- The following bug in *circulant.c* was found. I changed 

line 170
```
      while ('\n' != *fgets(b,64,file0)){           // column entries
```

to
```
      while (fgets(b,64,file0)){                    // column entries
        if (b[0] == '\n') break;

```

### Ctis

$ gcc -c complex.c
$ gcc -o ctis ctis.c -D REMOVE_NOISE=0 -D REPORT=1 -D 'FOCAL_PLANE="focal_plane"' -D PICTURE=1 -D 'PICTURE_OUT="reconstructed_image"' -D SEED=37 -D EPSILON=1e-8 -D MU=.01 -D LIM=2 -D RESTART=4 -D COMPUTE_PT=1 -D FTHREADS=1 -D THREADS=4 complex.o  -lfftw3_threads -lfftw3 -lm -lpthread
$ ./ctis

- In the command line flag passed to gcc I changed *-D REMOVE_NOISE=1* to *REMOVE_NOISE=0* and *-D FTHREADS=4* to *-D FTHREADS=1*

- In the *param* file I changed
```
#define GX_DIM        2048      // focal plane width <--- pixels
#define GY_DIM        2048      // focal plane height <--- pixels
#define WVDIM  ((int)   24)     // # of wavelengths system matrix calibrated at
#define XDIM   ((int)   81)     // field stop width <--- pixels
#define YDIM   ((int)   90)     // field stop height <--- pixels
#define DX     ((int)    1)
#define DY     ((int)    1)
#define GX     ((int) 2048)     // system matrix 'gap'
```

to 
```
#define GX_DIM        512      // focal plane width <--- pixels
#define GY_DIM        640      // focal plane height <--- pixels
#define WVDIM  ((int)   1)     // # of wavelengths system matrix calibrated at
#define XDIM   ((int)   132)     // field stop width <--- pixels
#define YDIM   ((int)   132)     // field stop height <--- pixels
#define DX     ((int)    1)
#define DY     ((int)    1)
#define GX     ((int) 512)     // system matrix 'gap'
```


- The folowing bug in the command line flag was found. I changed *-D EPSILON=(1e-8)* to *-D EPSILON=1e-8*

- The following bug in *ctis.c* was found. I changed 

line 179
```
    sprintf(b,"A_%d",i+5);
```

to
```
    sprintf(b,"A_%d",i);
```
