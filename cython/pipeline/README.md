Cython wrapper for MOAO {#PY_README}
====================================

This wrapper relies on the "common" code of the MOAO project, using only double precision (for now)
This wrapper intend to provide an easy way to create the input files for the MOAO pipeline


dependencies
------------
* common MOAO sources
* python      2.7
* cython      0.25.2
* astropy     1.3.1
* cfitsio 
* fftw3

note that installing Anaconda (https://conda.io/docs/user-guide/install/download.html) will satisfy the requirement for python, cython and astropy

Installation
------------
In setup.py, fill in the path 'LAPACK_ROOT', 'CFITSIO_LIB', FFTW_ROOT  and the name of the lapack library 'libLapack':
~~~{.py}
LAPACK_ROOT= #/path/to/lapack/root
CFITSIO_LIB= #/path/to/libcfitsio.so
FFTW_LIB=    #/path/to/libcfitsio.so
libLapack=   #lapack library ex "lapack", "openblas"
~~~


then run:
~~~
python setup build_ext -i
~~~

Generating the input files
--------------------------
A set of parameter files is generated with the functioncreate_file available in generate_file.py:
~~~
python -i generate_file.py
>>>createFiles(diamtel,nNgs,nLgs,lgs_cst)
~~~
This set of parameters consists in:
    * sys_params.txt (description of the optical system)
    * Dx.fits (command matrix)
    * abs2fi.fits (the influence function of the actuators)
    * idx.fits (the list of valid subapertures)
    * otf_tel.fits (he optical transfer function of the telescope)
