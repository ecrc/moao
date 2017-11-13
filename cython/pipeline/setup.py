#compile with: python setup.py build_ext -i
import  os
from os.path import join as pjoin
import codecs
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy


USE_SINGLE=0

LAPACK_ROOT= #/path/to/lapack/root
CFITSIO_LIB= #/path/to/libcfitsio.so
FFTW_LIB=    #/path/to/libcfitsio.so
libLapack=   #lapack library ex "lapack", "openblas"

LAPACK_INC=LAPACK_ROOT+"include/"
LAPACK_LIB=LAPACK_ROOT+"lib/"


#sources path
CPATH=os.path.dirname(os.path.dirname(os.getcwd()))+"/common"

#get numpy include path
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

os.environ["CC"] = "g++"



parFile = None
if not os.path.isfile("par.pxi"):
    parFile = codecs.open("par.pxi", mode='w', encoding='utf-8')
else:
    import warnings
    warnings.warn("par.pxi found, it will not be updated", Warning)


if(USE_SINGLE==1):
    REAL_T="-DUSE_SINGLE"
else:
    REAL_T="-DUSE_DOUBLE"

if parFile:
    parFile.write("DEF USE_SINGLE=%d # 0/1 \n" % USE_SINGLE)
    parFile.close()


ext_modules=[
    Extension("pipe",
              sources=["pipe.pyx",
                       "../../common/tomo_struct.c",
                       "../../common/noise.c",
                       "../../common/matcov.c",
                       "../../common/matcov_kernels.c",
                       "../../common/utils.c",
                       "../../common/intersample.c"
              ],
              extra_compile_args=[REAL_T,"-DUSE_LAPACK"],
              library_dirs=[CPATH,
                            LAPACK_LIB,
			    CFITSIO_LIB,
			    FFTW_LIB,
                            ],
              libraries=[libLapack,
                        "cfitsio",
                        "fftw3"
                        ],
              language='c++',
              runtime_library_dirs=[CPATH,
                                    LAPACK_LIB,
			            CFITSIO_LIB,
			            FFTW_LIB,
                                   ],
              include_dirs = [numpy_include, 
			      "../../common/",
                              CPATH,
                              LAPACK_INC
                              ]
              #extra_object=[]
    )
]

setup(
  name = "pipe",
  ext_modules = cythonize(ext_modules)
)
