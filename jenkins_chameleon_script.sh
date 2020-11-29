#!/bin/bash  -e

blas="mkl" # "openblas" or"mkl"
gpu_available="1" # "1" or "0"
use_gpu="1" #"1" or "0"

if [ "$gpu_available" == "0" ];then
    use_gpu="0"
fi

if [ "$blas" = "openblas" ]; then
    echo "using openblas"

    if [ "$gpu_available" == "1" ]; then
        module load ecrc-extras
        module load gcc/5.5.0
        module load cuda/8.0
        module load fftw/3.3.7-gcc-5.5.0
	module load cfitsio/3.420-gcc-5.5.0
	module load openblas/0.2.20-gcc-5.5.0-singlethread
        module load openmpi/3.0.0-gcc-5.5.0
        module load starpu/1.2.3-gcc-5.5.0-openblas-openmpi-3.0.0-cuda-8.0
        module load cmake
    else
        module load ecrc-extras
        module load gcc/5.5.0
        module load fftw/3.3.7-gcc-5.5.0
	module load cfitsio/3.420-gcc-5.5.0
	module load openblas/0.2.20-gcc-5.5.0-singlethread
	module load starpu/1.2.3-gcc-5.5.0-openblas-openmpi-3.0.0
        module load cmake
    fi

    CHAMELEON_BUILD=build_openblas
    CHAMELEON_BLAS="-DBLA_VENDOR=Open"
    export OMP_NUM_THREADS=1

elif [ "$blas" = "mkl" ];then
    echo "using mkl"

    if [ "$gpu_available" == "1" ]; then
        module load mkl/2018-initial
        module load ecrc-extras
        module load gcc/5.5.0
        module load fftw/3.3.7-gcc-5.5.0
        module load cfitsio/3.420-gcc-5.5.0
        module load cuda/8.0
        module load openmpi/3.0.0-gcc-5.5.0
        module load starpu/1.2.3-gcc-5.5.0-mkl-openmpi-3.0.0-cuda-8.0
        module load cmake
    else
        module load mkl/2018-initial
        module load ecrc-extras
        module load gcc/5.5.0
        module load fftw/3.3.7-gcc-5.5.0
        module load cfitsio/3.420-gcc-5.5.0
        module load starpu/1.2.3-gcc-5.5.0-mkl-openmpi-3.0.0
        module load cmake
    fi

    CHAMELEON_BUILD=build_mkl
    CHAMELEON_BLAS="-DBLA_VENDOR=Intel10_64lp_seq"
else
    echo "first argument must be 'openblas' or 'mkl'"
    return 
fi

# Check if we are already in MOAO repo dir or not.
echo "=================================================="
echo " GET moao_SOURCES"
echo "=================================================="
if git -C $PWD remote -v | grep -q 'https://github.com/ecrc/moao-dev.git'
then
    echo "already in repo"
    # we are, lets go to the top dir (where .git is)
    until test -d $PWD/.git ;
    do
        cd ..
    done;
else
    echo "clone repo"
    #we are not, we need to clone the repo
    git clone https://github.com/ecrc/moao-dev.git MOAO-git
    cd MOAO-git
fi
MOAO_ROOT=$PWD
echo "==================================================\n"

#move to chameleon-repo if exists
echo "=================================================="
echo " GET CHAMELEON SOURCES"
echo "=================================================="
#Use chameleon internal packages
export CHAMELEONDIR=$MOAO_ROOT"/chameleon-repo"
if [ ! -d $CHAMELEONDIR ]; then
    #need to clone the repo
    #clone chameleon
    mkdir $CHAMELEONDIR
    git clone https://gitlab.inria.fr/solverstack/chameleon.git $CHAMELEONDIR
else
    if git -C $CHAMELEONDIR remote -v | grep -q 'https://gitlab.inria.fr/solverstack/chameleon.git'
    then
        echo "already contains sources"
    else
        #folder exists but does not contain sources
        # clean and clone again (should not append)
        clean and clone
        cd $CHAMELEONDIR/..
        rm -rf $CHAMELEONDIR
    	mkdir $CHAMELEONDIR
    	git clone https://gitlab.inria.fr/solverstack/chameleon.git $CHAMELEONDIR
    fi
fi
# Update submodules
git submodule init
git submodule update
cd $CHAMELEONDIR
# Check if we are already in chameleon repo dir or not.
echo "==================================================\n"

#install chameleon
echo "=================================================="
echo " BUILDING CHAMELEON"
echo "=================================================="
CHAMELEON_INSTALL=$CHAMELEONDIR/$CHAMELEON_BUILD/install
mkdir -p $CHAMELEON_INSTALL
cd $CHAMELEONDIR/$CHAMELEON_BUILD

#if [ "$gpu_available" == "0" ];then
#    cmake .. -DCMAKE_INSTALL_PREFIX=$CHAMELEON_INSTALL -DCHAMELEON_USE_CUDA=OFF -DCHAMELEON_USE_MPI=OFF -DBUILD_SHARED_LIBS=ON $CHAMELEON_BLAS
#else
#    cmake .. -DCMAKE_INSTALL_PREFIX=$CHAMELEON_INSTALL -DCHAMELEON_USE_CUDA=ON  -DCHAMELEON_USE_MPI=OFF -DBUILD_SHARED_LIBS=ON $CHAMELEON_BLAS
#fi
#make -j
#make install
export PKG_CONFIG_PATH=$CHAMELEON_INSTALL/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$CHAMELEON_INSTALL/lib/:$LD_LIBRARY_PATH

echo "=================================================="
echo " BUILDING MOAO_CHAMELEON"
echo "=================================================="
#install MOAO_chameleon
cd $MOAO_ROOT/chameleon
if [ "$use_gpu" == "1" ] ; then
echo "USE_GPU"
module li
    CHAMELEON_BUILD=${CHAMELEON_BUILD}_gpu
    mkdir -p ${CHAMELEON_BUILD}
    cd ${CHAMELEON_BUILD}
cmake ../.. -Dproject=chameleon  -DCFITSIO_LIB=$CFITSIO_ROOT/lib -DCFITSIO_INC=$CFITSIO_ROOT/include -DFFTW_LIB=$FFTW_ROOT/lib -DFFTW_INC=$FFTW_ROOT/include -DUSE_INTERSAMPLE=ON -DUSE_GPU=ON -DCUDA_GENCODE="60 35 "
else
echo "USE_CPU"
    mkdir -p ${CHAMELEON_BUILD}
    cd ${CHAMELEON_BUILD}
cmake ../.. -Dproject=chameleon  -DCFITSIO_LIB=$CFITSIO_ROOT/lib -DCFITSIO_INC=$CFITSIO_ROOT/include -DFFTW_LIB=$FFTW_ROOT/lib -DFFTW_INC=$FFTW_ROOT/include -DUSE_INTERSAMPLE=ON -DUSE_GPU=OFF
fi
#compile MOAO_chameleon
sed -i -e '32s/float/double/g' ${MOAO_ROOT}/chameleon/check_moao_chameleon.cpp
make VERBOSE=1
export LD_LIBRARY_PATH=${MOAO_ROOT}/chameleon/${CHAMELEON_BUILD}/common:${MOAO_ROOT}/chameleon/${CHAMELEON_BUILD}/chameleon:$LD_LIBRARY_PATH

#get data for testing
DATAPATH=$MOAO_ROOT"/datafile/"
DATAFILE="jenkins_input"
mkdir -p $DATAPATH/$DATAFILE
cd $DATAPATH
wget --quiet --no-check-certificate  "https://drive.google.com/uc?export=download&id=0Bw6iRA3hQZNCVEtVRjA1Q2xwM00" -O $DATAFILE.tgz
tar -zxf $DATAFILE.tgz -C $DATAFILE
cd $DATAFILE
mkdir -p nx10908_wfs8_Nssp40/2Layers  nx5328_wfs8_Nssp28/2Layers  nx532_wfs6_Nssp10/2Layers
mkdir -p nx10908_wfs8_Nssp40/10Layers  nx5328_wfs8_Nssp28/10Layers  nx532_wfs6_Nssp10/10Layers

cp nx532_wfs6_Nssp10/psf_nx532_10layers.fits nx532_wfs6_Nssp10/10Layers/psf.fits
cp nx532_wfs6_Nssp10/psf_nx532_2layers.fits  nx532_wfs6_Nssp10/2Layers/psf.fits
cp nx5328_wfs8_Nssp28/psf_nx5328_10layers.fits nx5328_wfs8_Nssp28/10Layers/psf.fits
cp nx5328_wfs8_Nssp28/psf_nx5328_2layers.fits  nx5328_wfs8_Nssp28/2Layers/psf.fits
cp nx10908_wfs8_Nssp40/psf_nx10908_10layers.fits nx10908_wfs8_Nssp40/10Layers/psf.fits
cp nx10908_wfs8_Nssp40/psf_nx10908_2layers.fits  nx10908_wfs8_Nssp40/2Layers/psf.fits

#testing
echo "=================================================="
echo " RUNNING MOAO_CHAMELEON"
echo "=================================================="

nthreads=`lscpu | grep "^Thread" | cut -d: -f 2 | tr -d " "`
nsockets=`grep "^physical id" /proc/cpuinfo | sort | uniq | wc -l`
ncorepersocket=`grep "^core id" /proc/cpuinfo | sort  | uniq | wc -l`
totalcpu=`grep "^processor" /proc/cpuinfo | wc -l` # total ONLINE cpus
cputhreads=$(( nthreads * nsockets * ncorepersocket ))
cpucores=$(( nsockets * ncorepersocket ))

ngpu=`nvidia-smi -L | wc -l`

if [ "$use_gpu" == "1" ] ; then
    cd $MOAO_ROOT/chameleon/${CHAMELEON_BUILD}
    ncpu=$((cpucores - 1 - ngpu ))
else
    cd $MOAO_ROOT/chameleon/${CHAMELEON_BUILD}
    ncpu=$((cpucores-1))
fi


sys_path=${DATAPATH}/${DATAFILE}
atm_path=${DATAPATH}/${DATAFILE}/profiles
ts=2000
acc=1.e-3

echo "PWD:    $(pwd)"
echo "LD_LIB: $LD_LIBRARY_PATH"
for fs in $(cd ${sys_path} && ls -rd nx*); do
    sed -i -e '18s/\.//g' ${sys_path}/${fs}/sys-params.txt
    for af in $(cd ${atm_path} && ls -rd *);do
        out=${fs}_${af}
        echo $out
        echo "./check_moao --sys_path=${sys_path}/${fs}/ --atm_path=${atm_path}/${af}/ --sys_check=${sys_path}/${fs}/${af}/  --out_path=./ --ngpus=${ngpu} --ncores=${ncpu} --ts=${ts}" >log_${out} 2>err_${out}
        ./check_moao --sys_path=${sys_path}/${fs}/ --atm_path=${atm_path}/${af}/ --sys_check=${sys_path}/${fs}/${af}/  --out_path=./ --ngpus=${ngpu} --ncores=${ncpu} --ts=${ts} --accuracy=${acc} >>log_${out} 2>err_${out}
        mv psf.fits psf_${out}.fits
    done
done


#clean input data
#rm $DATAPATH/$DATAFILE.tgz
#rm -r $DATAPATH/$DATAFILE
