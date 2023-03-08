#!/bin/bash

INTEL=19.0.3.199
PRG=6.0.5
BUPCV=2019.4.2
UPCXXV=2019.3.2
#
module rm PrgEnv-intel
module rm PrgEnv-gnu
module rm PrgEnv-cray
module load PrgEnv-intel/$PRG
module load git
module load darshan
module load cmake

module swap intel/$INTEL
module load gcc

# use the new GASNet-EX bupc and upc++
module load bupc-narrow/$BUPCV
module load upcxx/$UPCXXV

export HIPMER_BUILD_OPTS="-DHIPMER_USE_UPC_ATOMIC_API=1"

module list

export HIPMER_ENV=cori
export CC=$(which cc)
export CXX=$(which CC)
export MPICC=$(which cc)
export MPICXX=$(which CC)
export HIPMER_BUILD_TYPE="Release"
export CORES_PER_NODE=${CORES_PER_NODE:=32}
export HYPERTHREADS=${HYPERTHREADS:=1}
#export HIPMER_BUILD_OPTS=""
export UPC_SHARED_HEAP_SIZE=${UPC_SHARED_HEAP_SIZE:=2500}
export HIPMER_FULL_BUILD=${HIPMER_FULL_BUILD:=0}
export HIPMER_UPC_ALLOCATOR=1
export USE_SBCAST=1
export HIPMER_NO_AVX512F=1
export HIPMER_BUILD_TEST=1
export CHECK_FREE_HUGEPAGES_MB=${CHECK_FREE_HUGEPAGES_MB:=0}
export MIN_NODES_FOR_HUGEPAGE_CHECK=${MIN_NODES_FOR_HUGEPAGE_CHECK:=2}
export MPICH_GNI_MALLOC_FALLBACK=enabled
export HUGETLB_MORECORE=no
export HIPMER_CRAY_MPI_HACK=1

# files to copy to HIPMER_INSTALL (located in .misc_deploy)
export HIPMER_BIN_SCRIPTS="sbatch_cori.sh sbatch_cori-jgi.sh"
# copy to HIPMER_INSTALL and use as upcrun -conf=
#export HIPMER_UPCRUN_CONF=upcrun.conf
export MPICH_GNI_NDREG_ENTRIES=1024

# These are the machine specific environmental variables for this install

export HIPMER_ENV=cori
unset HIPMER_ENV_SCRIPT

# runtime parameters
export HYPERTHREADS=${HYPERTHREADS:=1}

# machine specific method to execute mpi & upc; default memory per thread for UPC
export MPIRUN=${MPIRUN:=/usr/bin/srun}
export UPCRUN=${UPCRUN:=/global/common/cori/ftg/upc/$BUPCV/hsw/intel/PrgEnv-intel-$PRG-$INTEL-narrow/runtime/inst/bin/upcrun -q}

if [ -n "${HIPMER_PROBE_MEM}" ]
then
  unset PHYS_MEM_MB
  unset GASNET_PHYSMEM_MAX
  unset UPC_SHARED_HEAP_SIZE
  unset GASNET_PHYSMEM_NOPROBE
else
  #export PHYS_MEM_MB=${PHYS_MEM_MB:=489387}
  #export GASNET_PHYSMEM_MAX=${GASNET_PHYSMEM_MAX:=489387MB}
  #export UPC_SHARED_HEAP_SIZE=${UPC_SHARED_HEAP_SIZE:=2500}
  export GASNET_PHYSMEM_NOPROBE=${GASNET_PHYSMEM_NOPROBE:=1}
fi

export UPC_PTHREADS_OPT=${UPC_PTHREADS_OPT:=}

# Compilers
export CC=/opt/cray/pe/craype/2.5.15/bin/cc
export CXX=/opt/cray/pe/craype/2.5.15/bin/CC
export MPICC=/opt/cray/pe/craype/2.5.15/bin/cc
export MPICXX=/opt/cray/pe/craype/2.5.15/bin/CC
export UPCC=/global/common/cori/ftg/upc/$BUPCV/hsw/intel/PrgEnv-intel-$PRG-$INTEL-narrow/runtime/inst/bin/upcc
export CMAKE=/usr/common/software/cmake/3.3.2/bin/cmake

# build parameters
export HIPMER_BUILD_TYPE=Release
export HIPMER_BUILD_OPTS="-DHIPMER_USE_UPC_ATOMIC_API=1"
export HIPMER_POST_INSTALL=
export HIPMER_POSIX_SHM=${HIPMER_POSIX_SHM:=}
export HIPMER_SHARED_MEM_PCT=${HIPMER_SHARED_MEM_PCT:=48}
export HIPMER_UPC_ALLOCATOR=${HIPMER_UPC_ALLOCATOR:=1}

