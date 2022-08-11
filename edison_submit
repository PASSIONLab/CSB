#!/bin/bash -l
#PBS -N test_ccm
#PBS -q ccm_queue
#PBS -l mppwidth=24,walltime=16:00
#PBS -j oe

cd $PBS_O_WORKDIR
module load ccm
module load intel
export CRAY_ROOTFS=DSL
CILK_NWORKERS=12 numactl --cpunodebind=0 ./both $SCRATCH/Z5.N5.Nmax=6.p15.t6.fnond.cus nosym binary csc
CILK_NWORKERS=24 numactl --cpunodebind=0 ./both $SCRATCH/Z5.N5.Nmax=6.p15.t6.fnond.cus nosym binary csc
CILK_NWORKERS=12 numactl --cpunodebind=1 ./both $SCRATCH/Z5.N5.Nmax=6.p15.t6.fnond.cus nosym binary csc
CILK_NWORKERS=24 numactl --cpunodebind=1 ./both $SCRATCH/Z5.N5.Nmax=6.p15.t6.fnond.cus nosym binary csc
