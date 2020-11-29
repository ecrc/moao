#!/bin/bash
#SBATCH --account=k1217
#SBATCH --job-name=moao
#SBATCH --output=moao-%j.out
#SBATCH --error=moao-%j.err
#SBATCH --nodes=1
#SBATCH --time=00:30:00

#module swap PrgEnv-cray PrgEnv-intel
export CRAYPE_LINK_TYPE=dynamic
module load perftools-base
export MKL_NUM_THREADS=1 OMP_NUM_THREADS=1;

#grid_order -C -n 32 -c 4,8 -g $nprow,$npcol,$nbpbs > MPICH_RANK_ORDER
#export MPICH_RANK_REORDER_METHOD=3
#srun --ntasks=$ntasks --hint=nomultithread --ntasks-per-node=32 --ntasks-per-socket=16 ./moao_scalapack --nprow=4 --npcol=8 --nb=4 --filePath="/project/k1217/ltaiefh/codes/moao-dev/datafile/check/nx456_nLayers10_wfs6_Nssp10/"  --suffix="_nx456_nLayers10_wfs6_Nssp10" --maxrefine=1 --maxobs=1 --v
srun --ntasks=32 --hint=nomultithread --ntasks-per-node=32 --ntasks-per-socket=16 ./test_moao_scalapack --nprow=4 --npcol=8 --nb=4 --filePath="/project/k1217/ltaiefh/codes/moao-dev/datafile/check/nx456_nLayers10_wfs6_Nssp10/"  --suffix="_nx456_nLayers10_wfs6_Nssp10" --maxrefine=1 --maxobs=1 --output --v
echo "== Node lists:", $SLURM_JOB_NODELIST
