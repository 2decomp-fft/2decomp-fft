#!/bin/bash
set -euo pipefail

LOCAL_RANK=${OMPI_COMM_WORLD_LOCAL_RANK:-${MPI_LOCALRANKID:-${SLURM_LOCALID:-0}}}
export LOCAL_RANK

GPU_BIND_BACKEND=${GPU_BIND_BACKEND:-rocm}
SELECTED_GPU=${LOCAL_RANK}

case "${GPU_BIND_BACKEND}" in
  cuda)
    export CUDA_VISIBLE_DEVICES=${SELECTED_GPU}
    unset HIP_VISIBLE_DEVICES
    unset ROCR_VISIBLE_DEVICES
    ;;
  rocm)
    export ROCR_VISIBLE_DEVICES=${SELECTED_GPU}
    unset HIP_VISIBLE_DEVICES
    unset CUDA_VISIBLE_DEVICES
    ;;
  *)
    echo "[ERROR] GPU_BIND_BACKEND must be either 'rocm' or 'cuda' (got '${GPU_BIND_BACKEND}')" >&2
    exit 2
    ;;
esac

# After masking to a single visible GPU, OpenMP should use logical device 0.
export OMP_DEFAULT_DEVICE=0

echo "[LOG] local rank ${LOCAL_RANK}: backend=${GPU_BIND_BACKEND} physical_gpu=${SELECTED_GPU} CUDA=${CUDA_VISIBLE_DEVICES:-unset} HIP=${HIP_VISIBLE_DEVICES:-unset} ROCR=${ROCR_VISIBLE_DEVICES:-unset} OMP=${OMP_DEFAULT_DEVICE}"
echo ""

exec "$@"
