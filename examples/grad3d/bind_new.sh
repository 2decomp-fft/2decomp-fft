#!/bin/bash
set -euo pipefail

GLOBAL_RANK="${OMPI_COMM_WORLD_RANK:-${PMI_RANK:-${SLURM_PROCID:-?}}}"
LOCAL_RANK="${OMPI_COMM_WORLD_LOCAL_RANK:-${MPI_LOCALRANKID:-${SLURM_LOCALID:-0}}}"
LOCAL_SIZE="${OMPI_COMM_WORLD_LOCAL_SIZE:-${MPI_LOCALNRANKS:-${PMI_LOCAL_SIZE:-${SLURM_NTASKS_PER_NODE:-1}}}}"
NPROCS="${OMPI_COMM_WORLD_SIZE:-${PMI_SIZE:-${SLURM_NPROCS:-${SLURM_NTASKS:-1}}}}"

if [[ -z "${OMP_NUM_THREADS:-}" ]]; then
  export OMP_NUM_THREADS=1
fi

if [[ -z "${RANK_STRIDE:-}" ]]; then
  RANK_STRIDE=$((96 / LOCAL_SIZE))
fi

if [[ -z "${OMP_STRIDE:-}" ]]; then
  OMP_STRIDE=1
fi

if [[ -z "${CPU_SHIFT:-}" ]]; then
  CPU_SHIFT=0
fi

if [[ -z "${EXTRA_CPU:-}" ]]; then
  EXTRA_CPU=0
fi

if [[ -z "${GPUS_PER_NUMA:-}" ]]; then
  GPUS_PER_NUMA=6
fi

if [[ -z "${CORES_PER_NUMA:-}" ]]; then
  CORES_PER_NUMA=24
fi

# Determine the currently visible GPU list from existing env or Slurm metadata.
VISIBLE_GPUS="${ROCR_VISIBLE_DEVICES:-${HIP_VISIBLE_DEVICES:-${CUDA_VISIBLE_DEVICES:-${SLURM_STEP_GPUS:-${SLURM_JOB_GPUS:-}}}}}"
CPU_BASES=(0 24 48 72)

if [[ -n "${VISIBLE_GPUS}" ]]; then
  IFS=',' read -r -a GPU_LIST <<< "${VISIBLE_GPUS}"
  if (( ${#GPU_LIST[@]} > 1 )); then
    GPU_INDEX=$(( LOCAL_RANK % ${#GPU_LIST[@]} ))
    SELECTED_GPU="${GPU_LIST[GPU_INDEX]}"
  else
    SELECTED_GPU="${GPU_LIST[0]}"
  fi
  echo "[LOG] rank=${GLOBAL_RANK} local=${LOCAL_RANK}: bind from visible list -> GPU ${SELECTED_GPU}"
else
  SELECTED_GPU="${LOCAL_RANK}"
  echo "[LOG] rank=${GLOBAL_RANK} local=${LOCAL_RANK}: bind by local rank -> GPU ${SELECTED_GPU}"
fi

GPU_BIND_BACKEND="${GPU_BIND_BACKEND:-rocm}"
case "${GPU_BIND_BACKEND}" in
  cuda)
    export CUDA_VISIBLE_DEVICES="${SELECTED_GPU}"
    unset HIP_VISIBLE_DEVICES
    unset ROCR_VISIBLE_DEVICES
    ;;
  rocm)
    export ROCR_VISIBLE_DEVICES="${SELECTED_GPU}"
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

if [[ "${SELECTED_GPU}" =~ ^[0-9]+$ ]]; then
  GPU_NUMA=$((SELECTED_GPU / GPUS_PER_NUMA))
  GPU_LOCAL_INDEX=$((SELECTED_GPU % GPUS_PER_NUMA))
  CPU_SLOT_STRIDE=$((CORES_PER_NUMA / GPUS_PER_NUMA))
  if (( GPU_NUMA < ${#CPU_BASES[@]} )); then
    CPU_START=$((CPU_BASES[GPU_NUMA] + GPU_LOCAL_INDEX * CPU_SLOT_STRIDE + CPU_SHIFT))
  else
    CPU_START=$((RANK_STRIDE * LOCAL_RANK + CPU_SHIFT))
  fi
else
  CPU_START=$((RANK_STRIDE * LOCAL_RANK + CPU_SHIFT))
fi
CPU_STOP=$((CPU_START + OMP_NUM_THREADS * OMP_STRIDE - 1 + EXTRA_CPU))
export GOMP_CPU_AFFINITY="${CPU_START}-${CPU_STOP}:${OMP_STRIDE}"

echo "[LOG] rank=${GLOBAL_RANK} local=${LOCAL_RANK}: backend=${GPU_BIND_BACKEND} physical_gpu=${SELECTED_GPU} gpu_numa=${GPU_NUMA:-unknown} CUDA=${CUDA_VISIBLE_DEVICES:-unset} HIP=${HIP_VISIBLE_DEVICES:-unset} ROCR=${ROCR_VISIBLE_DEVICES:-unset} OMP=${OMP_DEFAULT_DEVICE}"
echo "[LOG] rank=${GLOBAL_RANK} local=${LOCAL_RANK}: cpu=${GOMP_CPU_AFFINITY} local_size=${LOCAL_SIZE}"

echo
exec "$@"
