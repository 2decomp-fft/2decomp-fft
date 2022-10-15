#!/bin/bash

export LOCAL_RANK=${OMPI_COMM_WORLD_LOCAL_RANK}
export CUDA_VISIBLE_DEVICES=${LOCAL_RANK}

echo "[LOG] local rank $LOCAL_RANK: bind to $CUDA_VISIBLE_DEVICES"
echo ""

$* 

#if [[ ${LOCAL_RANK} == 0 ]]; then
#        #nsys profile --trace=cuda,nvtx,mpi,openacc --cuda-memory-usage=true --gpu-metrics-device=0 $*
#        nsys profile --trace=cuda,nvtx,mpi,openacc --cuda-memory-usage=true $*
#else
#        $*
#fi
