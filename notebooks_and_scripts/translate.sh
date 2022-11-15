#!/bin/sh

DATA=${1}
MODEL=${2}
FILE=${3}
N_BEST=${4}

onmt_translate \
-model ${DATA}/${MODEL}.pt \
-src ${DATA}/${FILE} \
-output ${DATA}/retro_predictions_${MODEL}_top_${N_BEST}.txt \
-batch_size 64 -replace_unk -max_length 200 \
-gpu 0 -n_best ${N_BEST} -beam_size 10
