#!/bin/sh

DATA=${1}
FILE=${2}

onmt_preprocess \
-train_src ${DATA}/${FILE}.tagged_filtered.train.products_tokens \
-train_tgt ${DATA}/${FILE}.tagged_filtered.train.precursors_tokens \
-valid_src ${DATA}/${FILE}.tagged_filtered.validation.products_tokens \
-valid_tgt ${DATA}/${FILE}.tagged_filtered.validation.precursors_tokens \
-save_data ${DATA}/${FILE} \
-src_seq_length 1000 -tgt_seq_length 1000 \
-src_vocab_size 1000 -tgt_vocab_size 1000 -share_vocab