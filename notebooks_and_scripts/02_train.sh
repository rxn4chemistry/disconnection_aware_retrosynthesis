#!/bin/sh

DATA=${1}
SAVE_MODEL=${2}

mkdir ${DATA}/logs
onmt_train \
-data ${DATA}/${SAVE_MODEL} \
-save_model ${DATA}/${SAVE_MODEL} \
-seed 42 \
-gpu_ranks 0 \
-save_checkpoint_steps 5000 \
-keep_checkpoint 20 \
-train_steps 260000 \
-param_init 0 \
-param_init_glorot \
-max_generator_batches 32 \
-batch_size 6144 \
-batch_type tokens \
-normalization tokens \
-max_grad_norm 0 \
-accum_count 4 \
-optim adam \
-adam_beta1 0.9 \
-adam_beta2 0.998 \
-decay_method noam \
-warmup_steps 8000  \
-learning_rate 2 \
-label_smoothing 0.0 \
-report_every 1000  \
-valid_batch_size 8 \
-layers 4 \
-rnn_size 384 \
-word_vec_size 384 \
-encoder_type transformer \
-decoder_type transformer \
-dropout 0.1 \
-position_encoding -share_embeddings -global_attention general \
-global_attention_function softmax -self_attn_type scaled-dot \
-heads 8 -transformer_ff 2048 \
--tensorboard --tensorboard_log_dir ${DATA}/logs
