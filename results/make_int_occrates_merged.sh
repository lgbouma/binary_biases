#!/usr/bin/env bash

# general models
pdftk \
  int_occ_rate_vs_radius_model_3_Zsub2_0.50_rpu_22.5.pdf \
  int_occ_rate_vs_radius_logy_model_3_Zsub2_0.50_rpu_22.5.pdf \
  int_occ_rate_vs_radius_model_3_withtext_Zsub2_0.50_rpu_22.5.pdf \
  int_occ_rate_vs_radius_model_4_xcut.pdf \
  output int_occ_rate_merged.pdf

# as function of Z_2
pdftk \
  int_occ_rate_vs_radius_logs_model_3_Zsub2_0.00_rpu_22.5.pdf \
  int_occ_rate_vs_radius_logs_model_3_Zsub2_0.25_rpu_22.5.pdf \
  int_occ_rate_vs_radius_logs_model_3_Zsub2_0.50_rpu_22.5.pdf \
  int_occ_rate_vs_radius_model_3_withtext_Zsub2_0.00_rpu_22.5.pdf \
  int_occ_rate_vs_radius_model_3_withtext_Zsub2_0.25_rpu_22.5.pdf \
  int_occ_rate_vs_radius_model_3_withtext_Zsub2_0.50_rpu_22.5.pdf \
  int_occ_rate_vs_radius_model_3_rpu_22.5_manyZs.pdf \
  output int_occ_rate_merged_fn_of_Z2.pdf

#as function of rpu
pdftk \
  int_occ_rate_vs_radius_logy_model_3_Zsub2_0.50_rpu_15.0.pdf \
  int_occ_rate_vs_radius_logy_model_3_Zsub2_0.50_rpu_17.5.pdf \
  int_occ_rate_vs_radius_logy_model_3_Zsub2_0.50_rpu_20.0.pdf \
  int_occ_rate_vs_radius_logy_model_3_Zsub2_0.50_rpu_22.5.pdf \
  int_occ_rate_vs_radius_logy_model_3_Zsub2_0.50_rpu_25.0.pdf \
  int_occ_rate_vs_radius_model_3_withtext_Zsub2_0.50_rpu_15.0.pdf \
  int_occ_rate_vs_radius_model_3_withtext_Zsub2_0.50_rpu_17.5.pdf \
  int_occ_rate_vs_radius_model_3_withtext_Zsub2_0.50_rpu_20.0.pdf \
  int_occ_rate_vs_radius_model_3_withtext_Zsub2_0.50_rpu_22.5.pdf \
  int_occ_rate_vs_radius_model_3_withtext_Zsub2_0.50_rpu_25.0.pdf \
  output int_occ_rate_merged_fn_of_rpu.pdf

