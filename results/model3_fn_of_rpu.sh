#!/usr/bin/env bash

pdftk \
  rate_density_vs_radius_logs_model_3_Zsub2_0.50_rpu_15.0.pdf \
  rate_density_vs_radius_logs_model_3_Zsub2_0.50_rpu_17.5.pdf \
  rate_density_vs_radius_logs_model_3_Zsub2_0.50_rpu_20.0.pdf \
  rate_density_vs_radius_logs_model_3_Zsub2_0.50_rpu_22.5.pdf \
  rate_density_vs_radius_logs_model_3_Zsub2_0.50_rpu_25.0.pdf \
  rate_density_vs_radius_model_3_withtext_Zsub2_0.50_rpu_15.0.pdf \
  rate_density_vs_radius_model_3_withtext_Zsub2_0.50_rpu_17.5.pdf \
  rate_density_vs_radius_model_3_withtext_Zsub2_0.50_rpu_20.0.pdf \
  rate_density_vs_radius_model_3_withtext_Zsub2_0.50_rpu_22.5.pdf \
  rate_density_vs_radius_model_3_withtext_Zsub2_0.50_rpu_25.0.pdf \
  output rates_fn_of_rpu.pdf

