#!/usr/bin/env bash

pdftk \
  rate_density_vs_radius_logs_model_3_Zsub2_0.00_rpu_22.5.pdf \
  rate_density_vs_radius_logs_model_3_Zsub2_0.25_rpu_22.5.pdf \
  rate_density_vs_radius_logs_model_3_Zsub2_0.50_rpu_22.5.pdf \
  rate_density_vs_radius_model_3_withtext_Zsub2_0.00_rpu_22.5.pdf \
  rate_density_vs_radius_model_3_withtext_Zsub2_0.25_rpu_22.5.pdf \
  rate_density_vs_radius_model_3_withtext_Zsub2_0.50_rpu_22.5.pdf \
  output rates_fn_of_Z2.pdf

