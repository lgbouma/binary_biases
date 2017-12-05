#!/usr/bin/env bash

pdftk \
  rate_density_vs_radius_model_1_brokenx.pdf \
  rate_density_vs_radius_model_2.pdf \
  rate_density_vs_radius_logs_model_2.pdf \
  rate_density_vs_radius_model_3_Zsub2_0.5.pdf \
  rate_density_vs_radius_logs_model_3_Zsub2_0.5.pdf \
  rate_density_vs_radius_model_3_withtext_Zsub2_0.5.pdf \
  output rates_merged.pdf

# OUTDATED: "error cases" plots.

#pdftk \
#  errcases_rate_density_vs_radius_model_1_brokenx.pdf \
#  errcases_rate_density_vs_radius_model_2.pdf \
#  errcases_rate_density_vs_radius_logs_model_2.pdf \
#  errcases_rate_density_vs_radius_model_3_Zsub2_0.5.pdf \
#  errcases_rate_density_vs_radius_logs_model_3_Zsub2_0.5.pdf \
#  errcases_rate_density_vs_radius_model_3_withtext_Zsub2_0.5.pdf \
#  output errcases_merged.pdf


