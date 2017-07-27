#!/bin/bash
# -------------------------------------------------------------------------------
# Programmer(s):  David J. Gardner, Christopher J. Vogl @ LLNL
# -------------------------------------------------------------------------------
# Driver for plot results from KiD
# -------------------------------------------------------------------------------
# Last updated: 19 May 2017
# -------------------------------------------------------------------------------

RESULTSDIR=${S}/KiD
TEST="warm1"
LABEL="warm1_Test"

if [ 1 == 1 ]; then
    python KiD_plot.py \
        ${LABEL} cloud_mass_path \
        ${RESULTSDIR}/${TEST} \
        --LegendOutside \
        --stepsizes 1 5 10 15 30 60 120 300 600 900 1200

    python KiD_plot.py \
        ${LABEL} ice_mass_path \
        ${RESULTSDIR}/${TEST} \
        --LegendOutside \
        --stepsizes 1 5 10 15 30 60 120 300 600 900 1200

    python KiD_plot.py \
        ${LABEL} rain_mass_path \
        ${RESULTSDIR}/${TEST} \
        --LegendOutside \
        --stepsizes 1 5 10 15 30 60 120 300 600 900 1200

    python KiD_plot.py \
        ${LABEL} snow_mass_path \
        ${RESULTSDIR}/${TEST} \
        --LegendOutside \
        --stepsizes 1 5 10 15 30 60 120 300 600 900 1200
fi

if [ 1 == 1 ]; then
    python KiD_convergence.py \
        ${LABEL} \
        ${RESULTSDIR}/${TEST}/${TEST}_mg2_acme_v1beta_dt1.0_mstep1.nc \
        --outfiles \
        "${RESULTSDIR}/${TEST}/${TEST}_mg2_acme_v1beta_dt1.0_mstep5.nc" \
        "${RESULTSDIR}/${TEST}/${TEST}_mg2_acme_v1beta_dt1.0_mstep10.nc" \
        "${RESULTSDIR}/${TEST}/${TEST}_mg2_acme_v1beta_dt1.0_mstep15.nc" \
        "${RESULTSDIR}/${TEST}/${TEST}_mg2_acme_v1beta_dt1.0_mstep30.nc" \
        "${RESULTSDIR}/${TEST}/${TEST}_mg2_acme_v1beta_dt1.0_mstep60.nc" \
        "${RESULTSDIR}/${TEST}/${TEST}_mg2_acme_v1beta_dt1.0_mstep120.nc" \
        "${RESULTSDIR}/${TEST}/${TEST}_mg2_acme_v1beta_dt1.0_mstep300.nc" \
        "${RESULTSDIR}/${TEST}/${TEST}_mg2_acme_v1beta_dt1.0_mstep600.nc" \
        "${RESULTSDIR}/${TEST}/${TEST}_mg2_acme_v1beta_dt1.0_mstep900.nc" \
        "${RESULTSDIR}/${TEST}/${TEST}_mg2_acme_v1beta_dt1.0_mstep1200.nc" \
        --stepsizes 5 10 15 30 60 120 300 600 900 1200 \
        --WriteArrays --Warm
fi
