#!/usr/bin/env python3

from gwtools.data import load_event_data, get_event_waveforms, get_qtransform_data, get_snr_data
from gwtools.plotting import (
    plot_q_transform,
    plot_mass_ratio_distribution,
    plot_m1_over_m2_distribution,
    plot_chi_eff_distribution,
    plot_waveform_comparison,
    plot_matched_filter_snr
)

def main():
    # Load event data
    event_dict = load_event_data()

    # --- Part A: Plot histogram-style analyses ---
    plot_mass_ratio_distribution(event_dict)
    plot_m1_over_m2_distribution(event_dict)
    plot_chi_eff_distribution(event_dict)

    # --- Part B: Plot waveform comparisons for GW150914 ---
    event_name = "GW150914"
    gps, m1, m2, snr, distance, chi_eff, total_mass, final_mass = event_dict[event_name]
    print(f"{event_name}: m1 = {m1}, m2 = {m2}, distance = {distance}")
    plot_waveform_comparison(event_name, m1, m2, distance)

    # --- Part C: Plot Q-transform ---
    qdata, qgram, t0 = get_qtransform_data(event_name)
    plot_q_transform(qgram, t0)

    # --- Part D: Matched filter SNR time series ---
    snr_ts = get_snr_data()
    plot_matched_filter_snr(snr_ts)

if __name__ == "__main__":
    main()

