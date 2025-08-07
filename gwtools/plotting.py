import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from gwpy.timeseries import TimeSeries


def plot_q_transform(qspecgram, t0, freq_ylim=(20, 500), figsize=(8, 4)):
    """
    Plot a Q-transform spectrogram around a given time.
    """
    plot = qspecgram.plot(figsize=figsize)
    ax = plot.gca()
    ax.set_xscale('seconds')
    ax.set_yscale('log')
    ax.set_ylim(freq_ylim)
    ax.set_ylabel('Frequency [Hz]')
    ax.grid(True, axis='y', which='both')
    ax.colorbar(cmap='viridis', label='Normalized energy')
    ax.set_title(f"Q-transform around t0 = {t0}")
    plot.show()


def plot_mass_ratio_distribution(event_dict):
    """
    Plot histogram of mass ratios (q = m2 / m1) from event_dict.
    """
    mass_ratios = []
    for _, (_, m1, m2, *_rest) in event_dict.items():
        if m1 and m2:
            if m1 < m2:
                m1, m2 = m2, m1
            q = m2 / m1
            mass_ratios.append(q)

    plt.figure(figsize=(8, 5))
    plt.hist(mass_ratios, bins=20, color='skyblue', edgecolor='black')
    plt.xlabel("Mass ratio (q = m2 / m1)")
    plt.ylabel("Number of events")
    plt.title("Mass Ratio Distribution of GW Events")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.show()


def plot_m1_over_m2_distribution(event_dict):
    """
    Plot histogram of m1 / m2.
    """
    ratios = []
    for _, (_, m1, m2, *_rest) in event_dict.items():
        if m1 and m2 and not np.isnan(m1) and not np.isnan(m2):
            ratios.append(m1 / m2)

    plt.figure(figsize=(8, 5))
    plt.hist(ratios, bins=20, color='lightcoral', edgecolor='black')
    plt.xlabel("Mass ratio (m1 / m2)")
    plt.ylabel("Number of events")
    plt.title("Distribution of m1 / m2 for GW Events")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.show()


def plot_chi_eff_distribution(event_dict):
    """
    Plot distribution of effective spin χ_eff.
    """
    chi_values = [
        chi_eff for _, (_, _, _, _, _, chi_eff, _, _) in event_dict.items()
        if chi_eff is not None and not np.isnan(chi_eff)
    ]

    mean_val = np.mean(chi_values)
    median_val = np.median(chi_values)

    plt.figure(figsize=(8, 5))
    sns.histplot(chi_values, bins=20, color='lightgreen', edgecolor='black', stat='count', kde=False)
    sns.kdeplot(chi_values, color='darkgreen', linewidth=2, label="KDE")

    plt.axvline(0, color='red', linestyle='--', label="χ_eff = 0")
    plt.axvline(mean_val, color='blue', linestyle='-', label=f"Mean = {mean_val:.2f}")
    plt.axvline(median_val, color='purple', linestyle=':', label=f"Median = {median_val:.2f}")

    plt.xlabel("χ_eff")
    plt.ylabel("Number of events")
    plt.title("χ_eff Distribution of GW Events")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.legend()
    plt.show()


def plot_waveform_comparison(event_name, m1, m2, distance, approximants=None):
    """
    Plot waveform polarizations for different approximants for a given event.
    """
    from pycbc.waveform import get_td_waveform

    if approximants is None:
        approximants = ["IMRPhenomD", "IMRPhenomXAS", "IMRPhenomXPHM"]

    for approx in approximants:
        hp, hc = get_td_waveform(
            approximant=approx,
            mass1=m1,
            mass2=m2,
            distance=distance,
            inclination=0,
            delta_t=1.0 / 4096,
            f_lower=20
        )

        peak_idx = np.argmax(abs(hp))
        hp.start_time -= hp.sample_times[peak_idx]
        hc.start_time -= hc.sample_times[peak_idx]

        mask = (hp.sample_times >= -0.2) & (hp.sample_times <= 0.05)

        plt.figure(figsize=(10, 4))
        plt.plot(hp.sample_times[mask], hp[mask], 'k', label='Plus Polarization (h+)', linewidth=1.5)
        plt.plot(hc.sample_times[mask], hc[mask], 'gray', linestyle='--', label='Cross Polarization (h×)', linewidth=1.5)
        plt.xlabel("Time (s) relative to merger")
        plt.ylabel("Strain")
        plt.title(f"Theoretical Waveform ({approx}) for {event_name}")
        plt.legend(frameon=True)
        plt.grid(True, linestyle="--", alpha=0.6)
        plt.show()


def plot_matched_filter_snr(snr_ts, t0=1126259462.427, xlim=(1126259461, 1126259463)):
    """
    Plot SNR from matched filter.
    """
    plot = snr_ts.plot()
    ax = plot.gca()
    ax.set_xlim(*xlim)
    ax.set_epoch(t0)
    ax.set_ylabel('Signal-to-noise ratio (SNR)')
    ax.set_title('LIGO-Hanford signal-correlation for GW150914')
    plot.show()

