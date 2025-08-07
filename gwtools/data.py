from gwosc import datasets
from gwpy.timeseries import TimeSeries
from pycbc.waveform import get_fd_waveform

def load_event_data():
    """
    Returns a dictionary of gravitational wave events with parameters:
    (GPS time, m1, m2, SNR, distance, chi_eff, total_mass, final_mass)
    """
    return {
        "GW150914": (
            1126259462.4,  # GPS
            35.6,          # m1
            30.6,          # m2
            26.0,          # SNR
            440.0,         # distance [Mpc]
            -0.01,         # chi_eff
            66.2,          # total mass
            63.1           # final mass
        ),
        "GW231123_135430": (
            1384782888.6,
            137.0,
            103.0,
            22.6,
            2200.0,
            0.31,
            238.0,
            225.0
        ),
        "GW200322_091133": (
            1268903511.3,
            38.0,
            11.3,
            4.5,
            3500.0,
            0.27,
            50.0,
            48.0
        ),
        "GW200316_215756": (
            1268431094.1,
            13.1,
            7.8,
            10.3,
            1120.0,
            0.13,
            21.2,
            20.2
        ),
        "GW200311_115853": (
            1267963151.3,
            34.2,
            27.7,
            17.8,
            1170.0,
            -0.02,
            61.9,
            59.0
        )
    }

def get_event_waveforms(m1, m2, distance, delta_f=1.0 / 4, f_lower=20.0, f_final=2048.0):
    """
    Returns frequency-domain waveform (plus polarization only)
    """
    hp, _ = get_fd_waveform(approximant="IMRPhenomD",
                            mass1=m1,
                            mass2=m2,
                            delta_f=delta_f,
                            f_lower=f_lower,
                            f_final=f_final,
                            distance=distance)
    return hp

def get_qtransform_data(detector, gps_time, duration=0.5):
    """
    Returns the q-transform spectrogram around a GW event
    """
    data = TimeSeries.fetch_open_data(detector, gps_time - 16, gps_time + 16)
    qspecgram = data.q_transform(outseg=(gps_time - duration/2, gps_time + duration/2))
    return qspecgram

def get_snr_data(detector, gps_time, m1, m2, distance):
    """
    Returns the time-domain matched filter SNR around a GW event
    """
    from pycbc.filter import matched_filter

    # Load and preprocess data
    data = TimeSeries.fetch_open_data(detector, gps_time - 16, gps_time + 16)
    data = data.highpass(15.0)
    psd = data.psd(4)
    data = data.crop(gps_time - 2, gps_time + 2)

    # Generate waveform
    hp, _ = get_fd_waveform(approximant="IMRPhenomD",
                            mass1=m1,
                            mass2=m2,
                            delta_f=psd.df.value,
                            f_lower=20,
                            f_final=2048,
                            distance=distance)

    # Matched filter
    snr = matched_filter(hp, data.to_pycbc(), psd=psd.to_pycbc(), low_frequency_cutoff=15.0)
    return snr, data

