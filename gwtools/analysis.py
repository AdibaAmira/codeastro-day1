from gwosc import datasets
from gwpy.timeseries import TimeSeries
from gwpy.segments import Segment
from pycbc.waveform import get_fd_waveform
from pycbc.filter import matched_filter
import numpy as np


def fetch_data(detector="H1", gps_start=1126259446, gps_end=1126259478):
    """
    Fetches LIGO time series data from GWOSC.
    """
    data = TimeSeries.fetch_open_data(detector, gps_start, gps_end)
    return data


def fetch_event_data(event_name="GW150914", detector="H1", pre_buffer=30, post_buffer=30):
    """
    Fetch data around a GWOSC event.
    """
    gps = datasets.event_gps(event_name)
    data = TimeSeries.fetch_open_data(detector, gps - pre_buffer, gps + post_buffer)
    return data, gps


def compute_q_transform(data, t0, window=0.3):
    """
    Compute Q-transform zoomed in around a central time t0.
    """
    qspecgram = data.q_transform(outseg=(t0 - window / 2, t0 + window / 2))
    return qspecgram


def compute_qgram(data, gps, qrange=(4, 150), mismatch=0.35):
    """
    Computes Q-gram search around the GPS time.
    """
    search = Segment(gps - 0.25, gps + 0.25)
    qgram = data.q_gram(qrange=qrange, search=search, mismatch=mismatch)
    return qgram


def compute_matched_filter_snr(detector="H1", gps_start=1126259446, gps_end=1126259478,
                                m1=40, m2=32, f_lower=20):
    """
    Computes matched-filter SNR using PyCBC for a given time range and waveform.
    """
    data = TimeSeries.fetch_open_data(detector, gps_start, gps_end)
    high = data.highpass(15)

    psd = high.psd(4, 2)
    zoom = high.crop(1126259460, 1126259464)

    hp, _ = get_fd_waveform(
        approximant="IMRPhenomD",
        mass1=m1,
        mass2=m2,
        f_lower=f_lower,
        f_final=2048,
        delta_f=psd.df.value,
    )

    snr = matched_filter(hp, zoom.to_pycbc(), psd=psd.to_pycbc(),
                         low_frequency_cutoff=f_lower).numpy()

    # Convert to GWPy TimeSeries
    dt = zoom.dt.value
    epoch = zoom.t0.value
    snr_ts = TimeSeries(snr, dt=dt, epoch=epoch).abs()

    return snr_ts

