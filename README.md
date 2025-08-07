# codeastro-day1

This is a gravitational wave data analysis toolkit built using [GWPy](https://gwpy.github.io/) and [PyCBC](https://pycbc.org/).  
It enables matched filtering, Q-transform analysis, and waveform generation from open GWOSC data.

##Project Structure

codeastro-day1/
├── main.py # CLI interface for running analysis
├── requirements.txt # Dependencies with versions
├── README.md # This file
├── gwtools/ # Analysis module
│ ├── analysis.py # Signal processing and waveform tools
│ └── data.py # Event metadata and fetching



##Features

- Downloads LIGO data from [GWOSC](https://www.gw-openscience.org/)
- Performs:
  - Q-transform visualization
  - Matched filtering
- CLI support for choosing event, detector, duration

##How to Run

```bash
pip install -r requirements.txt
python main.py --event GW150914 --detector H1 --duration 32

