
# \section*{Aftershock decay rates and power-law fits (50 points)}
# For this assignment, submit a short report with figures documenting your result! You will also have to submit your Python code with comments.

# \section*{Required submission files:}
# \begin{enumerate}
#   \item Report (including the following figures: Prague seismicity map, aftershock rates, map of injection wells and earthquakes in OK, seismicity rate in Oklahoma)
# \end{enumerate}

# \section*{2) Python code}
# The vast majority of earthquakes are aftershocks, i.e. earthquakes that are caused by stress perturbations of preceding events. Thus, the study of aftershocks is of key importance in seismology and provides a critical component for operational earthquake forecasts and shortterm hazard assessments. Aftershocks are special in the sense that the average expected rate is well described by an empirical relationship called Omori's law:

# $$
# \mathrm{d} N / \mathrm{d} t=K /(c+t)^{p}
# $$

# where $\mathrm{d} N / \mathrm{d} t$ is the event rate after the mainshock, $t$ is time after the mainshock, and $K, c$ and $p$ are constants describing the general productivity and decay rate of the aftershock sequence. This empirical relationship provides the opportunity to predict the average expected rate of earthquakes from hours to years after the mainshock.

# Consider the record of seismic events within the central U.S., constrained be the following parameters:

# \section*{Table 1:}
# \begin{center}
# \begin{tabular}{|l|l|l|l|}
# \hline
# Time range & MAG range & Longitude in $^{\circ}$ & Latitude in $^{\circ}$ \\
# \hline
# $2010 / 01 / 01$ to 2012/12/31 & 2.0 to 7 & -100 to -95 & 34.3 to 37 \\
# \hline
# \end{tabular}
# \end{center}

# \begin{center}
# \includegraphics[max width=\textwidth]{2025_10_30_289cd0643bcc3ee09e54g-2}
# \end{center}

# Exemplary figure of the M5.7 Prague earthquake as well as fore and aftershocks in map view. The mainshock is highlighted by the yellow star and focal mechanism labeled by the letter B. Two additional M5 events are labeled A and C .\\
# from: Keranen, K.M., Savage, H.M., Abers, G. a., and Cochran, E.S., 2013, Potentially induced earthquakes in Oklahoma, USA: Links between wastewater injection and the 2011 Mw 5.7 earthquake sequence: Geology, v. 41, no. 6, p. 699-702, doi: 10.1130/G34045.1.

# \section*{Download:}
# \begin{enumerate}
#   \item Download the earthquake catalog from:\\
# \href{https://earthquake.usgs.gov/earthquakes/map/}{https://earthquake.usgs.gov/earthquakes/map/}, using the parameter specified in Table 1 (you can also go to: \href{https://earthquake.usgs.gov/earthquakes/search/}{https://earthquake.usgs.gov/earthquakes/search/}). You can either go to the website and input the parameters manually or you can use the python model os.system and the bash command 'wget'. You should end up with an ASCII data file that contains 364 events.
# \end{enumerate}

# \section*{Data processing:}
# \begin{enumerate}
#   \setcounter{enumi}{1}
#   \item Write a script that imports the earthquake catalog! This is a nice exercise because you are dealing with somewhat complex, multivariate data that includes date, time, magnitude type and regular floating point numbers. You can focus on loading the first 6 columns (Date to Mag). Useful commands for this exercise are numpy.loadtxt, numpy.genfromtxt and standard ASCII file imports with file\_obj = open( [file\_name], 'r)
# \end{enumerate}

# Another way to load ASCII tables is with Pandas, a Python package that is designed to mimic sql-like queries of complex relational databases. Pandas finds much application in machine learning and data science and you should be familiar with the basic syntax, such as:

# \begin{verbatim}
# >>>import pandas as pd
# >>>df = pd.read_csv( [file], use_cols =.., delim_whitespace=True)
# \end{verbatim}

# \begin{enumerate}
#   \setcounter{enumi}{2}
#   \item Select the data of interest:\\
# a. Look at the downloaded datafile and identify the M5.7 Prague main-shock in 2011\\
# b. Spatial Filtering: Select the area of interest using the following empirical scaling relation to identify the area of influence of a mainshock:
# \end{enumerate}

# $$
# R_{\max }=10^{(0.25 M-0.22)}[\mathrm{km}]^{*}
# $$

# There are several ways to filter the data to only include events within $R_{\text {max }}$ such as:\\
# i) create an array of Boolean values (syntax is identical to MATLAB, use numpy.logical\_and and numpy.logical\_or for upper/lower bound conditions)\\
# ii) use the numpy.where() function\\
# *From: Gardner, J.K., and Knopoff, L., 1974, Is the sequence of earthquakes in southern California, with aftershocks removed, Poissonian? Bull. Seismol. Soc. Am, v. 64, no. 5, p. 1363-1367.

# You can project the data from spherical into Cartesian or into an equal distance projection, e.g., using matplotlib basemap. Alternatively, you can directly use the haversine formula to compute distances between your lon/lat vectors from the seis\_utis module (\href{https://en.wikipedia.org/wiki/Haversine}{https://en.wikipedia.org/wiki/Haversine} formula).\\
# c. Time Filtering: Select only events after the mainshock.\\
# d. Convert date-time to days relative to the mainshock occurrence. You can either convert all datetime to decimal years or use the datetime function. The time difference is then:\\
# $t_{\mathrm{as}}=t_{\mathrm{eq}}-t_{\mathrm{MS}}$, where $t_{\mathrm{as}}$ are the aftershock times (in days) relative to the mainshock, $t_{\mathrm{MS}}$\\
# Note that Pandas will automatically recognize DateTime strings and create a TimeDelta object when subtracting the respective origin times of all events from the mainshock origin time.

# Aftershock rates and rate-decay fitting:\\
# 4. Compute the aftershock rate decay using the following equation for linear density estimates:

# Linear Density (forward looking window):

# $$
# \rho_{i}^{e q}=\frac{k}{t_{i+k}-t_{i}} ; \text { with } i=0 \text { to }(n-k)^{* *} ;
# $$

# Centered:

# $$
# \rho_{i+\frac{k}{2}}^{e q}=\frac{k}{t_{i+k}-t_{i}} ; \text { with } i=0 \text { to }(n-k) ;
# $$

# where $k$ is the sample size with higher values of $k$ resulting in smoother density estimates, and $t_{\mathrm{i}+\mathrm{k}}-t_{\mathrm{i}}$ is the distance between the $\mathrm{i}^{\text {th }}$ and the $\mathrm{i}+\mathrm{k}^{\text {th }}$ event. This is essentially a sliding sample window from your first to your $n-k^{\text {th }}$ event with a step size of 1 .\\
# **See: Silverman, B.W., 1986, DENSITY ESTIMATION FOR STATISTICS AND DATA ANALYSIS: Monographs on Statistics and Applied Probability, p. 1-22.\\
# 5. Compare the temporal decay from the above equation to a simple histogram of events after the mainshock! What are some advantages of using statistical density estimates vs. a binned representation of the data?\\
# 6. Fit the power-law portion of the time decay using:

# $$
# \mathrm{d} N / \mathrm{d} t \sim t^{p}
# $$

# Note that there are several steps involved in determining the fit, including: i) logtransformation of the data, ii) selecting the correct time range over which the data shows power-law behavior and iii) finally determining the fit in a least squares sense. You can use:

# \begin{verbatim}
#     >>>import scipt.stats
#     >>>slope, interc, R2, p, err = scipy.stats.linregress(
# np.log10(a_t), np.log10( a_rate) )
# \end{verbatim}

# \begin{enumerate}
#   \setcounter{enumi}{6}
#   \item Plot the fit and the data on double logarithmic scales! Report the aftershock decay exponent, $p$, and the R2-value.
#   \item Discuss your results, observations and implications!
# \end{enumerate}

#%%============================================================================
#                  Import modules and Set defaults parameters
# =============================================================================

import pandas as pd
import requests
import matplotlib.pyplot as plt
import pygmt
from pathlib import Path
from sklearn.metrics.pairwise import haversine_distances  

#%%============================================================================
#                 Input parameters
# =============================================================================

data_file = 'earthquake_catalog.csv'  # Path to the earthquake catalog file
startdate = '2010-01-01'
enddate = '2012-12-31'
min_mag = 2.0
max_mag = 7.0
lon_min = -100.0
lon_max = -95.0
lat_min = 34.3
lat_max = 37.0
format = 'csv'  # File format for download json, geojson, xml, csv

#%%============================================================================
#                 Manage directories
# =============================================================================

work_dir = Path.cwd()
data_dir = work_dir / 'data'
data_dir.mkdir(parents=True, exist_ok=True)

#%%============================================================================
#                 Get earthquake catalog
#==============================================================================

url = (f'https://earthquake.usgs.gov/fdsnws/event/1/query?'
       f'starttime={startdate}&endtime={enddate}&minmagnitude={min_mag}'
       f'&maxmagnitude={max_mag}&minlongitude={lon_min}&maxlongitude={lon_max}'
       f'&minlatitude={lat_min}&maxlatitude={lat_max}&format={format}')

response = requests.get(url)
data_path = data_dir / data_file
with open(data_path, 'wb') as file:
    file.write(response.content)
print(f'Downloaded earthquake catalog to {data_path}')

#%%============================================================================
#                 Load earthquake catalog
# =============================================================================

df = pd.read_csv(data_path, 
                 usecols=['time', 'latitude', 'longitude', 'depth', 'mag'],
                 parse_dates=['time']) # Parse 'time' column as datetime objects
                                        # it makes date-time operations easier
print(f'Loaded {len(df)} events from the catalog.')

print(df.head())

#%%============================================================================
#                 Identify mainshock and filter data
# =============================================================================