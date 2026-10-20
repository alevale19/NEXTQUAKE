# NEXTQUAKE

# Overview
A MATLAB framework for building an erthquake rupture forecast model taking into account time-dependance and multi-fault rupture inversions via Simulated Annealing.

# Get started – User Guide

## 1. Installation and overview

Install **MATLAB R2025b** and run `NEXTQUAKE.m`. A graphical user interface (GUI) will open, allowing you to specify the parameters required to construct an earthquake rupture forecast (ERF) model.
NEXTQUAKE first generates a catalogue of plausible earthquake ruptures for the study area. It then estimates time-independent (Poissonian) and time-dependent (Brownian Passage Time, BPT) rupture rates by inverting the specified constraints using a simulated annealing (SA) algorithm. Finally, it automatically saves figures and a log file containing information about the run. Optionally, it can also generate input files for subsequent probabilistic seismic hazard analysis (PSHA) using OpenQuake (**version TBD**).

The workflow is organized into three main stages:
1. **Rupture catalogue generation:** construct a catalogue of plausible single-fault and multi-fault ruptures.
2. **Inversion:** estimate rupture occurrence rates by fitting the specified constraints.
3. **Output generation:** save figures, log files, and, optionally, OpenQuake input files and additional rupture-catalogue products.

## 2. Rupture catalogue setup

### 2.1. Input fault database

Place a `.txt` file containing the geometric and kinematic parameters of the fault sections in the `input` folder.

For each fault section, specify the following fields:

| Field | Description | Unit |
|---|---|---|
| `FaultID` | Unique fault-section identifier | — |
| `FaultName` | Fault-section name | — |
| `Dip` | Fault dip angle | ° |
| `Rake` | Fault rake angle | ° |
| `Strike` | Fault strike angle | ° |
| `upperST` | Upper seismogenic depth | km |
| `lowerST` | Lower seismogenic depth | km |
| `CreepFraction` | Fraction of the fault slip accommodated by creep, between 0 and 1 | — |

Ensure that the `MainFaults_coords` folder contains one `.txt` file for each `FaultName`. Each file must contain the longitude and latitude coordinates, in decimal degrees, defining the fault's surface trace.
In the GUI, enter the input filename in the **Input Faults File** field, without the `.txt` extension.

### 2.2. Parameters controlling rupture catalogue generation

After loading the fault database, specify the parameters that control the construction of the rupture catalogue.

**Forecast Mw**

Defines the target magnitude used to determine the spatial discretization of fault sections into subsections. Select this parameter according to the objectives of the forecast and the characteristics of the study region.

**Fault mechanism**

Select the kinematic style used to determine the scaling relationships applied to rupture areas and magnitudes:

- `N`: Normal faulting.
- `R`: Reverse faulting.
- `SS`: Strike-slip faulting.
- `All`: All crustal kinematic styles.
- `SUB`: Subduction-zone faulting.

**Regional b-value**

Specifies the Gutenberg–Richter (GR) b-value used to construct a truncated regional magnitude–frequency distribution (MFD) based on the seismic moment rates of the faults in the database.
The resulting MFD can be used as the regional MFD constraint or as a reference distribution for comparison with the inversion results.
This parameter is optional. If no value is specified, the code uses a default value of 1. However, if a regional MFD input file is provided as a constraint, that file takes precedence over the MFD calculated using the b-value and fault moment rates.

**Maximum shallow gap and maximum deep gap**

Define the maximum allowable distance between subsections belonging to different fault sections for a connection to be considered admissible. The applicable threshold depends on the depth of the subsections.
The transition between the shallow- and deep-gap thresholds is controlled by `z lim` (km).
For example, if the maximum shallow gap is 5 km, the maximum deep gap is 10 km, and `z lim` is 7 km, the shallow-gap threshold applies to subsections shallower than 7 km, whereas the deep-gap threshold applies to subsections at greater depths.

**Maximum azimuth, rake, and dip differences**

Define the maximum allowable differences in azimuth, rake, and dip between fault sections when determining fault clusters.
These criteria control which fault sections are considered geometrically and kinematically compatible within the same fault system and can potentially participate in multi-fault ruptures.
For example, if the maximum azimuth difference is set to 40°, two fault sections whose azimuth difference exceeds 40° cannot be directly connected under this clustering criterion. Consequently, a multi-fault rupture cannot connect those sections directly through that link.

**Maximum aspect ratio**

Defines the maximum allowable ratio between rupture length along strike (\(L\)) and down-dip width (\(W\)):

\[
L/W
\]

This parameter is used to filter out candidate ruptures whose aspect ratios exceed the specified threshold.

**Maximum Mw**

Defines the maximum magnitude allowed in the rupture catalogue. For example, setting the maximum magnitude to 7.5 prevents the code from generating candidate ruptures with magnitudes greater than Mw 7.5.

**Output Catalogue Name**

Specifies the name of the generated rupture catalogue.
Generating the rupture catalogue can be computationally demanding. Once the catalogue has been generated and saved, it can be reused in subsequent runs without repeating the catalogue-generation stage. To do this, enter the existing catalogue filename in the **Pre-processed Rupture Catalog** field, allowing the inversion to be rerun using the previously generated catalogue.

## 3. Inversion data and constraints

### 3.1. Geological slip-rate constraint

The geological slip-rate constraint is mandatory.

Load a `.txt` file containing the following fields:

| Field | Description | Unit |
|---|---|---|
| `FaultName` | Fault-section name | — |
| `SRmin` | Minimum prescribed slip rate | mm/yr |
| `SRmax` | Maximum prescribed slip rate | mm/yr |

Ensure that `FaultName` is consistent across all input files used throughout the workflow.

### 3.2. Optional constraints

The following constraints can be included in the inversion when the corresponding information is available.

**Paleoseismic event rates**

Load a `.txt` file containing:

| Field | Description |
|---|---|
| `FaultName` | Fault-section name |
| `Lon` | Longitude of the paleoseismic site, in decimal degrees |
| `Lat` | Latitude of the paleoseismic site, in decimal degrees |
| `MeanRate` | Mean paleoseismic event rate |
| `MaxRate` | Maximum paleoseismic event rate |
| `MinRate` | Minimum paleoseismic event rate |

The code uses the site coordinates to assign the paleoseismic event-rate constraint to the corresponding fault subsection.

**Regional MFD**

Load a `.txt` file containing the following fields:

| Field | Description |
|---|---|
| `magBin` | Magnitude-bin value |
| `rate` | Incremental occurrence rate for that magnitude bin |

If this file is provided, it takes precedence over the regional MFD calculated internally from the regional b-value and the seismic moment rates of the faults.

**Fault MFD**

Load a `.txt` file containing:

| Field | Description |
|---|---|
| `FaultName` | Fault-section name |
| `bValue` | Target Gutenberg–Richter b-value |

The inversion penalizes departures of the resulting fault-specific magnitude distribution from the specified b-value.

**Segmentation constraint**

Use this constraint to penalize rupture scenarios involving selected pairs of fault sections.

Load a `.txt` file containing:

| Field | Description |
|---|---|
| `FaultA` | Name of the first fault section |
| `FaultB` | Name of the second fault section |

The code penalizes candidate ruptures that involve both specified fault sections. The effect depends on the segmentation-constraint weight; a finite penalty does not necessarily exclude these ruptures entirely.

**Point slip-rate constraint**

Load a `.txt` file containing:

| Field | Description |
|---|---|
| `FaultName` | Fault-section name |
| `Lon` | Longitude, in decimal degrees |
| `Lat` | Latitude, in decimal degrees |
| `SRmin` | Minimum prescribed slip rate, in mm/yr |
| `SRmax` | Maximum prescribed slip rate, in mm/yr |

The code assigns the specified slip-rate bounds to the subsection corresponding to the input coordinates. These bounds override the values assigned to that subsection through the mandatory fault-level slip-rate constraint.

### 3.3. Fault time settings

The fault time settings define the parameters required for the time-dependent BPT calculations.

Load a `.txt` file containing:

| Field | Description |
|---|---|
| `FaultName` | Fault-section name |
| `StartFrac` | Starting fraction of the fault section affected by the historical rupture |
| `EndFrac` | Ending fraction of the fault section affected by the historical rupture |
| `LastEventYr` | Year of the most recent event |
| `Thist` | Historical open interval |
| `CV` | Coefficient of variation of the recurrence-time distribution |

For example, a fault section with `StartFrac = 0` and `EndFrac = 1` is assumed to have ruptured over its entire length during the last event. If `LastEventYr = 1940`, the most recent event occurred in 1940. `Thist` can be set to `NaN` because the historical open interval is not required when the last-event year is known. A value of `CV = 0.5` specifies the coefficient of variation of the recurrence-time distribution.
If both LastEventYr and Thist are unknown and set to NaN, the code assumes a stationary renewal state.

**Forecast Time Window**

Specifies the duration of the forecast window, in years, over which the Poissonian and BPT probabilities are calculated.

**Year for Computing Elapsed Time**

Specifies the reference year used to calculate the elapsed time since the last earthquake.
For example, if LastEventYr = 1940 and the reference year is 2026, the elapsed time is:[2026-1940=86\text{ years}]

## 4. Simulated annealing settings

### 4.1. Constraint weights

The first four fields define the weights assigned to the corresponding constraints relative to the geological slip-rate constraint.
Weights can range from 0 to a sufficiently large value (e.g., (10^6)). A weight of 1 assigns the corresponding constraint a nominal weight equal to that of the slip-rate constraint. A weight of 0 effectively disables the corresponding penalty term.
The appropriate weights depend on the objective of the analysis and the relative importance assigned to each constraint.

### 4.2. Number of SA runs

Specifies the number of independent simulated annealing runs.
Because SA is a stochastic optimization algorithm, different runs can produce different rupture-rate solutions. A sufficiently large number of independent runs is recommended to characterize the variability among solutions; for example, more than 100 runs may be appropriate.
To record a random seed for each run, open getmodelparameters.m and uncomment the following line (currently around line 141):

params.seed = randi(100000);

This allows the generated seed to be recorded and used to reproduce a run, provided that the relevant random-number handling and execution conditions remain unchanged.

### 4.3. Number of iterations

Specifies the number of iterations performed by the simulated annealing algorithm.
The required number depends on the size and complexity of the problem. As an initial guideline, approximately 2,000 iterations per candidate rupture can be considered. This is a starting point rather than a universal requirement, and convergence and solution variability should be assessed for the specific application.

## 5. Outputs and additional tools

NEXTQUAKE automatically saves the generated figures and a log file containing information about each run in the output folder.

Additional tools are available for producing supplementary results:

plotLocalMFD(bestRates, invData, targets, subsectionData, params)
exportInputToOpenQuake(bestRates, outputs, invData, All_Rupture_Candidates, faultResults, params)
save_all_ruptures_SD_Batched(All_Final_Ruptures, subsectionData)

These functions can be used to:

plotLocalMFD: plot the magnitude–frequency distribution of individual fault sections.
exportInputToOpenQuake: export the optimized rupture rates as OpenQuake input files for subsequent PSHA calculations. Both Poissonian and BPT-based inputs can be exported.
save_all_ruptures_SD_Batched: save the rupture catalogue and associated subsection information in batches.

To generate fault-specific MFD plots, export OpenQuake input files, or save additional plots showing the candidate ruptures, uncomment the corresponding function calls in NEXTQUAKE.m.

**OpenQuake compatibility:** The export functionality has been tested with OpenQuake version (TBD).

# Contributors
Alessandro Valentini - University of Vienna, Department of Geology

# Licensing and citation
If you use this code, please cite our paper (...) and the Zenodo DOI (https://doi.org/10.5281/zenodo.23160969).
This project is licensed under the MIT License - see the LICENSE file for details.

# Contact
alessandro.valentini@univie.ac.at

# Acknowledgments
This work was carried out within the framework of the Austrian Science Fund (FWF) project NEXTQUAKE (Project no.:  PAT  2160424, Grant-DOI:10.55776/PAT2160424). 
