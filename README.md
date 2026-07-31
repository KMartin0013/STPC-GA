# STPC-GA

**Spatio-Temporal Principal Component Analysis for Gravimetry Applications**

STPC-GA is an open-source MATLAB toolbox for regional post-processing of
GRACE and GRACE Follow-On (GRACE-FO) Level-2 spherical-harmonic coefficient
(SHC) products. It provides a reproducible workflow for regional signal
localization, multi-product gap filling and ensemble reconstruction, and
empirical assessment of region-specific noise variability.

The workflow combines three sequential components:

1. **Spherical Slepian localization** to represent the gravity field using
   basis functions concentrated within a user-defined region.
2. **Multichannel Singular Spectrum Analysis (MSSA)** to fill data gaps and
   extract temporally coherent variability shared by multiple GRACE/GRACE-FO
   products.
3. **STPC screening and reconstruction** to separate retained signal
   reconstructed components (RCs) from removed noise RCs at a selected
   statistical significance level.

STPC-GA is intended as a regional GRACE/GRACE-FO processing toolbox. It
reduces leakage and excessive smoothing within the original SHC bandwidth;
it does **not** increase the intrinsic spatial resolution of the input
GRACE/GRACE-FO products.

---

## Contents

- [Main capabilities](#main-capabilities)
- [Interpretation of noise and significance levels](#interpretation-of-noise-and-significance-levels)
- [Requirements](#requirements)
- [Installation and startup](#installation-and-startup)
- [Input data and directory layout](#input-data-and-directory-layout)
- [Quick start](#quick-start)
- [R/S/M/N configuration](#rsmn-configuration)
- [Results and cache isolation](#results-and-cache-isolation)
- [Computational design](#computational-design)
- [Outputs](#outputs)
- [Known limitations](#known-limitations)
- [Citation](#citation)
- [License](#license)

---

## Main capabilities

- Regional conversion of GRACE/GRACE-FO SHCs to spherical Slepian
  coefficients.
- Application-dependent selection of the study-region buffer `R`:
  - positive/outward buffers for leakage-out-dominated terrestrial or island
    applications;
  - negative/inward buffers for leakage-in-dominated oceanic sub-basins;
  - simulation-based selection for regions with strong land-ocean or
    ice-ocean contrasts.
- Automatic or user-defined selection of:
  - Slepian truncation number `S`;
  - MSSA sliding-window length `M`;
  - retained MSSA mode number `N`.
- MSSA gap filling within and between GRACE and GRACE-FO periods.
- Multi-product reconstruction using GRACE/GRACE-FO products as separate
  MSSA channels.
- Region-specific empirical extraction of noise components at
  user-defined screening levels.
- Optional unfiltered, Gaussian-filtered, and DDK-filtered comparison
  products.
- Hydrological, cryospheric, and oceanic output units, including:
  - equivalent water height (EWH);
  - total mass change;
  - mass-term sea-level change.
- Configuration-isolated caches for reproducible parameter experiments.
- Cold-start and warm-start runtime and cache-status logs in
  `Benchmark_Log.csv` and `Benchmark_Log.mat`.
- Example workflows for the Mekong River Basin, Yangtze River Basin,
  Greenland Ice Sheet, and South China Sea.

---

## Interpretation of noise and significance levels

This section is important for interpreting STPC-GA outputs correctly.

### 1. What STPC-GA calls “noise”

For each spherical Slepian coefficient, MSSA produces a set of reconstructed
components. Let

- $X_{\mathrm{MSSA}}$ denote the reconstruction obtained from the available
  MSSA components before STPC screening; and
- $X_{\mathrm{STPC},p}$ denote the reconstruction obtained from the
  signal components retained by the STPC screening rule at significance
  level $p$.

The corresponding removed component is

$$
X_{\mathrm{noise},p}
=
X_{\mathrm{MSSA}}
-
X_{\mathrm{STPC},p}.
$$

In the software and output files, this difference is reported as the
**region-specific noise estimate** or **noise component**.

It should not automatically be interpreted as:

- pure white noise;
- a complete GRACE formal-error covariance;
- only north-south striping;
- proof that every removed component is non-geophysical.

The removed component can contain a mixture of measurement noise,
processing-center-specific variability, residual striping, leakage-related
variability, and temporally correlated or quasi-periodic errors.

### 2. What the K-S and Lilliefors tests do

STPC-GA uses Kolmogorov-Smirnov and Lilliefors tests as a distribution-based
screening rule for MSSA-derived RCs. These tests examine distributional
characteristics; they do **not** test temporal whiteness or independence.

Consequently, colored or quasi-periodic GRACE errors—such as residual
ocean-tide or atmosphere-ocean de-aliasing errors—may be retained together
with real geophysical variability. This is particularly relevant in coastal
and low-signal oceanic regions.

### 3. Meaning of the significance level `p`

The user-selected significance level `p` controls the statistical screening
of RCs and therefore changes the retained reconstruction and removed
noise-like component.

Users should compare several values of `p` and examine the stability of the
result. The manuscript examples use 5%, 10%, and 30% screening levels.

### 4. Noise thresholds for extreme-event analysis

For optional flood and drought identification, a basin-average,
detrended-and-deseasonalized TWS series can be compared with thresholds
derived from the standard deviation of the STPC noise-like time series, for
example $3\sigma$.

This threshold is an empirical, region-specific detection criterion. Event
interpretation should be checked against independent observations such as
land-surface-model TWS, river discharge, water level, precipitation, or
surface-water extent whenever such data are available.

### 5. Multi-product ensembles

Adding more GRACE products does not necessarily improve the result.
Performance also depends on:

- independence among products;
- data-center processing strategies;
- regional error characteristics;
- low-degree replacement and background-model consistency.

Combination products such as COST-G may contain information already
represented by individual analysis-center products. Users should therefore
test alternative product combinations rather than treating all available
products as independent, equally informative channels.

---

## Requirements

STPC-GA has been tested with **MATLAB R2022b**.

### MATLAB toolboxes

- **Statistics and Machine Learning Toolbox** — required for the statistical
  tests.
- **Parallel Computing Toolbox** — strongly recommended for gap filling and
  comparison workflows.

### Required Slepian dependencies

Download:

- [slepian_alpha](https://github.com/csdms-contrib/slepian_alpha)
- [slepian_bravo](https://github.com/csdms-contrib/slepian_bravo)
- [slepian_delta](https://github.com/csdms-contrib/slepian_delta)
- [slepian_zero](https://github.com/csdms-contrib/slepian_zero)

### Optional filtering dependency

Required only when DDK comparison products are requested:

- [GRACE-filter](https://github.com/strawpants/GRACE-filter)

### Optional plotting dependencies

Some diagnostic figures require:

- [m_map](https://www.eoas.ubc.ca/~rich/map.html)
- [hatchfill2](https://www.mathworks.com/matlabcentral/fileexchange/53593-hatchfill2)
- [special heatmap](https://www.mathworks.com/matlabcentral/fileexchange/125520-special-heatmap)

A compatible compact-subplot implementation is included, so
`tight_subplot` does not need to be installed separately.

---

## Installation and startup

Clone or download the repository and start it from MATLAB:

```matlab
projectRoot = 'D:\path\to\STPC-GA';
startup_STPC_GA(projectRoot);
```

If third-party dependencies are stored outside the repository:

```matlab
projectRoot = 'D:\path\to\STPC-GA';
dependencyRoot = 'D:\MATLAB\STPC_dependencies';
startup_STPC_GA(projectRoot, dependencyRoot);
```

The dependency root can also be specified using the
`STPC_DEPENDENCIES` environment variable.

The startup routine:

- adds STPC-GA and dependency folders to the MATLAB path;
- defines `IFILES` and `STPC_SRC`;
- creates missing writable runtime/cache directories.

Startup does **not** download:

- third-party software;
- GRACE/GRACE-FO SHCs;
- degree-1 or low-degree replacement coefficients;
- GIA models;
- regional boundary files.

---

## Input data and directory layout

The repository root is used as the default software and data root:

```text
STPC-GA/
├─ COASTS/                       regional boundary files
├─ Data/                         static supporting data
├─ EARTHMODELS/
├─ GIA/
├─ GRACE/
│  ├─ Degree1/
│  ├─ Degree2/
│  ├─ Originals/RL06/<product>/
│  └─ SlepianExpansions/         generated cache
├─ MOD/                          generated simulation models
├─ examples/
├─ external/                     user-installed dependencies
├─ src/
└─ startup_STPC_GA.m
```

Example `used_files.mat` files document expected inputs but do not replace
the original GRACE/GRACE-FO data.

The current release supports RL06 parsing workflows for the following
products:

- CSR;
- JPL;
- GFZ;
- ITSG;
- COST-G;
- Tongji-Grace2022.

Users remain responsible for checking:

- file naming and date conventions;
- maximum degree and order;
- degree-1, C20, and C30 replacement strategies;
- GIA corrections;
- atmosphere/ocean background-model compatibility;
- missing or incomplete months.

Temporal gravity-field products can be obtained from the
[International Centre for Global Earth Models (ICGEM)](https://icgem.gfz.de/sl/temporal).

---

## Quick start

After installing dependencies and input data, run one of the examples:

```matlab
run('examples/example_Mekong.m')
run('examples/example_Yangtze.m')
run('examples/example_greenland.m')
run('examples/example_SCS.m')
```

The examples locate the repository automatically and use:

```matlab
results = run_case_fast(config);
```

Use the compatibility implementation only when required:

```matlab
results = run_case(config);
```

A typical configuration is:

```matlab
config.R = "STPC";
config.S = "STPC";
config.M = "STPC";
config.N = "STPC";

config.turningNumber = 5;

config.resultRoot = fullfile(projectRoot, 'Results_STPCGA1_mekong');
config.figureRoot = fullfile(projectRoot, 'Figure_STPCGA1_mekong');
config.caseName = 'paper_default';

config.redo = false;
config.plotProcess = true;

results = run_case_fast(config);
```

Check the distributed example scripts for the exact field names supported by
the current release.

---

## R/S/M/N configuration

| Parameter | Meaning | Automatic mode | Other supported values |
|---|---|---|---|
| `config.R` | Regional buffer, in degrees | `"STPC"` | Numeric value |
| `config.S` | Number of retained Slepian functions | `"STPC"` | `"Shannon"`, `["Fixed","0.1"]`, `["Sel_Gau_default","500"]`, positive integer |
| `config.M` | MSSA reconstruction window | `"STPC"` | Positive integer |
| `config.N` | Number of retained MSSA modes | `"STPC"` | `"Wcorr"`, positive integer |

Example:

```matlab
config.R = "STPC";
config.S = ["Sel_Gau_default","500"];
config.M = 120;
config.N = "STPC";
config.turningNumber = 5;
```

Interpretation:

- `config.R = "STPC"` determines the regional buffer using the configured
  simulation procedure.
- `config.S = "Shannon"` uses the Shannon number.
- `config.S = ["Fixed","0.1"]` retains Slepian functions with concentration
  eigenvalues greater than 0.1.
- `config.S = ["Sel_Gau_default","500"]` applies selective 500-km Gaussian
  smoothing to higher-order Slepian components.
- Numeric `config.S` directly fixes the number of retained Slepian functions.
- `config.N = "Wcorr"` generates alternative `N` candidates using
  w-correlation; the final selection still follows the STPC
  candidate-screening workflow.

The candidate grid adapts to fixed parameters:

| `S` mode | `N` mode | Candidate grid |
|---|---|---|
| STPC | STPC or Wcorr | K × K |
| Fixed | STPC or Wcorr | 1 × K |
| STPC | Fixed | K × 1 |
| Fixed | Fixed | 1 × 1 |

When both `S` and `N` are fixed, the tuning search is bypassed.

`turningNumber` defaults to 5. Some legacy diagnostic layouts were designed
for five candidates, so custom values should be checked when
`plotProcess = true`.

---

## Results and cache isolation

Each run should define a shared result root and a unique `caseName`:

```matlab
config.resultRoot = fullfile(projectRoot, 'Results_STPCGA1_mekong');
config.figureRoot = fullfile(projectRoot, 'Figure_STPCGA1_mekong');
config.caseName = 'paper_default';
```

Configuration-dependent files are stored under:

```text
Results_STPCGA1_mekong/
└─ Configurations/
   └─ paper_default/
      ├─ Configuration.mat
      ├─ Basic_Information.mat
      ├─ Slepian_Information.mat
      ├─ MainSlep_*.mat
      ├─ MainMSSA_*.mat
      ├─ MainSTPC_V3_*.mat
      ├─ AddData/
      └─ Benchmark_Log.*
```

Internal names such as `MainSTPC_V3_*` are retained for cache-schema
compatibility. They are not the public software version number.

Changing a science-affecting parameter while reusing the same `caseName`
raises a configuration-mismatch error when `redo = false`.

Use either:

- a new `caseName`; or
- `redo = true` to intentionally recompute and replace the cached case.

Automatic-buffer simulations are stored under `<resultRoot>/Simulate` and
are reused only when the stored region and model metadata match the current
configuration.

---

## Computational design

STPC-GA avoids repeating the most expensive decompositions for every
candidate pair.

The workflow:

1. computes the required Slepian basis once for a region and configuration;
2. performs the maximum required MSSA decomposition once;
3. stores p-independent SVD and RC products;
4. evaluates candidate `S/N` combinations by selecting and summing existing
   components;
5. produces complete maps and diagnostic figures only after final parameter
   selection.

Warm runs reuse:

- regional Slepian bases;
- transformed Slepian coefficients;
- gap-filled series;
- MSSA decompositions;
- final STPC products;
- compatible simulation results.

Runtime logs separate:

- Slepian transformation;
- Gaussian/DDK comparison processing;
- main MSSA processing;
- STPC reconstruction.

Cold-run time, memory use, and disk requirements depend strongly on:

- SHC maximum degree;
- region size and geometry;
- number of products;
- number of parameter candidates;
- requested figures;
- enabled comparison filters.

---

## Outputs

Depending on the region type and enabled options, STPC-GA produces:

1. Gap-filled regional time series averaged across selected
   GRACE/GRACE-FO products.
2. Gridded regional reconstructions of:
   - equivalent water height;
   - total mass change;
   - mass-term sea-level change.
3. Signal-like STPC reconstructions at selected significance levels.
4. Region-specific MSSA-to-STPC noise-like residuals.
5. Optional unfiltered, Gaussian, and DDK comparison products.
6. Buffer, `S`, `M`, and `N` selection diagnostics.
7. Spatial and temporal eigenvalue diagnostics.
8. Cold/warm runtime and cache-status logs. Peak memory should be measured
   separately with the supplied benchmark workflow when needed.
9. Optional trend, seasonal, and extreme-event analysis products, depending
   on the selected example or analysis driver.

---

## Known limitations

- STPC-GA reduces leakage and smoothing but cannot exceed the intrinsic
  spatial resolution of the input SHCs.
- K-S and Lilliefors tests do not diagnose temporal whiteness.
- Colored and quasi-periodic errors may remain in the retained signal,
  especially in coastal and low-signal regions.
- The MSSA-to-STPC residual is an empirical noise-like estimate, not a full
  formal uncertainty covariance.
- Strongly dependent input products can unintentionally receive excess
  weight in a multi-product ensemble.
- Automatic turning points may be less distinct in weak-signal regions.
  Users should inspect diagnostic figures and compare automatic and fixed
  parameter choices.
- A fixed `N` must not exceed the available RC count. STPC-GA raises an
  error rather than silently truncating it.
- DDK processing requires both:
  - `GRACE-filter-master/src/matlab`;
  - `GRACE-filter-master/data/DDK`.
- COST-G and Tongji parsing requires users to verify product-specific
  metadata and correction strategies.
- Cold runs may require several gigabytes of memory and substantial disk
  space.
- `Benchmark_Log.*` records runtime and cache status, not automatic
  peak-memory measurements.
- The public release excludes manuscript-only plotting scripts, local
  validation drivers, generated caches, and thesis-only routines.

---

## Citation

Please cite the STPC-GA software paper when available, together with the
foundational methods relevant to the selected workflow.

1. Harig, C., & Simons, F. J. (2012). Mapping Greenland's mass loss in
   space and time. *Proceedings of the National Academy of Sciences*,
   109(49), 19934–19937.  
   https://doi.org/10.1073/pnas.1206785109

2. Ma, Z., Fok, H. S., Tenzer, R., & Chen, J. (2024). A novel Slepian
   approach for determining mass-term sea level from GRACE over the South
   China Sea. *International Journal of Applied Earth Observation and
   Geoinformation*, 132, 104065.  
   https://doi.org/10.1016/j.jag.2024.104065

3. Gauer, L. M., Chanard, K., & Fleitout, L. (2023). Data-driven gap
   filling and spatio-temporal filtering of the GRACE and GRACE-FO records.
   *Journal of Geophysical Research: Solid Earth*, 128(5), e2022JB025561.  
   https://doi.org/10.1029/2022JB025561

---

## License

STPC-GA is distributed under the repository [LICENSE](LICENSE). Third-party
dependencies retain their own licenses and are not redistributed by this
repository.
