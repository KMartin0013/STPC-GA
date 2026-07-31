# STPC-GA (Version 1) changes

## Overview

Compared with the previous public release, **STPC-GA (Version 1)** adds a
portable startup workflow, user-selectable `R/S/M/N` parameters,
configuration-isolated caches, accelerated MSSA processing, reusable STPC
decompositions, optional w-correlation candidates, expanded GRACE product
handling, and structured runtime logging.

The principal user-visible changes are:

- `R`, `S`, `M`, and `N` can be selected automatically or specified
  directly in each launch script.
- Shannon truncation, eigenvalue-threshold truncation, selective Gaussian
  smoothing, and fixed numeric S are available through one interface.
- Fixed S or N automatically reduces the tuning grid from K×K to 1×K,
  K×1, or 1×1.
- Each `caseName` has an independent result and figure directory, with a
  saved configuration record that prevents incompatible cache reuse.
- MSSA gap filling and comparison filters have accelerated processing
  paths, while STPC tuning reuses a single maximum decomposition.
- W-correlation can provide alternative N candidates without replacing the
  existing STPC final-selection criterion.
- Startup now creates writable Slepian and transformation-cache directories
  and supports dependencies stored outside the source tree.
- Cold/warm stage timings and cache status are recorded in CSV and MAT logs.

## Suggested Git commit sequence

The following commit grouping keeps code review focused. The detailed
function tables below can be copied into the corresponding commit bodies.

1. `feat(startup): add portable STPC-GA bootstrap and runtime directories`
2. `feat(config): add isolated user-selectable R/S/M/N workflows`
3. `perf(mssa): add cached and parallel MSSA processing paths`
4. `feat(stpc): cache tuning decompositions and add w-correlation candidates`
5. `feat(slepian): extend buffer, truncation, and product handling`
6. `docs(release): publish STPC-GA Version 1 examples and documentation`

## Added functions

### Startup and configuration

| Function | Commit description |
|---|---|
| `startup_STPC_GA` | Add the public Version 1 entry point; locate the project, add source/dependency paths, initialize environment variables, and create runtime directories. |
| `initialize_stpc_directories` | Create writable Slepian, rotation, GRACE transformation, and simulation-cache directories without populating static inputs. |
| `stpc_dependency_root` | Resolve external dependencies from an explicit argument, `STPC_DEPENDENCIES`, `external/`, or the legacy `src/required_softwares` layout. |
| `stpc_source_root` | Resolve the active source directory without hard-coding a development folder name. |
| `v3_normalize_config` | Normalize the public `R/S/M/N` interface, translate supported legacy fields, validate modes, and derive configuration-specific result directories. |
| `v3_configuration_snapshot` | Extract science-affecting configuration fields into a stable cache-consistency record. |
| `v3_validate_configuration` | Prevent incompatible parameter sets from reusing one `caseName`; permit deliberate replacement only with `redo=true`. |
| `stpc_configure_center_list` | Normalize input-center lists, build deterministic ensemble identifiers, and maintain the ensemble registry used by downstream products. |
| `stpc_note_suffix` | Convert an optional experiment note into a safe, consistent filename suffix. |
| `stpc_concat_results` | Align heterogeneous STPC and smoothing-result struct fields before preserving the historical flat result-array interface. |
| `stpc_benchmark_append` | Append structured stage timing and cache-status records to CSV and MAT benchmark logs. |
| `run_case_fast` | Add the main controller that uses the accelerated smoothing-comparison path while retaining the complete Slepian–MSSA–STPC workflow. |
| `tight_subplot` | Provide a local compact-axes implementation with explicit argument validation, removing a required plotting download. |

### MSSA

| Function | Commit description |
|---|---|
| `MSSA_fast` | Add a compatible MSSA implementation with economy/truncated SVD options, selectable eigenvalue output, and unchanged EOF/PC/RC dimensions. |
| `MSSA_gap_fitting_only_correct_fast` | Add configurable fast iterative gap filling with optional figures/tests and reduced intermediate retention. |
| `collect_and_average_institutions_smooth_fast` | Add the fast multi-center smoothing-comparison pipeline with reusable per-method caches. |
| `run_mssa_smooth_fast` | Add the public wrapper for accelerated unfiltered, Gaussian, and DDK comparisons. |

### Slepian and automatic buffer selection

| Function | Commit description |
|---|---|
| `buffer_selection_diagnostics` | Add deterministic conflict resolution for PCC- and RMSE-selected buffers using trend, amplitude, and phase errors. |

### STPC and w-correlation

| Function | Commit description |
|---|---|
| `prepare_stpc_decomposition_cache` | Compute the maximum p-independent MSSA/SVD decomposition once and expose slices for K×K, 1×K, K×1, and 1×1 searches. |
| `compute_mssa_wcorrelation` | Compute joint-channel weighted correlations using SSA trajectory multiplicity weights. |
| `select_wcorr_N_candidates` | Rank distinct retained/residual RC boundaries from block-contrast w-correlation scores. |
| `det_N_WCorr` | Generate alternative `N` candidates from w-correlation while leaving final selection to STPC screening. |
| `plot_mssa_wcorrelation` | Plot compact w-correlation matrices and candidate-score diagnostics. |
| `run_stpc_center_specific_selected` | Reuse selected S/N values to export center-specific STPC products without repeating the full search. |
| `stpc_inter_center_sensitivity_m` | Summarize center-specific temporal and spatial spread for optional sensitivity assessment. |

## Modified functions

### Controllers and configuration preparation

| Function | Commit description |
|---|---|
| `setup_paths` | Replace fixed in-source dependencies with portable discovery, add dependency paths, initialize writable directories, preserve STPC-GA compatibility-function precedence, and report optional DDK availability. |
| `prepare_basic_info` | Persist shared/configuration result roots, `caseName`, normalized center metadata, and the resolved R/S/M/N selection. |
| `prepare_slepian_info` | Translate normalized R/S choices into automatic/fixed buffer and truncation settings, including Shannon, threshold, selective Gaussian, and numeric S. |
| `prepare_mssa_info` | Carry normalized M/N selections and configurable `turningNumber` into MSSA processing. |
| `prepare_stpc_info` | Pass normalized N-selection metadata and search settings to the STPC stage. |
| `run_case` | Normalize and validate configuration cases, isolate outputs, record per-stage cold/warm timing, reuse versioned STPC caches, and safely combine output structs. |

### MSSA processing

| Function | Commit description |
|---|---|
| `build_mssa_paths_and_tags` | Build MSSA text/figure/cache paths from configuration directories and normalized filename notes. |
| `collect_and_average_institutions` | Use normalized center metadata and configuration-local paths while retaining ensemble averaging behavior. |
| `collect_and_average_institutions_smooth` | Align the compatibility smoothing path with normalized center metadata and isolated outputs. |
| `det_opt_M` | Support either a fixed reconstruction window or STPC window selection, cache candidate calculations, and limit plotting to enabled process output. |
| `MSSA_final_noCDF_freqsort` | Respect process-plot settings and make component plotting robust to the actually available modes. |
| `MSSA_gap_fitting_correct` | Parallelize independent channel/institution gap-filling work and reduce repeated setup in iterative MSSA reconstruction. |
| `MSSA_gap_fitting_only_correct` | Align the compatibility gap-filling path with current plotting and output conventions. |
| `run_gapfilling` | Use the normalized gap parameters and updated parallel gap-filling workflow. |
| `run_mssa_main` | Support fixed or automatically selected `M`, carry resolved metadata, and preserve configuration-local caches. |

### Slepian transformation and simulation

| Function | Commit description |
|---|---|
| `run_slepian_main` | Apply normalized S-selection modes, preserve full-basis products for later STPC choices, and isolate configuration-dependent caches. |
| `choose_S_from_eigenvalues` | Implement the public selection order: Shannon, STPC, eigenvalue threshold, selective Gaussian, or direct fixed S; reject unsupported modes explicitly. |
| `det_opt_buffer` | Add shared-simulation metadata checks, automatic PCC/RMSE conflict diagnostics, configurable buffer candidates, and configuration-aware plotting/caching. |
| `MODEL_Make_m` | Write generated model products under the project data root while loading the static coefficient script from `src/MOD`. |
| `MODEL_Checkerboard_Make_new_m` | Move generated checkerboard models out of the source tree and retain the static source coefficient definition. |
| `TWSt2slept_model` | Store generated simulation models and Slepian expansions under writable project-root cache directories. |
| `grace2plmt_m` | Add COSTG/Tongji RL06 parsing, degree validation, date extraction, and documented atmosphere/ocean-background fallback handling. |
| `gracedeg1_m` | Remove an interactive pause, refresh metadata after marker creation, and make missing-file checks safe for unattended runs. |
| `plm_DDK` | Resolve GRACE-filter through the configured dependency root and report incomplete code/data installations explicitly. |
| `periodic_analysis_m` | Add OLS covariance, degrees-of-freedom correction, and propagated trend/amplitude/phase uncertainty outputs. |

### STPC selection and reconstruction

| Function | Commit description |
|---|---|
| `det_S_TP` | Make the number of S turning-point candidates configurable and validate insufficient candidate domains. |
| `det_N_TP` | Make N turning-point candidates configurable, support reduced fixed-parameter grids, and report insufficient modes clearly. |
| `run_MSSA_decompose` | Reuse cached full RC/SVD products, apply significance tests per requested p, and enforce fixed-N bounds. |
| `run_stpc_pvalue` | Convert candidate screening to lightweight RMS-only evaluation, suppress redundant plots, and generate full products only for the selected pair. |
| `run_stpc_main` | Coordinate fixed/adaptive S and N, optional w-correlation candidates, one-time decomposition caching, reduced search grids, and final-only reconstruction. |
| `grids_reconstruction_EWHorMSL_RMS` | Accept cached reconstruction inputs and avoid rebuilding invariant spatial products. |
| `grids_reconstruction_EWHorMSL_RMS_sigres` | Separate algorithm and plotting timing while producing full selected signal/residual reconstructions. |
| `plot_mssa_decompose` | Adapt component layouts to the resolved S/N counts and avoid assumptions tied to a full 5×5 search. |
| `plot_stpc_spatial_temporal_RMS` | Plot candidate grids using their actual dimensions, including single-row and single-column fixed-parameter cases. |

## Updated launch scripts

| Script | Commit description |
|---|---|
| `examples/example_Mekong.m` | Replace local absolute paths and legacy numeric modes with portable startup, public R/S/M/N options, and isolated Version 1 outputs. |
| `examples/example_greenland.m` | Add the portable Version 1 ice-region configuration using the shared public launch structure. |
| `examples/example_SCS.m` | Add the portable Version 1 ocean configuration with negative buffer candidates and isolated outputs. |
| `examples/example_Yangtze.m` | Add the portable Version 1 land-region configuration using the shared public launch structure. |

## Other release files

- `.gitignore` prevents generated bases, models, results, figures, external
  dependencies, and MATLAB temporary files from entering future commits.
- `external/README.md` documents the supported dependency layouts.
- `README.md` now documents Version 1 installation, startup, R/S/M/N
  configuration, cache isolation, computational design, outputs, and
  limitations.
- `src/MOD/Coeff.mat` updates the static simulation coefficients; generated
  models themselves are not shipped.

## Compatibility notes

- Supported legacy `buffer_deg`, `M='Auto'`, and `S_choice` inputs are
  translated only when the new public field is absent. New launch scripts
  should use `R/S/M/N`.
- Legacy `S_choice=6` remains an alias for STPC selection. Removed modes 4
  and 5 report an error instead of changing behavior silently.
- A numeric fixed N is applied consistently to all Slepian coefficients,
  including the ocean IB term, and must not exceed available RCs.
- Reusing a `caseName` with changed science parameters requires
  `redo=true`; otherwise the run stops before consuming incompatible caches.
