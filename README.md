# STPC_GA (Version 2)
A Lite version of Spatio-Temporal Principal Component (STPC) filter for basin-scale hydrological, cryospheric and oceanographic applications of GRACE/GRACE-FO.  

By reconciling two eigenvalue problems in spherical Slepian method and Multichannel Singular Spectrum Analysis (MSSA) based on different noise significance levels, the STPC filter establishes a statistical framework for basin-scale analyses via automatically determined, region-specific physical parameters governing spatio-temporal localization, leakage-in/out effects differentiation, uncertainty minimization across data centers, and mitigation of north-south striping.

## Notice before use
1. The Slepian procedure of this code is based on the [foundational Software from the Simons Laboratories](https://geoweb.princeton.edu/people/simons/software.html), and [its add-on version](https://github.com/KMartin0013/Slepian_ocean_add-on) for ocean application (e.g., GAD, GIA, IB correction).
2. It is strongly recommended to familiarize yourself with the foundational work by <a href="https://polarice.geo.arizona.edu/">C. Harig</a> &amp; <a href="http://www.frederik.net">F. J. Simons</a> ([Ref 1](https://doi.org/10.1073/pnas.1206785109)), and our previous work on ocean ([Ref 2](https://doi.org/10.1016/j.jag.2024.104065)) before using this add-on code.
3. **As the code has not yet been fully optimised, it is advisable to reserve 5–8GB of space (for ocean cases maybe double the space) before running each case.** 


**Required software:**<br>
[slepian_alpha](https://github.com/csdms-contrib/slepian_alpha) (If a higher version of MATLAB, you must replace ls2cell.m with our version in the /src directory.)  
[slepian_bravo](https://github.com/csdms-contrib/slepian_bravo)  
[slepian_delta](https://github.com/csdms-contrib/slepian_delta)    
[slepian_zero](https://github.com/csdms-contrib/slepian_zero) (actually only guyotphysics.m)  
[GRACE-filter](https://github.com/strawpants/GRACE-filter) (for DDK filter)  

**Required software for plotting:**<br>
[m_map](https://www.eoas.ubc.ca/~rich/map.html)  
[tight_subplot](https://ch.mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w)  
[hatchfill2](https://ch.mathworks.com/matlabcentral/fileexchange/53593-hatchfill2)  
[special heatmap](https://ch.mathworks.com/matlabcentral/fileexchange/125520-special-heatmap/)  


## Output of this package
1. Area-weighted equivalent water height (EWH, unit: cm), mass changes (unit: Gt), mass-term sea level changes (unit: cm, in ocean cases) time series at different significance levels of noise.  
2. Gridded EWH, mass-term sea level changes (MSL, in ocean cases) at different significance levels of noise.
3. Noise level specific to the given regions. This can be obtained by subtracting the final STPC procedure results from the MSSA procedure results.  
4. Output 1 and 2 but for other filters (e.g., Gaussian filter, DDK filter, no filter).  

Notice: All output time series will be gap-filled ([Ref 3](https://doi.org/10.1029/2022JB025561)) and averaged across different data centers specified in the running code.  


## Usage of this package  
1. Download the required software (or codes) and this codes as mentioned above.  
2. Download the Spherical Harmonic coefficients and the associated correction files (e.g., degree 1, degree 20 correction) from different data centers.  
  2.1 Notice that for now this version supports the Release 05 (RL05) and Release 06 (RL06) data from the Center of Space Research (CSR) at the University of Texas, the Jet Propulsion Laboratory (JPL), GeoForschungsZentrum (GFZ), and the Institute of Geodesy at the Graz University of Technology (ITSG). If you use data centers other than these four, please revise the function 'grace2plmt_m.m' and 'gracedeg1_m' accordingly.  
  2.2 All necessary files for each data center is listed in '/GRACE/Originals/RL[05,06]/[CSR,GFZ,JPL,ITSG]0324/used_files.mat'.  
3. Three primary procedures are included in the STPC filter.  
  3.1 The first procedure 'run_slepian_main.m' calculates spherical Slepian functions and corresponding spherical Slepian coefficients, primarily following the guideline in [slepian_delta](https://github.com/csdms-contrib/slepian_delta). If not specified, the optimal buffer zone for the region R will be automatically determined.  
  3.2 The second procedure 'run_MSSA_main.m' fills the data gaps within and between GRACE and GRACE-FO datasets. The RC corresponding to different orders of spherical Slepian coefficients will be averaged over all useful products. If not specified, the optimal size of the sliding window M in the reconstruction procedure will be automatically determined.  
  3.3 The third procedure 'run_STPC_main.m' calculats spatial and temporal eigenvalues. The optimal combination of S and N would be determined for each given significance level p. Finally, the spatial maps of the mass variations and the area-weighted time series would be calculated based on different significance level p.  
  3.4 (optional) The fourth procedure 'run_mssa_smooth.m' calcualtes the averaged data applying other filters (e.g., Gaussian filter and DDK filter). No filter data can also be provided.  
4. Four cases have been provided in /examples.  
  'example_greenland.m' for a case study in the Greenland Ice Sheet.  
  'example_SCS.m' for a case study in the South China Sea.  
  'example_Yangtze.m' for a case study in the Yangtze River Basin.  
  'example_Mekong.m' for a case study in the Mekong River Basin.  

## Citation Information:
Please cite our work and the foundational work by C. Harig, F. J. Simons, and L. Gauer as appropriate if you find this package useful.  

[Ref 1] Harig, Christopher and Frederik J. Simons. 
Mapping Greenland's mass loss in space and time.
<i>Proc. Natl. Acad. Sc.</i>, 109(49), 19934-19937, 2012.
doi:https://doi.org/10.1073/pnas.1206785109

[Ref 2] Zhongtian, Ma, Hok Sum Fok, Robert Tenzer, Jianli Chen. 
A Novel Slepian Approach for Determining Mass-term Sea Level from GRACE over the South China Sea. <i>Int. J. Appl. Earth. Obs. Geoinf.</i>, 132, 104065, 2024.
doi:https://doi.org/10.1016/j.jag.2024.104065

[Ref 3] L. Gauer, K. Chanard, L. Fleitout. 
Data‐driven gap filling and spatio‐temporal filtering of the GRACE and GRACE‐FO records. <i>J. Geophys. Res. Solid Earth</i>, 128(5), e2022JB025561, 2023. doi:https://doi.org/10.1029/2022JB025561
