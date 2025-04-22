- Codes for estimating connected and disconnected reproduction numbers are in the folder `uk_spatial`.
- `uk_spatial/simulated_data` contains all the codes related to the simulated data.
- `uk_spatial/region_fitting` contains all the codes related to the case study: Regions of England and LTLAs of the North East region.

- *Fig. 2–3:*
  - `simulated_data/simulated_data_generation.R` – for generating simulated data
  - `simulated_data/simulated_fitting_region.R` – for estimating connected and disconnected reproduction numbers $R_t$ with simulated data 
  - `simulated_data/region_in_out_inf.R` – for estimating infections from three different sources

- *Fig. 4:*
  - `region_fitting/fitting_region.R` – for estimating connected and disconnected $R_t$ over the 9 regions of England

- *Fig. 5–7:*
  - `region_fitting/north_east_ltla_fitting.R` – for estimating connected and disconnected $R_t$ over the 12 LTLAs of the North East region of England
  - `region_fitting/ltla_ne_in_out.R` – for estimating infections from three different sources for the LTLAs of the North East region
