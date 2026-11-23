## Figures 2-4

The code above generates panels for Figures 2–4 in the main text, as well as the associated supplementary figures (see `fig2`, `fig3`, `fig4`). These figures were subsequently assembled and annotated in Adobe Illustrator; stylistic elements (e.g., font, line width, annotations) were adjusted for clarity and consistency, but the underlying data and graphical content were not altered.  

The file `utils.py` contains helper functions used for both data analysis and plotting.

**Requirements:** `python 3.12.2` or above (tested with `python 3.12.16`), `numpy 2.4.1` or above (tested with `numpy 2.5.5`), `pandas 3.0.0` or above (tested with `pandas 3.1.1`), `scipy 1.17.0` or above (tested with `scipy 1.18.1`), `matplotlib 3.10.8` or above (tested with `matplotlib 3.11.2`), `seaborn 0.13.2`.


#### Bootstrap resampling for confidence intervals

To estimate uncertainty in our _in vivo_ metrics (relative tumor burden and seeding, number of tumors, and clonal tumor sizes), we implemented a two-level nonparametric bootstrap procedure. First, mice were resampled with replacement 100 times to generate bootstrap cohorts, weighting each mouse by the number of neoplasms above a size threshold so that mice with higher burden contributed proportionally more (except for the bootstrap of tumor counts, where all mice were weighted equally). The size threshold was set at 100 cells for three-week cohorts and 1 cell for two-day cohorts. Within each resampled cohort, we then resampled tumors 100 times to account for clonal-level variability. For each bootstrap replicate, we computed the relevant statistics, and 95% confidence intervals were taken as the 2.5th and 97.5th percentiles of the bootstrap distributions. In the figures above, each plotted point typically corresponds to the mean of the data, and error bars indicate the bootstrap confidence intervals.

See `utils.bootstrap_burden_n_sizes()` and `bootstrapping.ipynb` for details.
