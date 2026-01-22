
#### Figures

The code above (`fig2/fig2.ipynb`, `fig3/fig3.ipynb`, `fig4/fig4.ipynb`) generates panels for Figures 2–4 in the main text, as well as the associated supplementary figures (`fig2/figs`, `fig3/figs`, `fig4/figs`). These figures were subsequently assembled and annotated in Adobe Illustrator; stylistic elements (e.g., font, line width, annotations) were adjusted for clarity and consistency, but the underlying data and graphical content were not altered. The file `utils.py` contains helper functions used for both data analysis and plotting. Python library dependencies are listed in `requirements.txt`.  

#### Bootstrap resampling for confidence intervals

To estimate uncertainty in our _in vivo_ metrics (relative tumor burden and seeding, number of tumors, and clonal tumor sizes), we implemented a two-level nonparametric bootstrap procedure. First, mice were resampled with replacement 100 times to generate bootstrap cohorts, weighting each mouse by the number of neoplasms above a size threshold so that mice with higher burden contributed proportionally more (except for the bootstrap of tumor counts, where all mice were weighted equally). The size threshold was set at 100 cells for three-week cohorts and 1 cell for two-day cohorts. Within each resampled cohort, we then resampled tumors 100 times to account for clonal-level variability. For each bootstrap replicate, we computed the relevant statistics, and 95% confidence intervals were taken as the 2.5th and 97.5th percentiles of the bootstrap distributions. In the figures above, each plotted point typically corresponds to the mean of the data, and error bars indicate the bootstrap confidence intervals.

See `utils.bootstrap_burden_n_sizes()` and `bootstrapping.py` for details.
