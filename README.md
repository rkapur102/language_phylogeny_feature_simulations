# language_phylogeny_feature_simulations

Code and data for research on modeling language evolution and feature change dynamics (diffusion, stability, internal change) through computer simulations.

A NOTE ON THE DATA:

The EURASIA data in the "populated_places" folder is from the Github repository associated with Wichmann (2017) (https://github.com/Sokiwi/LanguageFamilyExpansions) but originates from the GeoNames database (https://www.geonames.org/). The data is in the form of a zip file for size reasons and must be expanded to obtain the txt file. Then, all file paths in the simulation scripts must be filled in accordance with one's system. 

The Glottolog data used in calculating the real family dispersion distribution is also stored in this manner (the CSVs are compressed into a zip file) and also must be expanded for use, with file paths modified accordingly.

GUIDE TO OTHER FILES/FOLDERS:

"data_analysis" folder
----------------------
data_analysis.R - main data analysis and figure generation file. Contains the functions used to generate graphs for descriptive statistics. Also contains functions for plotting simulated languages, largest simulated families, language feature values, or populated places on a map

three_largest_families_visualization.R - code to plot three largest simulated families on map.

power_law_graphs.R - generates graphs relating to power law distribution.


"dispersion" folder
-------------------
real_language_dispersion.R - script used to calculate real language family dispersions from Glottolog data and calculate linear relationship with family size.

simulation_dispersion.R - script used to calculate dispersion of simulated families.  


"simulations" folder
--------------------
These are various types of simulations.

migration_only.R - Simulation used solely to examine birth, death, and migration.

structured.R - Simulation used to examine structured features. No descriptive statistics capabilities built in.

unstructured.R - Simulation used to examine unstructured features. No descriptive statistics capabilities built in.

unstructured_with_ds.R - Simulation used to examine unstructured features. Functions used to calculate descriptive statistics are included and used during the simulation.
