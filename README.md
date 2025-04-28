# Enhancing transferability of species distribution models – guidance through algorithm, variable, and ensemble modelling decisions

### Eugene Yu Hin Yau, Alice C Hughes, Timothy C Bonebrake
[![](https://img.shields.io/badge/Citation-Scientific%20Data-blue)](https://doi.org/link)
[![](https://img.shields.io/badge/Archive-figshare/25037645-blue)](https://doi.org/10.6084/m9.figshare.25037645)
[![](https://img.shields.io/badge/License-CC%20BY%204.0-blue)](https://creativecommons.org/licenses/by/4.0/)


Species distribution models (SDMs) are useful tools to estimate species’ distributions in unsampled environments. SDMs not only fill the persisting knowledge gap of species’ current distribution (Wallacean shortfall), but also allow modellers to estimate species distributions in novel environments (e.g., future climates). However, data available often fail to represent the full range of viable environments (fundamental niche) due to sampling biases (e.g., inaccessible regions), distribution delimited by dispersal barriers (e.g., ocean), and the availability of potentially suitable habitats. Therefore, many SDMs are likely based on truncated niche spaces and will have to extrapolate in unsampled environments, ultimately over- or under-estimating species’ distributions. 
We simulated 1,567 virtual species. Four occurrence datasets, each subsampled from different latitudinal range to manually truncate the sampled niche, were sampled from each of these virtual species. These datasets were used to explore the consequences of niche truncation on SDM transferability. We also tested SDM transferability under various combinations of SDM settings, including SDM algorithms, variable selection strategies, ensemble algorithms, and single model TSS filtering thresholds.
We found that niche truncation caused by 21.7% areas unsampled increased mean SDM inaccuracy by 41.5%, while the sampled niche can be truncated even when sampling from the full latitudinal range. SDMs tend to overpredict species distributions when extrapolating in unsampled environments, distorting our estimation of climate change response and invasion risk. We suggest using a set of SDM algorithms, tailored for modelling goals to minimise either false positives or false negatives, when extrapolating SDMs. We also recommend an inclusive approach in variable selection when ecological knowledge is insufficient to precisely identify the most relevant variables. Moreover, the choice of ensemble algorithms significantly impacts SDM transferability, thus, variation between ensemble model outputs should be reported. While niche truncation is best dealt with using data more representative of species’ full range of viable environments, decisions related to algorithm, variable and ensemble choices can improve model transferability.
 

<img align="left" src="https://github.com/eugeneyau/Enhancing-SDM-transferability/blob/main/readme_images/trunc%20range%20diff.png" width=900>    

<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />





# Table of Contents

- Cleaned butterfly occurrence records from [GBIF](https://doi.org/10.15468/dd.nvw5wr), [Borneo Butterfly Distribution Database (B2D2)](https://www-users.york.ac.uk/~jkh6/), [Spatial occurrence data for the animals of Bangladesh derived from Facebook](https://doi.pangaea.de/10.1594/PANGAEA.948104), and other credible published sources. 
  - [`Tropical Asian Butterfly Occurrence Database`](https://doi.org/10.6084/m9.figshare.25037645)
- Distribution maps of tropical Asian butterflies as predicted by species distribution models(SDMs) can be downloaded as separate raster files or one single PDF file from our [`Figshare repository`](https://doi.org/10.6084/m9.figshare.25037645).
- R script used to construct SDMs:
  - [`Code/SDM`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/tree/main/Code/SDM)
     - [`R Markdown file`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/blob/main/Code/SDM/RMD_TropicalAsia_Bfy_SDM.Rmd)
     - [`R script`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/blob/main/Code/SDM/TropicalAsia_Bfy_SDM.R)
- Please download files essential for running our R scripts from our [`Figshare repository (SDMsupp_files.zip)`](https://doi.org/10.6084/m9.figshare.25037645) before running our codes.
- Additional R scripts used to clean our dataset, prepare for distribution modeling, analyze SDM outputs, and buffer occurrence points of unmodelled species:
  - [`Code/Supplementary`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/tree/main/Code/Supplementary)
     - [`Harmonize species names`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/blob/main/Code/Supplementary/Bfy_Data_supp_update_sp_name.R)
     - [`Clean dataset family names`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/blob/main/Code/Supplementary/Bfy_Data_supp_update_family_name.R)
     - [`Identify possible biogeographic range of dispersal for each species`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/blob/main/Code/Supplementary/Bfy_Data_supp_id_landmass_mask.R)
     - [`Calculate spatial sampling effort`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/blob/main/Code/Supplementary/Bfy_Data_supp_get_dens_ras.R)
     - [`Map diversity patterns by stacking single species distributions`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/blob/main/Code/Supplementary/Bfy_Data_supp_plot_alpha_diversity.R)
     - [`SDM performance evaluation`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/blob/main/Code/Supplementary/Bfy_Data_supp_eval_summary%20(PO).R)
     - [`SDM variable importance`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/blob/main/Code/Supplementary/Bfy_Data_supp_var_imp_analysis.R)
     - [`Buffer occurrence points of unmodelled species`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/blob/main/Code/Supplementary/Bfy_Data_supp_unmodelled%20species_point_richness.R)
     - [`Compare SDM results with that of Daru (2024)`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/blob/main/Code/Supplementary/Bfy_Data_supp_validate_daru2024_map.R)
- JavaScript code used in Google Earth Engine to extract and filter Landsat data for use as SDM variable (NDVImean), view in [`Google Earth Engine`](https://code.earthengine.google.com/7e1c649f06f22536419886e34a14d830) or download code from here:[`Code/Variables`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/blob/main/Code/Variables/GEE_NDVImean.txt)

