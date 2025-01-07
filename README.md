# ensemble-Species-Distribution-Modeling
Forecasting hotspots of climatic suitability for grassland restoration under climate change in North America
Santosh Kumar Rana1*, Jessica Lindstrom2, Melissa A. Lehrer1, Marissa Ahlering3, Jill Hamilton1#

1 Department of Ecosystem Science and Management, Pennsylvania State University, State College, PA 16802, USA
2 Department of Biological Sciences, North Dakota State University, Fargo, ND, 58102, USA
3 The Nature Conservancy, 1101 West River Parkway, Suite 200, Minneapolis, MN, 55415, USA


#Corresponding author: jvh6349@psu.edu (Jill Hamilton)
*Present address: Arkansas Biosciences Institute, Arkansas State University, Jonesboro, AR 72401, USA

Abstract 
Local species-climate relationships are often considered in restoration management. However, as climate change disrupts species-climate relationships, identifying factors that influence the climatic niche now and into the future for individual species, functional groups, and communities will be increasingly important for restoration. This involves identifying hotspots of climatic suitability to target restoration efforts.
We identified bioclimatic factors influencing the distribution of 26 species and their associated functional groups commonly used in grassland restoration in North America using ensemble species distribution modeling (SDM). We predicted their climatic niche under current and future (2050) climates and identified hotspots of climatic suitability for diverse species and functional groups. These hotspots were then overlaid with estimates of landscape connectivity and protected status to quantify potential suitability for restoration now and into the future.  
Temperature and precipitation during warmer quarters largely influenced grassland species’ distribution. Climatically suitable hotspots were identified in Minnesota, North Dakota, and South Dakota, with projected northward shifts under future climate scenarios. Overlaying these hotspots with estimates of landscape connectivity and protected status revealed limited connectivity and protection, highlighting regions to prioritize for restoration and conservation efforts.
Leveraging an understanding of species-climate relationships, this research emphasizes the importance of quantifying connectivity and protected status across aggregated hotspots suitable climate for restoration and conservation. Identifying these hotspots now and into the future can be used to prioritize restoration efforts, ensuring long-term maintenance of functional ecosystems across grassland communities.

Keywords: Biomod2, climate-smart, grassland prairies, range shift, resilient & connected network, restoration 



This repository contains R scripts and data files for:
1. **Spatial Thinning for Rarefaction of Occurrence Points** using the `spThin` R package.
2. **Variance Inflation Factor (VIF) Analysis** to address multicollinearity among bioclimatic variables.
3. **Ensemble Species Distribution Modeling** for 26 species.

## Repository Structure
- `scripts/`: R scripts for spatial thinning, VIF analysis, and SDM modeling.
- `data/`: Data files, including bioclimatic variable values and input CSV for SDM.

## How to Use
1. Clone the repository:
   ```bash
   git clone https://github.com/your_username/Spatial-Thinning-and-SDM.git

