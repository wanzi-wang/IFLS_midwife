# thesis_data
This repo stores datasets and codes for my (master's) thesis

# Data
## Raw rainfall data Source:
- [University of Delaware Center for Climatic Research “Terrestrial Precipitation: 1900–2017 Gridded Monthly Time Series (1900–2017)](https://psl.noaa.gov/data/gridded/data.UDel_AirT_Precip.html)
- Monthly values for 1900/01 - 2017/12 V5.01
- Precipitation data information (OPeNDAP URL) is available on this [webpage](https://psl.noaa.gov/thredds/catalog/Datasets/udel.airt.precip/catalog.html?dataset=Datasets/udel.airt.precip/precip.mon.total.v501.nc)
- Each coordinate represents one IFLS district centroid.

## Processed Indonesian district coordinates data source:
- Maccini, Sharon, and Yang, Dean. Replication data for: Under the Weather: Health, Schooling, and Economic Consequences of Early-Life Rainfall. Nashville, TN: American Economic Association [publisher], 2009. Ann Arbor, MI: Inter-university Consortium for Political and Social Research [distributor], 2019-10-12. https://doi.org/10.3886/E113296V1

## Process flow (tentative):
1. Set up the coordinates dataset
2. Fetch rainfall data through the OPeNDAP URL and match it with district coordinates
3. Merge the coordinates and rainfall data

## Midwife, health, cognitive, educational outcomes data source: 
[Indonesia Family Life Survey (IFLS)](https://www.rand.org/well-being/social-and-behavioral-policy/data/FLS/IFLS.html)

Planning to use waves 1, 2, and 5.

Goal:
Identify the monthly rainfall deviations for each person during the estimated gestation period and the first year of life.

Calculate the long-run district-average-monthly rainfall.
