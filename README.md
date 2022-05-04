# Article: Effects of Data Aggregation(Buffer-Techniques)with Emerging Data on Bicycle Volume Estimation

## Authors 
#Md Mintu Miah1, 
#Stephen P Mattingly2, 
#Kate Kyung Hyun3, 
#Joseph Broach4

1 Doctoral Research Assistant, Department of Civil Engineering, the University of Texas at Arlington, TX 76019, USA, Corresponding Email: mdmintu.miah@mavs.uta.edu, ORCID: 0000-0001-6073-3896

2 Professor, Department of Civil Engineering, the University of Texas at Arlington, TX 76019, USA, Email: mattingly@uta.edu, ORCID: 0000-0001-6515-6813

3 Assistant Professor, Department of Civil Engineering, the University of Texas at Arlington, TX 76019, USA, Email: kate.hyun@uta.edu, ORCID: ORCID: 0000-0001-7432-8058

4 Adjunct Research Associate, Toulan School of Urban Studies and Planning, Portland State Univ., Portland, OR 97201. Email: jbroach@pdx.edu , ORCID: https://orcid.org/0000-0001-7753-501X. 

### Abstract

Data preparation/collection is the prerequisite of developing a bicycle direct demand model. Average Annual Daily Bicycle Traffic (AADBT) is the direct input for safety measures, crash analysis, infrastructure planning, and design. The majority of the past and current studies used the GIS tool to extract the variables at different buffer levels to develop the model which is time-consuming. However, very few researchers developed automatic data extraction scripts to extract the data and tested a wide range of buffer sizes and types of effects on AADBT estimation. The contribution of this study includes: 1) developed automatic python scripts to extract a wide range of variables at the link, distance, and buffer level (Euclidian and Network), 2) generated and tested large set of variables for both Air and Network buffers, 3) used OpenStreetMap data to standardize the data collection for six different geographies  (Portland, Eugene, Bend, Boulder, Charlotte, and Dallas) of the USA, 4) used emerging data Strava, and StreetLight with traditional data, 5) developed 10 folds-5 repeats cross-validated Poisson regression model to estimate AADBT.  The results indicate that 0.5 miles of air and network combined buffers perform best for the generalized model which is followed by all combined buffer techniques. Network buffer works better than Air buffer when all sizes of buffers are combined, however, both air and network buffer have similar performance for 0.5-mile buffer size. Overall, the generalized model obtained goodness of fit, R2 0.75. The city-specific models indicate that the buffer size and types are influenced by local characteristics of the geography. However, it is found that none of the cities requires more than a 0.5-mile buffer to obtain the best performance of the model. The Portland model works best with a 0.1-mile air buffer, Eugene for 0.1 miles Air+Network, Dallas, charlotte, Bend, and Boulder perform best with a 0.25-mile air buffer only. Overall air buffer works better than network buffers for the aforementioned best buffer sizes. Generalized models have better performance compared to city-specific models due to model learning variation of diverse network characteristics. The outcomes of this research will help planners and modelers to understand the infrastructures and factors that influence the bicycle activities and the sizes and types of buffers that need to be considered to extract the variables to construct a direct demand model.
Keywords: Strava, StreetLight, Bicycle, Buffer, Euclidian, Network, AADBT, Direct Demand Model, 

### Graphical Abstract:
![image](https://user-images.githubusercontent.com/60245323/166619115-ee0a5efd-570a-40ed-9645-35c39a72ddd5.png)

### Requirements
To install requirements: It is best to create an Anaconda environment and execute the following command:/
```bash
pip install -r requirements.txt
```
Library: An initial step creates a separate environment in desktop python to install necessary libraries such as Geopanda, panda, Network, osmnx, shapely, and matplotlib. Geopanda library reads and processes geospatial data and can be downloaded from the link (https://medium.com/analytics-vidhya/fastest-way-to-install-geopandas-in-jupyter-notebook-on-windows-8f734e11fa2b). Networkx and osmnx libraries extract data from OpenStreetMap such as node density, intersection density, number of lanes, and bicycle facilities. Data Inputs: After installing Geopandas Environment in python, the developed python script reads necessary data inputs including counter locations, Strava network, Bikeshare, OSM (shapefiles from BBBike), LEHD job data, NED raster elevation data, and NHGIS data. Geocode Type: A user needs to enter a geocode for data aggregations. In general, any geopanda file uses its original EPSG which could be replaced with a specific local EPSG. For example, all land use shapefiles for Portland, OR are set to EPSG 4326 (default) but should be changed to EPSG 2838 because this code is specifically assigned to Portland in the GIS environment. This local EPSG ensures the extraction of an accurate size of buffers for data aggregation. Note that the developed script already sets this EPSG information for six study locations. A user could use any given buffer size to generate and aggregate the data. Variables defining networks/facility density such as the number of major generators (e.g., schools, colleges, and universities) and network features (e.g., intersection density, number of lanes, and bicycle facilities) were extracted using direct overpass API through osmnx.  The research team also uses BBBike source to supplement land use datasets. Slope data was extracted using a National Elevation Dataset (NED) raster image. This study extracts an elevation at each node of the link and calculates its slope. The research team assumed that a bridge slope is zero. Weather data was collected from weather underground using direct overpass API. The collected data was automatically aggregated as a single data frame and exported in CSV format for further analysis and modeling. The final data frames are passed to the KeplerGL program to create dynamic and interactive visual maps. KeplerGL is a high-performance web-based application for large-scale geospatial data visual exploration. KeplerGL visualizes data in 2D or 3D space as shown in the Figures 3 and 4 respectively. A user selects variables to visualize and saves an interactive map in HTML format to share with others.
The data is available in 'data/'

## Project Structure
This project is split into multiple python and R programming file:
1. 'Network Buffer Creation Code/' folder provides automatic python scripts that will create the network buffer for each count station.
2. 'Data Extraction Code' folder provides the whole complete python scripts for six regions (Portland, Eugene, Bend, Dallas, Chalotte and Boulder)to extract the data. Air(Eucilidian) and Network buffer have separate python scripts. We have total 12 python scripts to extract the data, 2 scripts for each region. 
3. 'Variable Selection Code/' provides the python scripts that were used to reduce the list of variables for each modeling based on correlation and VIF
4. The Poisson Regression 10 folds-5 repeats model was run using R programming. All of the modeling code is available in folder 'Modelling R Code/'
5. Finally, the Author applied Portland best buffer model network-wide using 0.1 mile air buffer. The application data generation code is available in 'Application code/' folder
6. Modeling data is available in 'data/' folder


