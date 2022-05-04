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
To install requirements: It is best to create an Anaconda environment and execute the following command:
pip install -r requirements.txt

The data is available in 'data/'
