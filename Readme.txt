Readme.txt
--
  Program name: SiteSimilarityCharacterizationUnderSite-specificData
  Developers and contact E-mail: Liang Han, hanliang2023@mail.usts.edu.cn, sxrhanliang@126.com;
Wengang Zhang, zhangwg@cqu.edu.cn.
  Software: This program has been used in MATLAB 2021a and MATLAB 2022b. Maybe it can also be run in neighboring 
versions of MATLAB.
  Hardware: This program has been successfully run in my personal computer (PC). My PC configuration is: 12th Gen 
Intel(R) Core(TM) i5-12500H CPU @ 2.50 GHz 16.0 GB RAM.
--
  This program is developed to characterize site similarity under the site-specific data scenario, which includes
two modules: (1) statistical uncertainty characterization of borebole data, (2) Site similarity characterization
under site-specific data. 
  The 1st module is performed by the program named as "StatisticalUncertaintyCharacterizationForBoreholeData", 
which characterizes the statistical uncertainty of the site-specific data in the site similarity characterization. 
If the user do not want to consider the statistical uncertainty contained in the site-specific data, this module can 
be ignored. Nevertheless, the sub-program named as "InitialParametersEstimation_master" in 1st module can be used to
estimate the SOF values (scale of fluctuation) of borehole data by the maximum likelihood method (MLM), which is required
by site similarity characterization.
  The 2nd module is the main program, which is named as "SiteSimilarityCharacterizationUnderSite-specificData_main".
The users the characterize the site similarity under the site-specific data scenario by following the user manual.
--
  In the example, the data for study is mainly from the CLAY/10/7490 database avaliable in the TC304 website 
(http://140.112.12.21/issmge/tc304.htm?=1). Specially, the target data (i.e. geo-material parametric dataset in the Onsøy 
site in Norway) is from the reference of Sharma et al. (2022).
--
  To help the users to operate the program, the output files and test data are provided in the folders named as "Output files"
and "Test data". In these two files, the files are storaged in two folders for statistical uncertainty characterization and 
site similarity characterization, respectively, and these folders have been marked using text information distinctively.
--
Tips:
  When you copy and paste the code files into a folder in the MATLAB path such as the folder named as "bin",  if it reminds you
that "File not found in current folder or MATLAB path", you can type a space in the blank region of the code and save it.
--
Acknowledgements:
  The developers would like to thank the members of the TC304 Committee on Engineering Practice of Risk Assessment & Management
of the International Society of Soil Mechanics and Geotechnical Engineering for developing the 304dB database (http://140.112.12.21/issmge/tc304.htm?=1) used in this study and making it available for scientific inquiry.
--
References:
Sharma, A., Ching, J.Y., Phoon, K.K., 2022. A Hierarchical Bayesian Similarity Measure for Geotechnical Site Retrieval. Journal
  of Engineering Mechanics 148.
