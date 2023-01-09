Data processing on 3D skin samples

Matlab Version：2020a

Algorithm flow：

1.batchIntensity3D.m ：Decompress 3D.OCT files

Select 3D data path：Just select the ‘Data’ folder
Select the ‘Data’ folder, and then 3D.data type file will be generated in the Data folder and saved in the ‘Intensity3D’ folder
2.BatchSkinDataProcess3DV2.m ：Process 3D.data data

3.skinDataProcessV4.m ： Algorithmically extract sample surface from sample data

4.saveSkinDataImageV2_5 ：Image the processing results of the sample

5.jiaozhiceng : Segmentation of the cuticle on 2D sample data

6.jiaozhiceng2 : Calculating parameters for artificial skin texture features

7.LastDR : Calculation of artificial skin scattering coefficients

Select 3D.data path：Just select the ‘Intensity 3D’ folder
Select the 'Intensity 3D' folder for processing, then the processing result of the algorithm is in the 'Day13' folder, and the imaging result is in the 'result'.
