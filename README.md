# LuminosityOptimization

Thesis Project 

@Giulia Faletti

## Lists of files present in the directory

  - [LoadData](https://github.com/GiuliaFaletti/LuminosityOptimization/blob/main/LoadData.py): Module that defines different functions able to extract data from FillData.xlsx and TurnAroundData.xlsx. It is important to say that the loaded data have been previously cut considering only the turnaround times that can be used for statistical purposes, and the fills defined "physics fills";
  - [LuminosityOptimization](https://github.com/GiuliaFaletti/LuminosityOptimization/blob/main/LuminosityOptimization.py): Module that defines different functions that evaluate the optimization model parameters, like the optimized fill times;
  - [Models](https://github.com/GiuliaFaletti/LuminosityOptimization/blob/main/Models.py)
  - [CreatingVariableBins](https://github.com/GiuliaFaletti/LuminosityOptimization/blob/main/CreatingVariableBins.py)
  - [OnlineOptimize](https://github.com/GiuliaFaletti/LuminosityOptimization/blob/main/OnlineOptimize.py)
  - [FillsLuminosityEvolution](https://github.com/GiuliaFaletti/LuminosityOptimization/blob/main/FillsLuminosityEvolution.py)
  - [FillsLuminosityEvolution_parallel](https://github.com/GiuliaFaletti/LuminosityOptimization/blob/main/FillsLuminosityEvolution_parallel.py)
  - [TurnArounData](https://github.com/GiuliaFaletti/LuminosityOptimization/blob/main/TurnAroundData.xlsx)
  - [FillData](https://github.com/GiuliaFaletti/LuminosityOptimization/blob/main/FillData.xlsx)
  - [Previous_Year_Data](https://github.com/GiuliaFaletti/LuminosityOptimization/blob/main/Previous_Year_Data.xlsx)
  - [Current_Year_Data](https://github.com/GiuliaFaletti/LuminosityOptimization/blob/main/Current_Year_Data.xlsx)
  - [ATLAS](https://github.com/GiuliaFaletti/LuminosityOptimization/blob/main/ATLAS.zip):
  
      -------------------------- Input Floders ---------------------------------------------------

       - ATLAS_fill_2016: Atlas lumi files for 2016, whose detailed description is available on   
         [https://lpc.web.cern.ch/MassiFileDefinition_v2.htm].
       - ATLAS_fill_2017: Atlas lumi files for 2017, whose detailed description is available on
          [https://lpc.web.cern.ch/MassiFileDefinition_v2.htm];
       - ATLAS_fill_2018: Atlas lumi files for 2018, whose detailed description is available on
                        [https://lpc.web.cern.ch/MassiFileDefinition_v2.htm];
       - ATLAS_summary_2016: Atlas summary files for 2016, whose detailed description is available on
                        [https://lpc.web.cern.ch/MassiFileDefinition_v2.htm];
       - ATLAS_summary_2017: Atlas summary files for 2016, whose detailed description is available on
                        [https://lpc.web.cern.ch/MassiFileDefinition_v2.htm];
       - ATLAS_summary_2018: Atlas summary files for 2016, whose detailed description is available on 
                        [https://lpc.web.cern.ch/MassiFileDefinition_v2.htm].
                       
       -------------------------- Output Floders ---------------------------------------------------

       - FillsLuminosityEvolution2016: Plots of the Luminosity evolution of 2016;
       - FillsLuminosityEvolution2017: Plots of the Luminosity evolution of 2017;
       - FillsLuminosityEvolution2018: Plots of the Luminosity evolution of 2018;
       - OptimalFillsLuminosityEvolution2016: Hypotetical fill plots of the Luminosity evolution of 2016;
       - OptimalFillsLuminosityEvolution2017: Hypotetical fill plots of the Luminosity evolution of 2017;
       - OptimalFillsLuminosityEvolution2018: Hypotetical fill plots of the Luminosity evolution of 2018.
