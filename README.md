# LuminosityOptimization

Thesis Project 

@Giulia Faletti

## Lists of files present in the directory

  - [Optimize](https://github.com/GiuliaFaletti/LuminosityOptimization/blob/main/Optimize.py): Complete user friendly code that halp the user to go through different phases of the thesis work;
  - [OnlineOptimize](https://github.com/GiuliaFaletti/LuminosityOptimization/blob/main/OnlineOptimize.py): Code that analysing the turnaround times is able to return the optimized values of fill time;
  - [LoadData](https://github.com/GiuliaFaletti/LuminosityOptimization/blob/main/LoadData.py): Module that defines different functions able to extract data from FillData.xlsx and TurnAroundData.xlsx. It is important to say that the loaded data have been previously cut considering only the turnaround times that can be used for statistical purposes, and the fills defined "physics fills";
  - [LuminosityOptimization](https://github.com/GiuliaFaletti/LuminosityOptimization/blob/main/LuminosityOptimization.py): Module that defines different functions that evaluate the optimization model parameters, like the optimized fill times;
  - [Models](https://github.com/GiuliaFaletti/LuminosityOptimization/blob/main/Models.py): Module that defines the fit functions for the statistical analysis of the turnaround times;
  - [CreatingVariableBins](https://github.com/GiuliaFaletti/LuminosityOptimization/blob/main/CreatingVariableBins.py): Module that define the algorithm to create variable bins in the histograms;
  - [FillsLuminosityEvolution](https://github.com/GiuliaFaletti/LuminosityOptimization/blob/main/FillsLuminosityEvolution.py): Code that extracts from ATLAS_fill_{year} the Luminosity and Time arrays and plots them on the y and x-axis, saving them in FillsLuminosityEvolution{year};
  - [FillsLuminosityEvolution_parallel](https://github.com/GiuliaFaletti/LuminosityOptimization/blob/main/FillsLuminosityEvolution_parallel.py): Parallelized code that extracts from ATLAS_fill_{year} the Luminosity and Time arrays and plots them on the y and x-axis, saving them in FillsLuminosityEvolution{year};
  - [TurnArounData](https://github.com/GiuliaFaletti/LuminosityOptimization/blob/main/TurnAroundData.xlsx): Excel file with sorted statistical samples of turn around times (sample16, sample17, sample18).
  - [FillData](https://github.com/GiuliaFaletti/LuminosityOptimization/blob/main/FillData.xlsx): Excel file with turnaround times and fill times for each year (t16, tf16, ta17, tf17, ta18, tf18), fills numbers (NrFill_2016, NrFill_2017, Nr_fill2018) and statistical samples of turn around times (sample16, sample17, sample18);
  - [Previous_Year_Data](https://github.com/GiuliaFaletti/LuminosityOptimization/blob/main/Previous_Year_Data.xlsx): Excel file to insert the previous year data and machine parameters in code _OnlineOptimize.py_;
  - [Current_Year_Data](https://github.com/GiuliaFaletti/LuminosityOptimization/blob/main/Current_Year_Data.xlsx): Excel file to insert the current analyzed year data and machine parameters in code _OnlineOptimize.py_;
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
