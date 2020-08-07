# Io_Tidal
The software pertaining to the publication Spencer et al (2020), "Tidal controls on the lithospheric thickness and topography of Io from magmatic segregation and volcanism modelling".

The dynamical model in this work is a generalisation of the model presented in Spencer, Katz, and Hewitt (2020) "Magmatic intrusions control Io's crustal thickness". That used the io_compaction code, which can be found at https://github.com/Spencer-space/io_compaction

This repository contains:
IoTidal.m
One-dimensional dynamics model, generalised form of io_compaction

TidalHeatingCalc.m
Tidal heating model, used to calculate 3D heating or heating at a given latitude and longitude

IoTidal_script.m
Wrapper script used to run the coupled dynamic and tidal codes. Also conducts the isostasy calculations.

CSV files for figure 1
These files are tables of the forms 1000*4, where each column in a specified latitude and longitude, and the row is the vector of points from the CMB to the surface.

CSV files for figures 2 and 3
These files are matrices of plotted data. Rows are longitudes 0 to 360, columns are latitude from 90 to -90