#Where to Start

##Input
The program requires different input files. The program control is defined in the file ZAoptcps.con. The other input files are specified in the control file.

###ZAoptcps.con
format of example control file:

```
d:\SouthAfricaNCEP\ERA40\ERA40 7bina.dat 1958
500.0 10.0 100.0 0.0015
1980 01 01 00
1999 12 31 12
d:\SouthAfricaNCEP\ERA40\ZaERA40Avg.dat d:\SouthAfricaNCEP\ERA40\ZaERA40Std.dat
2
5.0
1 2 3 10 11 12 0 0 0 0 0 0
1.0
456789000000
-10.0 0.0
-50.0 -50.0
CPEA09ER40.def
4
2
1
2 d:\SouthAfricaNCEP\ERA40\CSAG2CP\CSAG2CP\CADATS2006.dat 48 51
```
The lines of the file are described one by one.

*d:\SouthAfricaNCEP\ERA40\ERA40 7bina.dat*<br>
This is the file containing the binary dataset with the 700 hPa data.<br>
*1958*<br>
Starting year for the binary data (used for record number calculation). 1958 for ERA 40 and 1979 for ERA interim<br>
*500.0 10.0 100.0 0.0015*<br>
Weights for the different parts of the objective function.<br>
*1980 01 01 00*<br>
Starting date for optimization<br>
*1999 12 31 00*<br>
Ending date of the optimization.<br>
*d:\SouthAfricaNCEP\ERA40\ZaERA40Avg.dat*<br>
Datafile containing the smoothed mean grids of daily 700hPa values. (annual cycle)<br>
*d:\SouthAfricaNCEP\ERA40\ZaERA40Std.dat*<br>
Datafile containing the smoothed standard deviation of the 700hPa values. (annual cycle)<br>
*2*<br>
Number of seasons to be considered.
The following lines are repeated for each season:<br>
*5.0*<br>
Weight of the objective functions corresponding to season 1. Wet (or strong signal) season should have a higher weight.<br>
*1 2 3 10 11 12 0 0 0 0 0 0*<br>
Months corresponding to season - 12 numbers filled with 0-s.<br>
. . .<br>
After each season is declared:<br>
*-10.0 0.0*<br>
Latitude and longitude of pressure grid upper left corner (negative directions are South and East).<br>
*-50.0 -50.0*<br>
Latitude and longitude of pressure grid lower right corner.<br>
*CPEA09ER40.def*<br>
Initial rule definition file. The number of rules to be obtained by classification is the number of
rules in this file.<br>
*4*<br>
Number of columns in the signal (precipitation) file.<br>
*2*<br>
Number of columns to be used for optimization.<br>
*1*<br>
Index of the first column to be used.<br>
*3*<br>
Index of the second column to be used.<br>
. . .<br>
repeated for each column to be used.<br>*
*d:\SouthAfricaNCEP\ERA40\CSAG2CP\CSAG2CP\CADAT$$$$.dat*<br>
 File containing the signal (precipitation statistics).