Optimize
========

PROGRAM_NAME is a simulation / optimization procedure that estimates evolutionary distance between two sequences, mean of indel length, and rate of insertion and rate of deletion per substitution, given raw alignment sequences in a FASTA file format. 

In order to run PROGRAM_NAME, you need to install R on your computer.  Download the most recent version of R from http://www.r-project.org/ and follow installation instructions. The library "optparse" also needs to be downloaded once R is installed.  To install this package, run the command as follows: 

install.packages("optparse", repose="http://R-Forge.R-Project.org") 

You will also needs to download and install Ngila Release-### or greater from:

  http://scit.us/projects/ngila/ 

You will also need to download and install Dawg Release-### or greater from: 

	http://scit.us/projects/dawg/


After Dawg and Ngila have been installed, the path for these programs can be added to the user path.  If you prefer not to add to the path, please edit the variable dawg_bin in Optimize.R and ngila_bin in ngila_wrapper.cpp to the correct path of the programs.  
PROGRAM_NAME runs from the command line.  To open the command prompt from Windows, you can use the start menu (e.g. Start -> All Programs -> Accessories -> Command Prompt).   Change the directory to the location of the files, e.g.

cd %USERPROFILE%/Path/to/files

Compile the ngila_wrapper.cpp by typing: 

make 

This should produce the object file ngila_wrapper.o and executable NgilaWrapper.  If at any time you would like to remove the object file and the executable, type the following command: 

make clean

It should be noted however, that if you want to use the PROGRAM_NAME properly, ngila_wrapper.cpp should be compiled correctly first. 

Help for the command line arguments of Optimize.R can be accessed by typing: 

./Optimize.R –h 
