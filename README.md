# DevTools-Analysis
The analysis portion of the DevTools software suite.


The parts that I've added are how to group everything for fitting. The various "FWLiteKinFit" files produce the shifts analyzed in the fitting procedure. To run them, execute the file "ALL_SUBMIT.sh" which will submit them to bsub. After they finish, mv the root files to the BSUB/ directory to use the following scripts.

After they are finished running, then you need to group the files. Execute the "CombineRooDataSetsData.sh" file to group the data correctly. For the signal files, there are two ways to group the results from the FWLiteKinFit jobs. For RooDataSets, execute "rootMacro_CombineToOneRooDataSetFile.C" in a ROOT terminal. For histograms, run "rootMacro_CombineToOneHistFile.C" in a root terminal. Now you have files that can use the RooFit code.
