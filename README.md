# DevTools-Analysis
The analysis portion of the DevTools software suite.

For all FWLite submissions, please use the submitToBSUB.sh command in the following way (which is also shown in the file "ALL_SUBMIT.sh":

```
./submitToBSUB.sh FWLiteKinFit_BTagUP SUBMIT MiniAOD_SIG_h125a6_MuonSelectionOnly_SEP24 8nh "_TEST"
```
In this script, I grab Cross Sections and SummedWeights values from other locations in my lxplus account. It's fairly convoluted and not well setup, but it works and I haven't had much reason to change it. If you have trouble getting it working, then just create your own file or implement your own method withwhich to grab the values and input into the FWLite files.

The parts that I've added are how to group everything for fitting. The various "FWLiteKinFit" files produce the shifts analyzed in the fitting procedure. To run them, execute the file "ALL_SUBMIT.sh" which will submit them to bsub. After they finish, mv the root files to the BSUB/ directory to use the following scripts. 

After they are finished running, then you need to group the files. Execute the "CombineRooDataSetsData.sh" file to group the data correctly. For the signal files, there are two ways to group the results from the FWLiteKinFit jobs. For RooDataSets, execute "rootMacro_CombineToOneRooDataSetFile.C" in a ROOT terminal. For histograms, run "rootMacro_CombineToOneHistFile.C" in a root terminal. Now you have files that can use the RooFit code.

To Make the comparison between our cleaned reconstruction and the standard tau reconstruction, use the file called "FakeRateAnalysis.py". You will need to run this over the two samples of the different reconstruction techniques. To create the plot, use this file "FakeRateEfficiencyPlotter.py". The paths will need to be reconfigured, but the framework should still work.

