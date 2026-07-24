This repository is for the HL-LHC track tracker.

It has the following code

1. CMSSW analyser to extract HL-LHC outer tracker geometry to a csv file
2. Code to make 3D model by convert csv file to a step file or glb file


# CMSSW analyser

Extracts HL-LHC outer tracker D110 geometry to `D110_T35_outer_tracker_sensor_surfaces.csv` file.

```
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsrel CMSSW_20_1_0_pre1 # el9
cd CMSSW_20_1_0_pre1/src
cmsenv

git clone https://github.com/jaebak/L1_track_finder.git .

scram b -j8

cd ExtractGeometry/OuterTrackerGeometry/test
cmsRun dumpD110OuterTrackerStacks_cfg.py
```

# Code to make 3D model
