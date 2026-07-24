This repository is for the HL-LHC track tracker.
It has the following code

1. CMSSW analyser to extract HL-LHC outer tracker geometry to a csv file.
2. Code to make 3D model of HL-LHC outer tracker by using csv file.
3. Code to make R-Z view of HL-LHC outer tracker by using csv file.

The output of the codes are below
- `D110_T35_outer_tracker_sensor_surfaces.csv`
- `tracker_segmented_glb_files` folder
- `pictures` folder

The csv file and the 3D model can be found below
- `lxplus:~jaebak/work_public/L1_track_finder/D110_T35_outer_tracker_sensor_surfaces.csv`
- `lxplus:~jaebak/work_public/L1_track_finder/tracker_segmented_glb_files`

# CMSSW analyser to extract HL-LHC outer tracker geometry

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

# Code to make 3D model of HL-LHC outer tracker

I run the below code on my Mac where I use `uv` for managing python packages.

More information on `uv` can be found below:

https://docs.astral.sh/uv/getting-started/installation/

```
cd ConvertTo3DModel
uv venv
source .venv/bin/activate
uv pip install cadquery

# Copy D110_T35_outer_tracker_sensor_surfaces.csv to ConvertTo3DModel folder

# It is instructive to see the cell grid of the sensors.
# But to reduce the rendering burden, the cell grid is simplified as below
# - PS module pixel cell (0.100 x  1.467 mm) grid: 960  x 32 -> 1 x 8 (simplified)
# - PS module strip cell (0.100 x 23.554 mm) grid: 960  x 2  -> 1 x 2 (simplified)
# - 2S module strip cell (0.090 x 50.250 mm) grid: 1016 x 2  -> 1 x 2 (simplified)

mkdir tracker_segmented_glb_files
# Note that the first command will take some time due to the cadquery python library downloading files.
./draw_tracker_segmented_cadquery.py D110_T35_outer_tracker_sensor_surfaces.csv tracker_segmented_glb_files/D110_T35_TOB_layer1_surface.glb --subdet TOB --layer 1 --z-horizontal
./draw_tracker_segmented_cadquery.py D110_T35_outer_tracker_sensor_surfaces.csv tracker_segmented_glb_files/D110_T35_TOB_layer2_surface.glb --subdet TOB --layer 2 --z-horizontal
./draw_tracker_segmented_cadquery.py D110_T35_outer_tracker_sensor_surfaces.csv tracker_segmented_glb_files/D110_T35_TOB_layer3_surface.glb --subdet TOB --layer 3 --z-horizontal
./draw_tracker_segmented_cadquery.py D110_T35_outer_tracker_sensor_surfaces.csv tracker_segmented_glb_files/D110_T35_TOB_layer4_surface.glb --subdet TOB --layer 4 --z-horizontal
./draw_tracker_segmented_cadquery.py D110_T35_outer_tracker_sensor_surfaces.csv tracker_segmented_glb_files/D110_T35_TOB_layer5_surface.glb --subdet TOB --layer 5 --z-horizontal
./draw_tracker_segmented_cadquery.py D110_T35_outer_tracker_sensor_surfaces.csv tracker_segmented_glb_files/D110_T35_TOB_layer6_surface.glb --subdet TOB --layer 6 --z-horizontal

./draw_tracker_segmented_cadquery.py D110_T35_outer_tracker_sensor_surfaces.csv tracker_segmented_glb_files/D110_T35_TID_disk1_surface.glb --subdet TID --layer 1 --z-horizontal
./draw_tracker_segmented_cadquery.py D110_T35_outer_tracker_sensor_surfaces.csv tracker_segmented_glb_files/D110_T35_TID_disk2_surface.glb --subdet TID --layer 2 --z-horizontal
./draw_tracker_segmented_cadquery.py D110_T35_outer_tracker_sensor_surfaces.csv tracker_segmented_glb_files/D110_T35_TID_disk3_surface.glb --subdet TID --layer 3 --z-horizontal
./draw_tracker_segmented_cadquery.py D110_T35_outer_tracker_sensor_surfaces.csv tracker_segmented_glb_files/D110_T35_TID_disk4_surface.glb --subdet TID --layer 4 --z-horizontal
./draw_tracker_segmented_cadquery.py D110_T35_outer_tracker_sensor_surfaces.csv tracker_segmented_glb_files/D110_T35_TID_disk5_surface.glb --subdet TID --layer 5 --z-horizontal

# One can open the 3D glb models with 3D viewers such as Godot.

```

# Code to make R-Z view of HL-LHC outer tracker
```
# Using above python setup
source .venv/bin/activate

./plot_cms_tracker_rz.py D110_T35_outer_tracker_sensor_surfaces.csv pictures/tracker_rz.pdf
```
