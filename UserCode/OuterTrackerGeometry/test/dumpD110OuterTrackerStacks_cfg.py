import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C17I13M9_cff import (
    Phase2C17I13M9,
)

process = cms.Process(
    "DUMPOT",
    Phase2C17I13M9,
)

process.load(
    "Configuration.Geometry."
    "GeometryExtendedRun4D110Reco_cff"
)

# Use the nominal ideal XML geometry rather than applying
# conditions-based alignment corrections.
process.trackerGeometry.applyAlignment = False

process.source = cms.Source(
    "EmptySource"
)

process.maxEvents = cms.untracked.PSet(
    input=cms.untracked.int32(1)
)

process.dumpOuterTrackerStacks = cms.EDAnalyzer(
    "DumpOuterTrackerStacks",
    outputFile=cms.string(
        "D110_T35_outer_tracker_sensor_surfaces.csv"
    ),
)

process.path = cms.Path(
    process.dumpOuterTrackerStacks
)
