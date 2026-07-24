#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/ESGetToken.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/GeometrySurface/interface/Bounds.h"
#include "DataFormats/GeometrySurface/interface/RectangularPlaneBounds.h"
#include "DataFormats/GeometrySurface/interface/Surface.h"
#include "DataFormats/GeometrySurface/interface/TrapezoidalPlaneBounds.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"

#include "Geometry/CommonTopologies/interface/GeomDetType.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

using LocalCorners = std::array<LocalPoint, 4>;
using GlobalCorners = std::array<GlobalPoint, 4>;

struct RadialRange {
  double minimum;
  double maximum;
};

struct ShapeDimensions {
  std::string shape;

  // Maximum full width. For a rectangle this is the ordinary width.
  double width;

  // For a rectangle these are both equal to width.
  double bottomWidth;
  double topWidth;

  double length;
  double thickness;
};

/*
 * Local sensor corner ordering:
 *
 *          local +y
 *             ^
 *             |
 *
 *       3-----------2
 *       |           |
 *       |           |
 *       0-----------1  ---> local +x
 *
 * For trapezoids:
 *   corners 0,1 = bottom edge at local -y
 *   corners 2,3 = top edge at local +y
 */

template <typename Point>
double cylindricalR(Point const& point) {
  return std::hypot(
      static_cast<double>(point.x()),
      static_cast<double>(point.y()));
}

template <typename Point>
double radius3D(Point const& point) {
  const double x = point.x();
  const double y = point.y();
  const double z = point.z();

  return std::sqrt(x * x + y * y + z * z);
}

template <typename Point>
double globalPhi(Point const& point) {
  return std::atan2(
      static_cast<double>(point.y()),
      static_cast<double>(point.x()));
}

ShapeDimensions getShapeDimensions(Bounds const& bounds) {
  if (auto const* rectangle =
          dynamic_cast<RectangularPlaneBounds const*>(&bounds)) {
    return ShapeDimensions{
        "rectangle",
        rectangle->width(),
        rectangle->width(),
        rectangle->width(),
        rectangle->length(),
        rectangle->thickness()};
  }

  if (auto const* trapezoid =
          dynamic_cast<TrapezoidalPlaneBounds const*>(&bounds)) {
    /*
     * TrapezoidalPlaneBounds::parameters() returns:
     *
     *   [0] half bottom-edge width
     *   [1] half top-edge width
     *   [2] half thickness
     *   [3] half apothem / half length
     */
    const auto parameters = trapezoid->parameters();

    const double bottomWidth = 2.0 * parameters[0];
    const double topWidth = 2.0 * parameters[1];

    return ShapeDimensions{
        "trapezoid",
        trapezoid->width(),
        bottomWidth,
        topWidth,
        trapezoid->length(),
        trapezoid->thickness()};
  }

  throw std::runtime_error(
      "Unsupported sensor bounds type. Expected "
      "RectangularPlaneBounds or TrapezoidalPlaneBounds.");
}

LocalCorners makeLocalCorners(
    Bounds const& bounds,
    float localZ) {
  if (auto const* rectangle =
          dynamic_cast<RectangularPlaneBounds const*>(&bounds)) {
    const float halfWidth =
        0.5f * rectangle->width();

    const float halfLength =
        0.5f * rectangle->length();

    return LocalCorners{
        LocalPoint(
            -halfWidth,
            -halfLength,
            localZ),
        LocalPoint(
            +halfWidth,
            -halfLength,
            localZ),
        LocalPoint(
            +halfWidth,
            +halfLength,
            localZ),
        LocalPoint(
            -halfWidth,
            +halfLength,
            localZ)};
  }

  if (auto const* trapezoid =
          dynamic_cast<TrapezoidalPlaneBounds const*>(&bounds)) {
    const auto parameters = trapezoid->parameters();

    const float halfBottomWidth = parameters[0];
    const float halfTopWidth = parameters[1];
    const float halfLength = parameters[3];

    return LocalCorners{
        LocalPoint(
            -halfBottomWidth,
            -halfLength,
            localZ),
        LocalPoint(
            +halfBottomWidth,
            -halfLength,
            localZ),
        LocalPoint(
            +halfTopWidth,
            +halfLength,
            localZ),
        LocalPoint(
            -halfTopWidth,
            +halfLength,
            localZ)};
  }

  throw std::runtime_error(
      "Unsupported sensor bounds type while creating corners.");
}

GlobalCorners makeGlobalCorners(
    Surface const& surface,
    LocalCorners const& localCorners) {
  return GlobalCorners{
      surface.toGlobal(localCorners[0]),
      surface.toGlobal(localCorners[1]),
      surface.toGlobal(localCorners[2]),
      surface.toGlobal(localCorners[3])};
}

/*
 * Distance in the global x-y plane between the beam axis and
 * the line segment connecting points a and b.
 */
double distanceToSegmentXY(
    GlobalPoint const& a,
    GlobalPoint const& b) {
  const double ax = a.x();
  const double ay = a.y();

  const double dx = b.x() - a.x();
  const double dy = b.y() - a.y();

  const double lengthSquared =
      dx * dx + dy * dy;

  if (lengthSquared <=
      std::numeric_limits<double>::epsilon()) {
    return std::hypot(ax, ay);
  }

  double fraction =
      -(ax * dx + ay * dy) / lengthSquared;

  fraction = std::clamp(
      fraction,
      0.0,
      1.0);

  const double closestX =
      ax + fraction * dx;

  const double closestY =
      ay + fraction * dy;

  return std::hypot(
      closestX,
      closestY);
}

/*
 * Test whether the beam axis lies within the x-y projection
 * of the quadrilateral.
 */
bool beamAxisInsideProjectedFace(
    GlobalCorners const& corners) {
  bool inside = false;

  for (std::size_t i = 0,
                   j = corners.size() - 1;
       i < corners.size();
       j = i++) {
    const double xi = corners[i].x();
    const double yi = corners[i].y();

    const double xj = corners[j].x();
    const double yj = corners[j].y();

    const bool crossesYZero =
        ((yi > 0.0) != (yj > 0.0));

    if (!crossesYZero) {
      continue;
    }

    const double xAtYZero =
        xi + (xj - xi) * (-yi) / (yj - yi);

    if (xAtYZero > 0.0) {
      inside = !inside;
    }
  }

  return inside;
}

/*
 * Calculate the exact cylindrical-R interval covered by the
 * quadrilateral face.
 *
 * Maximum R occurs at a corner.
 * Minimum R can occur:
 *   - at a corner;
 *   - on an edge;
 *   - at R=0 if the projected face contains the beam axis.
 */
RadialRange calculateFaceRadialRange(
    GlobalCorners const& corners) {
  double maximumR = 0.0;

  for (auto const& corner : corners) {
    maximumR =
        std::max(
            maximumR,
            cylindricalR(corner));
  }

  if (beamAxisInsideProjectedFace(corners)) {
    return RadialRange{
        0.0,
        maximumR};
  }

  double minimumR =
      std::numeric_limits<double>::max();

  for (std::size_t index = 0;
       index < corners.size();
       ++index) {
    const std::size_t next =
        (index + 1) % corners.size();

    minimumR =
        std::min(
            minimumR,
            distanceToSegmentXY(
                corners[index],
                corners[next]));
  }

  return RadialRange{
      minimumR,
      maximumR};
}

std::string detectorTypeName(
    TrackerGeometry::ModuleType detectorType) {
  using ModuleType = TrackerGeometry::ModuleType;

  switch (detectorType) {
    case ModuleType::Ph2PSP:
      return "Ph2PSP";

    case ModuleType::Ph2PSS:
      return "Ph2PSS";

    case ModuleType::Ph2SS:
      return "Ph2SS";

    case ModuleType::UNKNOWN:
      return "UNKNOWN";

    default:
      return "OTHER";
  }
}

std::string sensorTechnology(
    TrackerGeometry::ModuleType detectorType) {
  using ModuleType = TrackerGeometry::ModuleType;

  switch (detectorType) {
    case ModuleType::Ph2PSP:
      return "macro_pixel";

    case ModuleType::Ph2PSS:
      return "strip";

    case ModuleType::Ph2SS:
      return "strip";

    default:
      return "unknown";
  }
}

std::string stackModuleType(
    TrackerGeometry::ModuleType lowerType,
    TrackerGeometry::ModuleType upperType) {
  using ModuleType = TrackerGeometry::ModuleType;

  const bool lowerIsPS =
      lowerType == ModuleType::Ph2PSP ||
      lowerType == ModuleType::Ph2PSS;

  const bool upperIsPS =
      upperType == ModuleType::Ph2PSP ||
      upperType == ModuleType::Ph2PSS;

  const bool hasPSPixel =
      lowerType == ModuleType::Ph2PSP ||
      upperType == ModuleType::Ph2PSP;

  const bool hasPSStrip =
      lowerType == ModuleType::Ph2PSS ||
      upperType == ModuleType::Ph2PSS;

  if (lowerIsPS &&
      upperIsPS &&
      hasPSPixel &&
      hasPSStrip) {
    return "PS";
  }

  if (lowerType == ModuleType::Ph2SS &&
      upperType == ModuleType::Ph2SS) {
    return "2S";
  }

  return "unknown";
}

void writeFace(
    std::ofstream& output,
    std::uint32_t stackRawId,
    char const* subdetector,
    unsigned int layerOrDisk,
    unsigned int side,
    char const* sensorRole,
    DetId const& sensorId,
    std::string const& geomTypeName,
    TrackerGeometry::ModuleType detectorType,
    std::string const& moduleType,
    Surface const& surface,
    char const* faceName,
    float localZ) {
  const Bounds& bounds =
      surface.bounds();

  const ShapeDimensions dimensions =
      getShapeDimensions(bounds);

  const LocalCorners localCorners =
      makeLocalCorners(
          bounds,
          localZ);

  const GlobalCorners globalCorners =
      makeGlobalCorners(
          surface,
          localCorners);

  const RadialRange radialRange =
      calculateFaceRadialRange(
          globalCorners);

  auto const& center =
      surface.position();

  const std::string detectorTypeString =
      detectorTypeName(detectorType);

  const std::string technology =
      sensorTechnology(detectorType);

  for (std::size_t cornerIndex = 0;
       cornerIndex < localCorners.size();
       ++cornerIndex) {
    const LocalPoint& local =
        localCorners[cornerIndex];

    const GlobalPoint& global =
        globalCorners[cornerIndex];

    output
        << stackRawId << ','
        << subdetector << ','
        << layerOrDisk << ','
        << side << ','
        << sensorRole << ','
        << sensorId.rawId() << ','
        << geomTypeName << ','
        << detectorTypeString << ','
        << technology << ','
        << moduleType << ','
        << dimensions.shape << ','
        << faceName << ','
        << cornerIndex << ','

        << dimensions.width << ','
        << dimensions.bottomWidth << ','
        << dimensions.topWidth << ','
        << dimensions.length << ','
        << dimensions.thickness << ','

        << center.x() << ','
        << center.y() << ','
        << center.z() << ','
        << cylindricalR(center) << ','
        << radius3D(center) << ','
        << globalPhi(center) << ','

        << local.x() << ','
        << local.y() << ','
        << local.z() << ','

        << global.x() << ','
        << global.y() << ','
        << global.z() << ','
        << cylindricalR(global) << ','
        << radius3D(global) << ','
        << globalPhi(global) << ','

        << radialRange.minimum << ','
        << radialRange.maximum
        << '\n';
  }
}

void writeSensor(
    std::ofstream& output,
    std::uint32_t stackRawId,
    char const* subdetector,
    unsigned int layerOrDisk,
    unsigned int side,
    char const* sensorRole,
    DetId const& sensorId,
    std::string const& geomTypeName,
    TrackerGeometry::ModuleType detectorType,
    std::string const& moduleType,
    Surface const& surface) {
  const float halfThickness =
      0.5f * surface.bounds().thickness();

  /*
   * The actual CMSSW detector surface is the local-z = 0
   * central plane.
   */
  writeFace(
      output,
      stackRawId,
      subdetector,
      layerOrDisk,
      side,
      sensorRole,
      sensorId,
      geomTypeName,
      detectorType,
      moduleType,
      surface,
      "mid_plane",
      0.0f);

  /*
   * These two planes delimit the geometrical thickness
   * represented by Surface::bounds().
   */
  writeFace(
      output,
      stackRawId,
      subdetector,
      layerOrDisk,
      side,
      sensorRole,
      sensorId,
      geomTypeName,
      detectorType,
      moduleType,
      surface,
      "local_z_minus",
      -halfThickness);

  writeFace(
      output,
      stackRawId,
      subdetector,
      layerOrDisk,
      side,
      sensorRole,
      sensorId,
      geomTypeName,
      detectorType,
      moduleType,
      surface,
      "local_z_plus",
      +halfThickness);
}

}  // namespace

class DumpOuterTrackerStacks
    : public edm::one::EDAnalyzer<> {
public:
  explicit DumpOuterTrackerStacks(
      edm::ParameterSet const& config)
      : outputFileName_(
            config.getParameter<std::string>(
                "outputFile")),
        geometryToken_(
            esConsumes<
                TrackerGeometry,
                TrackerDigiGeometryRecord>()),
        topologyToken_(
            esConsumes<
                TrackerTopology,
                TrackerTopologyRcd>()) {}

  void analyze(
      edm::Event const&,
      edm::EventSetup const& eventSetup) override {
    /*
     * Normally the configuration runs one EmptySource event.
     * This guard prevents accidental rewriting if more events
     * are configured.
     */
    if (hasWritten_) {
      return;
    }

    hasWritten_ = true;

    auto const& geometry =
        eventSetup.getData(geometryToken_);

    auto const& topology =
        eventSetup.getData(topologyToken_);

    std::ofstream output(outputFileName_);

    if (!output.is_open()) {
      throw std::runtime_error(
          "Could not open output file: " +
          outputFileName_);
    }

    output << std::setprecision(15);

    output
        << "stack_rawid,"
        << "subdet,"
        << "layer_or_disk,"
        << "side,"
        << "sensor_role,"
        << "sensor_rawid,"
        << "geom_type_name,"
        << "detector_type,"
        << "sensor_technology,"
        << "module_type,"
        << "shape,"
        << "face,"
        << "corner,"

        << "width_cm,"
        << "bottom_width_cm,"
        << "top_width_cm,"
        << "length_cm,"
        << "thickness_cm,"

        << "center_global_x_cm,"
        << "center_global_y_cm,"
        << "center_global_z_cm,"
        << "center_global_r_cm,"
        << "center_global_r3d_cm,"
        << "center_global_phi_rad,"

        << "local_x_cm,"
        << "local_y_cm,"
        << "local_z_cm,"

        << "global_x_cm,"
        << "global_y_cm,"
        << "global_z_cm,"
        << "global_r_cm,"
        << "global_r3d_cm,"
        << "global_phi_rad,"

        << "face_min_global_r_cm,"
        << "face_max_global_r_cm"
        << '\n';

    unsigned int numberOfStacks = 0;
    unsigned int numberOfPSSensors = 0;
    unsigned int numberOf2SSensors = 0;
    unsigned int numberOfUnknownSensors = 0;

    for (auto const* detectorUnit :
         geometry.detUnits()) {
      const DetId sensorId =
          detectorUnit->geographicalId();

      /*
       * Phase-2 Outer Tracker:
       *
       *   TOB = barrel
       *   TID = endcap
       */
      const bool isBarrel =
          sensorId.subdetId() ==
          StripSubdetector::TOB;

      const bool isEndcap =
          sensorId.subdetId() ==
          StripSubdetector::TID;

      if (!isBarrel && !isEndcap) {
        continue;
      }

      /*
       * Begin with the lower sensor so each physical stack
       * is processed exactly once.
       */
      if (!topology.isLower(sensorId)) {
        continue;
      }

      const DetId lowerId =
          sensorId;

      const DetId upperId =
          topology.partnerDetId(lowerId);

      if (upperId.rawId() == 0) {
        edm::LogWarning(
            "DumpOuterTrackerStacks")
            << "No partner DetId found for lower sensor "
            << lowerId.rawId();

        continue;
      }

      auto const* lowerDetector =
          geometry.idToDetUnit(lowerId);

      auto const* upperDetector =
          geometry.idToDetUnit(upperId);

      if (lowerDetector == nullptr ||
          upperDetector == nullptr) {
        throw std::runtime_error(
            "Could not resolve both sensors for stack " +
            std::to_string(
                topology.stack(lowerId)));
      }

      const TrackerGeometry::ModuleType lowerType =
          geometry.getDetectorType(lowerId);

      const TrackerGeometry::ModuleType upperType =
          geometry.getDetectorType(upperId);

      const std::string moduleType =
          stackModuleType(
              lowerType,
              upperType);

      const std::uint32_t stackRawId =
          topology.stack(lowerId);

      const char* subdetector =
          isBarrel ? "TOB" : "TID";

      const unsigned int layerOrDisk =
          topology.layer(lowerId);

      const unsigned int side =
          topology.side(lowerId);

      const std::string lowerGeomTypeName =
          lowerDetector->type().name();

      const std::string upperGeomTypeName =
          upperDetector->type().name();

      writeSensor(
          output,
          stackRawId,
          subdetector,
          layerOrDisk,
          side,
          "lower",
          lowerId,
          lowerGeomTypeName,
          lowerType,
          moduleType,
          lowerDetector->surface());

      writeSensor(
          output,
          stackRawId,
          subdetector,
          layerOrDisk,
          side,
          "upper",
          upperId,
          upperGeomTypeName,
          upperType,
          moduleType,
          upperDetector->surface());

      for (TrackerGeometry::ModuleType detectorType :
           {lowerType, upperType}) {
        if (detectorType ==
                TrackerGeometry::ModuleType::Ph2PSP ||
            detectorType ==
                TrackerGeometry::ModuleType::Ph2PSS) {
          ++numberOfPSSensors;
        } else if (
            detectorType ==
            TrackerGeometry::ModuleType::Ph2SS) {
          ++numberOf2SSensors;
        } else {
          ++numberOfUnknownSensors;
        }
      }

      ++numberOfStacks;
    }

    output.close();

    edm::LogPrint(
        "DumpOuterTrackerStacks")
        << "Analyzer version: dimensions and sensor types\n"
        << "Wrote " << numberOfStacks
        << " Outer Tracker stacks to "
        << outputFileName_ << "\n"
        << "PS sensors: " << numberOfPSSensors << "\n"
        << "2S sensors: " << numberOf2SSensors << "\n"
        << "Unknown sensors: "
        << numberOfUnknownSensors;
  }

private:
  const std::string outputFileName_;

  const edm::ESGetToken<
      TrackerGeometry,
      TrackerDigiGeometryRecord>
      geometryToken_;

  const edm::ESGetToken<
      TrackerTopology,
      TrackerTopologyRcd>
      topologyToken_;

  bool hasWritten_ = false;
};

DEFINE_FWK_MODULE(DumpOuterTrackerStacks);
