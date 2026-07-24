#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Tuple

import cadquery as cq
import numpy as np
from cadquery.occ_impl.assembly import toCAF
from OCP.Message import Message_ProgressRange
from OCP.RWGltf import RWGltf_CafWriter
from OCP.TCollection import TCollection_AsciiString
from OCP.TColStd import TColStd_IndexedDataMapOfStringString


# CMSSW geometry is in cm; keep all CAD construction in mm.
CM_TO_MM = 10.0

# glTF/GLB defines linear distances in metres.  The custom exporter below
# tells Open CASCADE that one internal model unit is one millimetre.
MM_TO_M = 0.001
GLTF_OUTPUT_UNIT_M = 1.0
GLTF_LINEAR_TOLERANCE_MM = 1.0e-3
GLTF_ANGULAR_TOLERANCE_RAD = 0.1

Point3D = Tuple[float, float, float]
FaceCorners = Dict[int, Point3D]


SEGMENTATION_BY_TYPE = {
    "Ph2PSP": 8,
    "Ph2PSS": 2,
    "Ph2SS": 2,
}


@dataclass(frozen=True)
class SensorKey:
    stack_rawid: int
    sensor_rawid: int
    subdet: str
    layer_or_disk: int
    sensor_role: str
    detector_type: str


@dataclass(frozen=True)
class ProjectionInfo:
    maximum_mm: float
    rms_mm: float


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Draw the Phase-2 Outer Tracker with visible column "
            "segmentation: Ph2PSP x32, Ph2PSS x2, and Ph2SS x2."
        )
    )

    parser.add_argument(
        "input_csv",
        type=Path,
        help="Sensor-surface CSV produced by the CMSSW analyzer",
    )

    parser.add_argument(
        "output",
        type=Path,
        help="Output file: .step, .stp, .gltf, or .glb",
    )

    parser.add_argument(
        "--mode",
        choices=("surface", "solid"),
        default="surface",
        help=(
            "surface: draw segmented sensor mid-planes; "
            "solid: draw segmented sensor volumes"
        ),
    )

    parser.add_argument(
        "--subdet",
        choices=("TOB", "TID"),
        help="Restrict output to TOB or TID",
    )

    parser.add_argument(
        "--layer",
        type=int,
        help="Restrict output to one layer or disk",
    )

    parser.add_argument(
        "--stack-rawid",
        type=int,
        help="Restrict output to one stack",
    )

    parser.add_argument(
        "--sensor-rawid",
        type=int,
        help="Restrict output to one sensor",
    )

    parser.add_argument(
        "--sensor-role",
        choices=("lower", "upper"),
        help="Restrict output to lower or upper sensors",
    )

    parser.add_argument(
        "--detector-type",
        action="append",
        choices=("Ph2PSP", "Ph2PSS", "Ph2SS"),
        help=(
            "Restrict output to a detector type. "
            "Repeat this option to select multiple types."
        ),
    )

    parser.add_argument(
        "--every",
        type=int,
        default=1,
        help="Keep every Nth sensor for a lighter preview",
    )

    parser.add_argument(
        "--max-sensors",
        type=int,
        help="Maximum number of sensors to build",
    )

    parser.add_argument(
        "--gap-mm",
        type=float,
        default=0.0,
        help=(
            "Optional physical gap between neighboring bands. "
            "Default: 0, so the original sensor footprint is preserved."
        ),
    )

    parser.add_argument(
        "--report-projection-above-mm",
        type=float,
        default=None,
        help=(
            "Report sensors whose maximum plane-projection "
            "correction exceeds this value"
        ),
    )

    parser.add_argument(
        "--reject-projection-above-mm",
        type=float,
        default=None,
        help=(
            "Reject sensors whose maximum projection correction "
            "exceeds this value"
        ),
    )

    parser.add_argument(
        "--z-horizontal",
        action="store_true",
        help=(
            "Rotate the model by +90 degrees around global Y so "
            "the original global Z axis becomes the new +X axis"
        ),
    )

    return parser.parse_args()


def row_passes_filters(
    row: dict[str, str],
    args: argparse.Namespace,
) -> bool:
    if args.subdet is not None and row["subdet"] != args.subdet:
        return False

    if (
        args.layer is not None
        and int(row["layer_or_disk"]) != args.layer
    ):
        return False

    if (
        args.stack_rawid is not None
        and int(row["stack_rawid"]) != args.stack_rawid
    ):
        return False

    if (
        args.sensor_rawid is not None
        and int(row["sensor_rawid"]) != args.sensor_rawid
    ):
        return False

    if (
        args.sensor_role is not None
        and row["sensor_role"] != args.sensor_role
    ):
        return False

    if (
        args.detector_type is not None
        and row["detector_type"] not in args.detector_type
    ):
        return False

    if args.mode == "surface":
        return row["face"] == "mid_plane"

    return row["face"] in {
        "local_z_minus",
        "local_z_plus",
    }


def rotate_z_to_horizontal(point: Point3D) -> Point3D:
    """
    Rotate +90 degrees around the global Y axis.

    Coordinate mapping:
        original X -> new -Z
        original Y -> new  Y
        original Z -> new +X

    Thus the original CMS Z/beam axis becomes horizontal along +X.
    """
    x, y, z = point

    return (
        z,
        y,
        -x,
    )


def read_sensor_faces(
    csv_path: Path,
    args: argparse.Namespace,
) -> dict[SensorKey, dict[str, FaceCorners]]:
    sensors: dict[
        SensorKey,
        dict[str, FaceCorners],
    ] = defaultdict(lambda: defaultdict(dict))

    with csv_path.open(newline="") as input_file:
        reader = csv.DictReader(input_file)

        required_columns = {
            "stack_rawid",
            "sensor_rawid",
            "subdet",
            "layer_or_disk",
            "sensor_role",
            "detector_type",
            "face",
            "corner",
            "global_x_cm",
            "global_y_cm",
            "global_z_cm",
        }

        available_columns = set(reader.fieldnames or [])
        missing_columns = required_columns - available_columns

        if missing_columns:
            raise ValueError(
                "Input CSV is missing required columns: "
                + ", ".join(sorted(missing_columns))
            )

        for row in reader:
            if not row_passes_filters(row, args):
                continue

            detector_type = row["detector_type"].strip()

            if detector_type not in SEGMENTATION_BY_TYPE:
                continue

            key = SensorKey(
                stack_rawid=int(row["stack_rawid"]),
                sensor_rawid=int(row["sensor_rawid"]),
                subdet=row["subdet"].strip(),
                layer_or_disk=int(row["layer_or_disk"]),
                sensor_role=row["sensor_role"].strip(),
                detector_type=detector_type,
            )

            corner_number = int(row["corner"])

            original_point = (
                CM_TO_MM * float(row["global_x_cm"]),
                CM_TO_MM * float(row["global_y_cm"]),
                CM_TO_MM * float(row["global_z_cm"]),
            )

            point = (
                rotate_z_to_horizontal(original_point)
                if args.z_horizontal
                else original_point
            )

            sensors[key][row["face"]][corner_number] = point

    return {
        key: dict(faces)
        for key, faces in sensors.items()
    }


def ordered_corners(
    corners: FaceCorners,
    description: str,
) -> list[Point3D]:
    expected = {0, 1, 2, 3}
    available = set(corners)

    if available != expected:
        raise ValueError(
            f"{description} contains corners {sorted(available)}, "
            "but corners 0, 1, 2, and 3 are required"
        )

    return [corners[index] for index in range(4)]


def as_numpy(points: list[Point3D]) -> np.ndarray:
    array = np.asarray(points, dtype=float)

    if array.shape != (4, 3):
        raise ValueError(
            f"Expected four 3D points, received {array.shape}"
        )

    if not np.all(np.isfinite(array)):
        raise ValueError("Coordinates contain NaN or infinity")

    return array


def as_point_list(array: np.ndarray) -> list[Point3D]:
    return [
        (
            float(point[0]),
            float(point[1]),
            float(point[2]),
        )
        for point in array
    ]


def normalize(vector: np.ndarray) -> np.ndarray:
    length = float(np.linalg.norm(vector))

    if length <= 1.0e-14:
        raise ValueError("Cannot normalize a zero-length vector")

    return vector / length


def fit_best_plane(
    points: list[Point3D],
) -> tuple[np.ndarray, np.ndarray]:
    array = as_numpy(points)
    centroid = array.mean(axis=0)
    centered = array - centroid

    _, singular_values, right_vectors = np.linalg.svd(
        centered,
        full_matrices=False,
    )

    if singular_values[1] <= 1.0e-12:
        raise ValueError(
            "Sensor corners are degenerate or nearly collinear"
        )

    normal = normalize(right_vectors[-1])

    return centroid, normal


def project_points_to_plane(
    points: list[Point3D],
    plane_point: np.ndarray,
    plane_normal: np.ndarray,
) -> tuple[list[Point3D], ProjectionInfo]:
    array = as_numpy(points)
    normal = normalize(plane_normal)

    signed_distances = (
        array - plane_point
    ) @ normal

    projected = (
        array
        - signed_distances[:, np.newaxis]
        * normal[np.newaxis, :]
    )

    absolute_distances = np.abs(signed_distances)

    info = ProjectionInfo(
        maximum_mm=float(absolute_distances.max()),
        rms_mm=float(
            math.sqrt(
                float(np.mean(signed_distances**2))
            )
        ),
    )

    return as_point_list(projected), info


def project_quad_to_best_plane(
    points: list[Point3D],
) -> tuple[list[Point3D], ProjectionInfo]:
    centroid, normal = fit_best_plane(points)

    return project_points_to_plane(
        points,
        centroid,
        normal,
    )


def project_parallel_solid_faces(
    minus_points: list[Point3D],
    plus_points: list[Point3D],
) -> tuple[
    list[Point3D],
    list[Point3D],
    ProjectionInfo,
]:
    minus_array = as_numpy(minus_points)
    plus_array = as_numpy(plus_points)

    minus_centroid = minus_array.mean(axis=0)
    plus_centroid = plus_array.mean(axis=0)

    centroid_separation = plus_centroid - minus_centroid

    if float(np.linalg.norm(centroid_separation)) <= 1.0e-12:
        combined = np.vstack(
            [minus_array, plus_array]
        )

        combined_centroid = combined.mean(axis=0)
        centered = combined - combined_centroid

        _, _, right_vectors = np.linalg.svd(
            centered,
            full_matrices=False,
        )

        common_normal = normalize(right_vectors[-1])
    else:
        common_normal = normalize(centroid_separation)

    projected_minus, minus_info = project_points_to_plane(
        minus_points,
        minus_centroid,
        common_normal,
    )

    projected_plus, plus_info = project_points_to_plane(
        plus_points,
        plus_centroid,
        common_normal,
    )

    combined_rms = math.sqrt(
        0.5
        * (
            minus_info.rms_mm**2
            + plus_info.rms_mm**2
        )
    )

    return (
        projected_minus,
        projected_plus,
        ProjectionInfo(
            maximum_mm=max(
                minus_info.maximum_mm,
                plus_info.maximum_mm,
            ),
            rms_mm=combined_rms,
        ),
    )


def distance(
    first: Point3D,
    second: Point3D,
) -> float:
    return math.sqrt(
        (second[0] - first[0]) ** 2
        + (second[1] - first[1]) ** 2
        + (second[2] - first[2]) ** 2
    )


def interpolate(
    first: Point3D,
    second: Point3D,
    fraction: float,
) -> Point3D:
    return (
        first[0] + fraction * (second[0] - first[0]),
        first[1] + fraction * (second[1] - first[1]),
        first[2] + fraction * (second[2] - first[2]),
    )


def subdivide_quad(
    points: list[Point3D],
    divisions: int,
    gap_mm: float,
) -> list[list[Point3D]]:
    """
    Subdivide along the measurement-column/local-y direction.

    Original corner order:

        3 -------- 2
        |          |
        |          |
        0 -------- 1

    The subdivision interpolates along:
        left edge:  corner 0 -> corner 3
        right edge: corner 1 -> corner 2
    """
    if len(points) != 4:
        raise ValueError("Exactly four corners are required")

    if divisions < 1:
        raise ValueError("divisions must be at least 1")

    left_length = distance(points[0], points[3])
    right_length = distance(points[1], points[2])

    average_length = 0.5 * (left_length + right_length)

    if average_length <= 1.0e-12:
        raise ValueError("Sensor has zero length")

    gap_fraction = gap_mm / average_length

    if gap_fraction >= 1.0 / divisions:
        raise ValueError(
            "Requested gap is too large for the subdivision"
        )

    bands: list[list[Point3D]] = []

    for index in range(divisions):
        start_fraction = index / divisions
        end_fraction = (index + 1) / divisions

        # Put half the requested gap on either side of an
        # interior division boundary.
        if index > 0:
            start_fraction += 0.5 * gap_fraction

        if index < divisions - 1:
            end_fraction -= 0.5 * gap_fraction

        lower_left = interpolate(
            points[0],
            points[3],
            start_fraction,
        )

        lower_right = interpolate(
            points[1],
            points[2],
            start_fraction,
        )

        upper_right = interpolate(
            points[1],
            points[2],
            end_fraction,
        )

        upper_left = interpolate(
            points[0],
            points[3],
            end_fraction,
        )

        bands.append(
            [
                lower_left,
                lower_right,
                upper_right,
                upper_left,
            ]
        )

    return bands


def make_wire(points: list[Point3D]) -> cq.Wire:
    return cq.Wire.makePolygon(
        points,
        close=True,
    )


def make_face(points: list[Point3D]) -> cq.Face:
    wire = make_wire(points)
    return cq.Face.makeFromWires(wire)


def make_surface_segments(
    faces: dict[str, FaceCorners],
    divisions: int,
    gap_mm: float,
) -> tuple[list[cq.Face], ProjectionInfo]:
    if "mid_plane" not in faces:
        raise ValueError(
            "Sensor does not contain a mid_plane face"
        )

    original = ordered_corners(
        faces["mid_plane"],
        "mid_plane",
    )

    projected, projection_info = (
        project_quad_to_best_plane(original)
    )

    bands = subdivide_quad(
        projected,
        divisions,
        gap_mm,
    )

    shapes = [
        make_face(band)
        for band in bands
    ]

    return shapes, projection_info


def make_solid_segments(
    faces: dict[str, FaceCorners],
    divisions: int,
    gap_mm: float,
) -> tuple[list[cq.Solid], ProjectionInfo]:
    required = {
        "local_z_minus",
        "local_z_plus",
    }

    missing = required - set(faces)

    if missing:
        raise ValueError(
            "Sensor is missing required faces: "
            + ", ".join(sorted(missing))
        )

    original_minus = ordered_corners(
        faces["local_z_minus"],
        "local_z_minus",
    )

    original_plus = ordered_corners(
        faces["local_z_plus"],
        "local_z_plus",
    )

    (
        projected_minus,
        projected_plus,
        projection_info,
    ) = project_parallel_solid_faces(
        original_minus,
        original_plus,
    )

    minus_bands = subdivide_quad(
        projected_minus,
        divisions,
        gap_mm,
    )

    plus_bands = subdivide_quad(
        projected_plus,
        divisions,
        gap_mm,
    )

    solids: list[cq.Solid] = []

    for minus_band, plus_band in zip(
        minus_bands,
        plus_bands,
        strict=True,
    ):
        minus_wire = make_wire(minus_band)
        plus_wire = make_wire(plus_band)

        solid = cq.Solid.makeLoft(
            [minus_wire, plus_wire],
            ruled=True,
        )

        solids.append(solid)

    return solids, projection_info


def segment_color(
    key: SensorKey,
    segment_index: int,
) -> cq.Color:
    alternate = segment_index % 2

    if key.detector_type == "Ph2PSP":
        # Macro-pixel columns: orange/red.
        if alternate == 0:
            return cq.Color(0.95, 0.45, 0.15, 1.0)

        return cq.Color(0.75, 0.20, 0.10, 1.0)

    if key.detector_type == "Ph2PSS":
        # PS strip segments: blue.
        if alternate == 0:
            return cq.Color(0.20, 0.55, 0.95, 1.0)

        return cq.Color(0.10, 0.30, 0.75, 1.0)

    # 2S strip segments: green.
    if alternate == 0:
        return cq.Color(0.25, 0.80, 0.35, 1.0)

    return cq.Color(0.10, 0.55, 0.20, 1.0)


def make_sensor_name(key: SensorKey) -> str:
    return (
        f"{key.subdet}_"
        f"L{key.layer_or_disk}_"
        f"{key.detector_type}_"
        f"stack{key.stack_rawid}_"
        f"{key.sensor_role}_"
        f"sensor{key.sensor_rawid}"
    )


def validate_output_extension(
    output_path: Path,
) -> None:
    supported = {
        ".step",
        ".stp",
        ".gltf",
        ".glb",
    }

    extension = output_path.suffix.lower()

    if extension not in supported:
        raise ValueError(
            f"Unsupported output extension '{extension}'. "
            "Use .step, .stp, .gltf, or .glb"
        )


def export_gltf_in_metres(
    assembly: cq.Assembly,
    output_path: Path,
) -> None:
    """
    Export an assembly whose internal coordinates are millimetres.

    glTF requires metres. CadQuery's standard assembly exporter changes
    Z-up to glTF's Y-up, but it does not expose the source length unit.
    This mirrors CadQuery's exporter and configures Open CASCADE's
    coordinate-system converter to transform millimetres to metres.
    """
    binary = output_path.suffix.lower() == ".glb"
    original_location = assembly.loc

    try:
        # Match CadQuery's standard glTF axis conversion:
        # right-handed +Z up -> right-handed +Y up.
        assembly.loc *= cq.Location(
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            -90.0,
        )

        _, document = toCAF(
            assembly,
            coloredSTEP=True,
            mesh=True,
            tolerance=GLTF_LINEAR_TOLERANCE_MM,
            angularTolerance=GLTF_ANGULAR_TOLERANCE_RAD,
        )

        writer = RWGltf_CafWriter(
            TCollection_AsciiString(str(output_path)),
            binary,
        )

        converter = writer.ChangeCoordinateSystemConverter()
        converter.SetInputLengthUnit(MM_TO_M)
        converter.SetOutputLengthUnit(GLTF_OUTPUT_UNIT_M)

        succeeded = writer.Perform(
            document,
            TColStd_IndexedDataMapOfStringString(),
            Message_ProgressRange(),
        )
    finally:
        # Do not leave the in-memory assembly rotated if exporting fails.
        assembly.loc = original_location

    if not succeeded:
        raise RuntimeError(
            f"Open CASCADE failed to export '{output_path}'"
        )


def export_assembly(
    assembly: cq.Assembly,
    output_path: Path,
) -> None:
    validate_output_extension(output_path)

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    extension = output_path.suffix.lower()

    if extension in {".gltf", ".glb"}:
        export_gltf_in_metres(
            assembly,
            output_path,
        )
        return

    # STEP stays in the script's native millimetre units.
    assembly.export(str(output_path))


def main() -> None:
    args = parse_arguments()

    if args.every < 1:
        raise ValueError("--every must be at least 1")

    if (
        args.max_sensors is not None
        and args.max_sensors < 1
    ):
        raise ValueError(
            "--max-sensors must be at least 1"
        )

    if args.gap_mm < 0.0:
        raise ValueError("--gap-mm cannot be negative")

    if (
        args.report_projection_above_mm is not None
        and args.report_projection_above_mm < 0.0
    ):
        raise ValueError(
            "--report-projection-above-mm cannot be negative"
        )

    if (
        args.reject_projection_above_mm is not None
        and args.reject_projection_above_mm < 0.0
    ):
        raise ValueError(
            "--reject-projection-above-mm cannot be negative"
        )

    sensors = read_sensor_faces(
        args.input_csv,
        args,
    )

    ordered_sensors = sorted(
        sensors.items(),
        key=lambda item: (
            item[0].subdet,
            item[0].layer_or_disk,
            item[0].detector_type,
            item[0].stack_rawid,
            item[0].sensor_role,
            item[0].sensor_rawid,
        ),
    )

    ordered_sensors = ordered_sensors[:: args.every]

    if args.max_sensors is not None:
        ordered_sensors = ordered_sensors[
            : args.max_sensors
        ]

    if not ordered_sensors:
        raise RuntimeError(
            "No sensors matched the requested filters"
        )

    print(
        f"Found {len(ordered_sensors)} sensors "
        "after applying filters"
    )

    tracker_assembly = cq.Assembly(
        name="D110_T35_OuterTracker_Segmented"
    )

    sensors_built = 0
    segments_built = 0
    sensors_failed = 0

    for key, faces in ordered_sensors:
        divisions = SEGMENTATION_BY_TYPE[
            key.detector_type
        ]

        try:
            if args.mode == "surface":
                segment_shapes, projection_info = (
                    make_surface_segments(
                        faces,
                        divisions,
                        args.gap_mm,
                    )
                )
            else:
                segment_shapes, projection_info = (
                    make_solid_segments(
                        faces,
                        divisions,
                        args.gap_mm,
                    )
                )

            if (
                args.report_projection_above_mm is not None
                and projection_info.maximum_mm
                > args.report_projection_above_mm
            ):
                print(
                    "Projection warning: "
                    f"sensor={key.sensor_rawid}, "
                    f"type={key.detector_type}, "
                    f"maximum={projection_info.maximum_mm:.12g} mm, "
                    f"rms={projection_info.rms_mm:.12g} mm"
                )

            if (
                args.reject_projection_above_mm is not None
                and projection_info.maximum_mm
                > args.reject_projection_above_mm
            ):
                raise ValueError(
                    "Projection correction exceeds limit: "
                    f"{projection_info.maximum_mm:.12g} mm"
                )

            sensor_name = make_sensor_name(key)

            sensor_assembly = cq.Assembly(
                name=sensor_name
            )

            for segment_index, shape in enumerate(
                segment_shapes
            ):
                segment_name = (
                    f"{sensor_name}_"
                    f"column{segment_index:02d}"
                )

                sensor_assembly.add(
                    shape,
                    name=segment_name,
                    color=segment_color(
                        key,
                        segment_index,
                    ),
                )

                segments_built += 1

            tracker_assembly.add(
                sensor_assembly,
                name=sensor_name,
            )

            sensors_built += 1

        except Exception as error:
            sensors_failed += 1

            print(
                "Failed to construct "
                f"sensor={key.sensor_rawid}, "
                f"stack={key.stack_rawid}, "
                f"type={key.detector_type}: "
                f"{error}"
            )

        if sensors_built > 0 and sensors_built % 250 == 0:
            print(
                f"Constructed {sensors_built} sensors "
                f"and {segments_built} segments"
            )

    if sensors_built == 0:
        raise RuntimeError(
            "No valid sensor shapes were constructed"
        )

    print(
        f"Exporting {sensors_built} sensors "
        f"with {segments_built} segments "
        f"to {args.output}"
    )

    export_assembly(
        tracker_assembly,
        args.output,
    )

    print("Export complete")

    if sensors_failed:
        print(
            f"{sensors_failed} sensors could not be constructed"
        )


if __name__ == "__main__":
    main()
