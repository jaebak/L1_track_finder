#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Tuple

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator


CM_TO_MM = 10.0

Point3D = Tuple[float, float, float]
Point2D = Tuple[float, float]
FaceCorners = Dict[int, Point3D]


@dataclass(frozen=True)
class SensorKey:
    stack_rawid: int
    sensor_rawid: int
    subdet: str
    layer_or_disk: int
    sensor_role: str
    detector_type: str


# Colors chosen to resemble the supplied CMS-style r-z view.
COLOR_BY_TYPE = {
    "Ph2SS": "#ff3b1f",
    "Ph2PSS": "#0068b7",
    "Ph2PSP": "#f4b400",
}


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create a CMS tracker r-z projection from the same sensor-surface "
            "CSV used by the CadQuery model."
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
        help="Output plot: .png, .pdf, .svg, or another Matplotlib format",
    )
    parser.add_argument(
        "--z-side",
        choices=("positive", "negative", "both"),
        default="positive",
        help=(
            "Which detector half to draw. For 'negative', -z is displayed "
            "as a positive horizontal coordinate. Default: positive"
        ),
    )
    parser.add_argument(
        "--subdet",
        action="append",
        choices=("TOB", "TID"),
        help="Restrict the plot to one or more subdetectors",
    )
    parser.add_argument(
        "--layer",
        type=int,
        action="append",
        help="Restrict the plot to one or more layers/disks",
    )
    parser.add_argument(
        "--sensor-role",
        action="append",
        choices=("lower", "upper"),
        help="Restrict the plot to one or more sensor roles",
    )
    parser.add_argument(
        "--detector-type",
        action="append",
        choices=("Ph2PSP", "Ph2PSS", "Ph2SS"),
        help="Restrict the plot to one or more detector types",
    )
    parser.add_argument(
        "--phi-center-deg",
        type=float,
        default=None,
        help=(
            "Optionally keep a narrow azimuthal slice centred at this phi. "
            "When omitted, all phi positions are collapsed into r-z"
        ),
    )
    parser.add_argument(
        "--phi-half-width-deg",
        type=float,
        default=5.0,
        help="Half-width of the optional phi slice. Default: 5 degrees",
    )
    parser.add_argument(
        "--deduplicate-mm",
        type=float,
        default=0.5,
        help=(
            "Merge projected modules with nearly identical r-z bounding boxes. "
            "Set to 0 to disable. Default: 0.5 mm"
        ),
    )
    parser.add_argument(
        "--line-width",
        type=float,
        default=0.9,
        help="Sensor outline width. Default: 0.9",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.95,
        help="Sensor line opacity. Default: 0.95",
    )
    parser.add_argument(
        "--z-max-mm",
        type=float,
        default=None,
        help="Optional horizontal-axis maximum in mm",
    )
    parser.add_argument(
        "--r-max-mm",
        type=float,
        default=None,
        help="Optional vertical-axis maximum in mm",
    )
    parser.add_argument(
        "--eta-max",
        type=float,
        default=4.0,
        help="Largest pseudorapidity boundary label. Default: 4.0",
    )
    parser.add_argument(
        "--eta-step",
        type=float,
        default=0.2,
        help="Pseudorapidity boundary-label spacing. Default: 0.2",
    )
    parser.add_argument(
        "--no-eta-axis",
        action="store_true",
        help="Do not draw pseudorapidity ticks along the top/right boundary",
    )
    parser.add_argument(
        "--title",
        default=None,
        help="Optional plot title",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=220,
        help="Raster output resolution. Default: 220 dpi",
    )

    return parser.parse_args()


def row_passes_filters(row: dict[str, str], args: argparse.Namespace) -> bool:
    if row["face"].strip() != "mid_plane":
        return False

    detector_type = row["detector_type"].strip()
    if detector_type not in COLOR_BY_TYPE:
        return False

    if args.subdet is not None and row["subdet"].strip() not in args.subdet:
        return False

    if args.layer is not None and int(row["layer_or_disk"]) not in args.layer:
        return False

    if (
        args.sensor_role is not None
        and row["sensor_role"].strip() not in args.sensor_role
    ):
        return False

    if (
        args.detector_type is not None
        and detector_type not in args.detector_type
    ):
        return False

    return True


def read_mid_planes(
    csv_path: Path,
    args: argparse.Namespace,
) -> dict[SensorKey, FaceCorners]:
    sensors: dict[SensorKey, FaceCorners] = defaultdict(dict)

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

        missing = required_columns - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                "Input CSV is missing required columns: "
                + ", ".join(sorted(missing))
            )

        for row in reader:
            if not row_passes_filters(row, args):
                continue

            key = SensorKey(
                stack_rawid=int(row["stack_rawid"]),
                sensor_rawid=int(row["sensor_rawid"]),
                subdet=row["subdet"].strip(),
                layer_or_disk=int(row["layer_or_disk"]),
                sensor_role=row["sensor_role"].strip(),
                detector_type=row["detector_type"].strip(),
            )

            corner = int(row["corner"])
            sensors[key][corner] = (
                CM_TO_MM * float(row["global_x_cm"]),
                CM_TO_MM * float(row["global_y_cm"]),
                CM_TO_MM * float(row["global_z_cm"]),
            )

    return dict(sensors)


def ordered_corners(corners: FaceCorners, key: SensorKey) -> list[Point3D]:
    expected = {0, 1, 2, 3}
    available = set(corners)
    if available != expected:
        raise ValueError(
            f"Sensor {key.sensor_rawid} contains corners {sorted(available)}, "
            "but corners 0, 1, 2, and 3 are required"
        )
    return [corners[index] for index in range(4)]


def wrapped_angle_difference_deg(first: float, second: float) -> float:
    return abs((first - second + 180.0) % 360.0 - 180.0)


def sensor_centroid(points: Iterable[Point3D]) -> Point3D:
    points_list = list(points)
    count = len(points_list)
    return (
        sum(point[0] for point in points_list) / count,
        sum(point[1] for point in points_list) / count,
        sum(point[2] for point in points_list) / count,
    )


def passes_phi_slice(points: list[Point3D], args: argparse.Namespace) -> bool:
    if args.phi_center_deg is None:
        return True

    x, y, _ = sensor_centroid(points)
    phi_deg = math.degrees(math.atan2(y, x))
    return (
        wrapped_angle_difference_deg(phi_deg, args.phi_center_deg)
        <= args.phi_half_width_deg
    )


def project_point_to_rz(point: Point3D, z_side: str) -> Point2D:
    x, y, z = point
    radius = math.hypot(x, y)

    if z_side == "negative":
        horizontal = -z
    else:
        horizontal = z

    return horizontal, radius


def sensor_is_on_requested_side(points: list[Point3D], z_side: str) -> bool:
    _, _, centroid_z = sensor_centroid(points)
    if z_side == "positive":
        return centroid_z >= 0.0
    if z_side == "negative":
        return centroid_z <= 0.0
    return True


def deduplication_key(points: list[Point2D], tolerance_mm: float) -> tuple[int, ...] | None:
    if tolerance_mm <= 0.0:
        return None

    z_values = [point[0] for point in points]
    r_values = [point[1] for point in points]
    values = (
        min(z_values),
        max(z_values),
        min(r_values),
        max(r_values),
        sum(z_values) / len(z_values),
        sum(r_values) / len(r_values),
    )
    return tuple(round(value / tolerance_mm) for value in values)


def padded_limit(maximum: float, explicit: float | None, fraction: float) -> float:
    if explicit is not None:
        if explicit <= 0.0:
            raise ValueError("Axis maxima must be positive")
        return explicit
    if maximum <= 0.0:
        return 1.0
    return maximum * (1.0 + fraction)


def add_eta_boundary_axis(
    ax: plt.Axes,
    z_min: float,
    z_max: float,
    r_max: float,
    eta_max: float,
    eta_step: float,
) -> None:
    """
    Label constant-pseudorapidity rays where they leave the top/right frame.

    For positive z and r:
        z / r = sinh(eta)
    """
    if z_min < 0.0:
        return
    if eta_step <= 0.0 or eta_max < 0.0:
        raise ValueError("eta-step must be positive and eta-max non-negative")

    tick_dx = 0.018 * (z_max - z_min)
    tick_dy = 0.045 * r_max
    label_dx = 0.026 * (z_max - z_min)
    label_dy = 0.064 * r_max

    count = int(math.floor(eta_max / eta_step + 1.0e-9))
    eta_values = [index * eta_step for index in range(count + 1)]

    for eta in eta_values:
        if eta == 0.0:
            top_z = 0.0
        else:
            top_z = r_max * math.sinh(eta)

        if z_min <= top_z <= z_max:
            ax.plot(
                [top_z, top_z + tick_dx],
                [r_max, r_max + tick_dy],
                color="black",
                linewidth=0.55,
                clip_on=False,
            )
            ax.text(
                top_z,
                r_max + label_dy,
                f"{eta:.1f}",
                ha="center",
                va="bottom",
                fontsize=8,
                clip_on=False,
            )
            continue

        if eta <= 0.0:
            continue

        right_r = z_max / math.sinh(eta)
        if 0.0 <= right_r <= r_max:
            ax.plot(
                [z_max, z_max + tick_dx],
                [right_r, right_r + 0.35 * tick_dy],
                color="black",
                linewidth=0.55,
                clip_on=False,
            )
            ax.text(
                z_max + label_dx,
                right_r,
                f"{eta:.1f}",
                ha="left",
                va="center",
                fontsize=8,
                clip_on=False,
            )

    ax.text(
        z_max + 0.045 * (z_max - z_min),
        -0.03 * r_max,
        r"$\eta$",
        ha="center",
        va="top",
        fontsize=15,
        clip_on=False,
    )


def build_plot(args: argparse.Namespace) -> None:
    sensors = read_mid_planes(args.input_csv, args)
    if not sensors:
        raise RuntimeError("No mid-plane sensors matched the requested filters")

    projected_items: list[tuple[SensorKey, list[Point2D]]] = []
    seen: set[tuple[str, tuple[int, ...]]] = set()

    for key, corners in sorted(
        sensors.items(),
        key=lambda item: (
            item[0].subdet,
            item[0].layer_or_disk,
            item[0].detector_type,
            item[0].stack_rawid,
            item[0].sensor_role,
            item[0].sensor_rawid,
        ),
    ):
        points = ordered_corners(corners, key)

        if not sensor_is_on_requested_side(points, args.z_side):
            continue
        if not passes_phi_slice(points, args):
            continue

        projected = [project_point_to_rz(point, args.z_side) for point in points]

        signature = deduplication_key(projected, args.deduplicate_mm)
        if signature is not None:
            typed_signature = (key.detector_type, signature)
            if typed_signature in seen:
                continue
            seen.add(typed_signature)

        projected_items.append((key, projected))

    if not projected_items:
        raise RuntimeError("No sensors remained after z-side/phi-slice filtering")

    all_z = [point[0] for _, points in projected_items for point in points]
    all_r = [point[1] for _, points in projected_items for point in points]

    if args.z_side == "both":
        z_extent = max(abs(min(all_z)), abs(max(all_z)))
        if args.z_max_mm is not None:
            z_extent = args.z_max_mm
        z_min = -1.03 * z_extent
        z_max = 1.03 * z_extent
    else:
        z_min = 0.0
        z_max = padded_limit(max(all_z), args.z_max_mm, 0.04)

    r_max = padded_limit(max(all_r), args.r_max_mm, 0.04)

    fig, ax = plt.subplots(figsize=(14.0, 5.2))

    for key, points in projected_items:
        closed = points + [points[0]]
        ax.plot(
            [point[0] for point in closed],
            [point[1] for point in closed],
            color=COLOR_BY_TYPE[key.detector_type],
            linewidth=args.line_width,
            alpha=args.alpha,
            solid_capstyle="butt",
            solid_joinstyle="miter",
        )

    ax.set_xlim(z_min, z_max)
    ax.set_ylim(0.0, r_max)
    ax.set_xlabel("z [mm]", loc="right")
    ax.set_ylabel("r [mm]", loc="top")

    if args.title:
        ax.set_title(args.title)

    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.tick_params(which="major", direction="in", length=7, width=0.7)
    ax.tick_params(which="minor", direction="in", length=4, width=0.55)
    ax.grid(False)

    legend_handles = [
        Line2D(
            [0],
            [0],
            color=COLOR_BY_TYPE[detector_type],
            linewidth=2.0,
            label=detector_type,
        )
        for detector_type in COLOR_BY_TYPE
        if any(key.detector_type == detector_type for key, _ in projected_items)
    ]
    ax.legend(
        handles=legend_handles,
        loc="upper right",
        frameon=False,
        ncol=len(legend_handles),
        fontsize=8,
    )

    if not args.no_eta_axis and args.z_side != "both":
        add_eta_boundary_axis(
            ax=ax,
            z_min=z_min,
            z_max=z_max,
            r_max=r_max,
            eta_max=args.eta_max,
            eta_step=args.eta_step,
        )

    fig.subplots_adjust(
        left=0.075,
        right=0.87 if not args.no_eta_axis and args.z_side != "both" else 0.97,
        bottom=0.14,
        top=0.82 if not args.no_eta_axis and args.z_side != "both" else 0.93,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=args.dpi, facecolor="white")
    plt.close(fig)

    print(
        f"Plotted {len(projected_items)} unique projected sensors "
        f"to {args.output}"
    )


def main() -> None:
    args = parse_arguments()

    if args.phi_half_width_deg <= 0.0:
        raise ValueError("--phi-half-width-deg must be positive")
    if args.deduplicate_mm < 0.0:
        raise ValueError("--deduplicate-mm cannot be negative")
    if args.line_width <= 0.0:
        raise ValueError("--line-width must be positive")
    if not 0.0 < args.alpha <= 1.0:
        raise ValueError("--alpha must be in the interval (0, 1]")
    if args.dpi < 1:
        raise ValueError("--dpi must be positive")

    build_plot(args)


if __name__ == "__main__":
    main()
