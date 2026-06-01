#!/usr/bin/env python
"""Classify coastline line segments by nearest GCC coastal point."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Any

try:
    import geopandas as gpd
    import pandas as pd
    from pyproj import CRS
    from shapely.geometry import GeometryCollection, LineString, MultiLineString, MultiPoint
    from shapely.ops import voronoi_diagram
    from shapely.validation import make_valid
except ImportError as exc:  # pragma: no cover - exercised before dependencies exist
    raise SystemExit(
        "Missing geospatial dependencies. Install them with: "
        "python -m pip install -r requirements.txt"
    ) from exc


CLASS_LIKELY = "habitat protection likely relevant"
CLASS_POSSIBLE = "habitat protection possibly relevant"
CLASS_UNLIKELY = "habitat protection unlikely to matter"
CLASS_INSUFFICIENT = "insufficient data"

FIELD_CANDIDATES = {
    "coast_type": [
        "coast_type",
        "coastal_type",
        "shore_type",
        "shoreline_type",
        "typology",
        "class",
        "classification",
    ],
    "slope": ["slope", "beach_slope", "profile_slope"],
    "elevation": ["elevation", "elev", "height", "z", "max_elevation"],
    "wave": ["wave", "hs", "h_s", "significant_wave", "water_level", "surge"],
    "armoring": ["armoring", "armouring", "armor", "armour", "defense", "defence", "hard"],
    "vegetation": ["vegetation", "vegetated", "veg", "forest", "mangrove", "marsh"],
    "sandy": ["sandy", "sand", "beach"],
    "rocky": ["rocky", "rock", "cliff"],
}

COMPATIBLE_TERMS = {
    "barrier",
    "beach",
    "delta",
    "dune",
    "estuary",
    "lagoon",
    "mangrove",
    "marsh",
    "mud",
    "muddy",
    "reef",
    "sand",
    "sandy",
    "seagrass",
    "tidal flat",
    "vegetated",
    "wetland",
}

UNLIKELY_TERMS = {
    "armored",
    "armoured",
    "artificial",
    "cliff",
    "defended",
    "engineered",
    "hard",
    "rock",
    "rocky",
    "seawall",
    "steep",
}

EXPOSURE_TERMS = {"exposed", "high", "moderate", "storm", "surge", "wave"}


@dataclass(frozen=True)
class FieldMap:
    coast_type: str | None
    slope: str | None
    elevation: str | None
    wave: str | None
    armoring: str | None
    vegetation: str | None
    sandy: str | None
    rocky: str | None

    def selected_columns(self) -> list[str]:
        return list(dict.fromkeys(column for column in self.__dict__.values() if column is not None))


@dataclass(frozen=True)
class Classification:
    label: str
    reason: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Clip a coastline line dataset to an AOI, split it by Voronoi cells "
            "derived from GCC coastal points, and classify each line segment."
        )
    )
    parser.add_argument("--aoi", required=True, help="Polygon vector AOI path.")
    parser.add_argument("--aoi-layer", help="Optional AOI layer name.")
    parser.add_argument("--coastline", required=True, help="Coastline line vector path.")
    parser.add_argument("--coastline-layer", help="Optional coastline layer name.")
    parser.add_argument("--gcc", required=True, help="GCC point/transect vector path.")
    parser.add_argument("--gcc-layer", help="Optional GCC layer name.")
    parser.add_argument("--gcc-x-field", help="Longitude/X field for CSV GCC input. Defaults to lon/longitude/x.")
    parser.add_argument("--gcc-y-field", help="Latitude/Y field for CSV GCC input. Defaults to lat/latitude/y.")
    parser.add_argument("--gcc-crs", default="EPSG:4326", help="CRS for CSV GCC coordinates.")
    parser.add_argument("--output", required=True, help="Output vector path, preferably .gpkg.")
    parser.add_argument("--output-layer", default="shoreline_segments")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite output path if it exists.")

    parser.add_argument("--gcc-id-field", help="GCC identifier field to carry into the output.")
    parser.add_argument("--coast-type-field")
    parser.add_argument("--slope-field")
    parser.add_argument("--elevation-field")
    parser.add_argument("--wave-field")
    parser.add_argument("--armoring-field")
    parser.add_argument("--vegetation-field")
    parser.add_argument("--sandy-field")
    parser.add_argument("--rocky-field")

    parser.add_argument(
        "--metric-crs",
        help="Metric CRS for Voronoi, clipping, length, and distance. Defaults to local UTM or EPSG:8857.",
    )
    parser.add_argument(
        "--gcc-buffer-m",
        type=float,
        default=50_000,
        help="Include GCC points within this buffer around the AOI before tessellation.",
    )
    parser.add_argument(
        "--max-gcc-distance-m",
        type=float,
        default=10_000,
        help="Mark segments farther than this from their GCC point as insufficient data. Use 0 to disable.",
    )
    parser.add_argument(
        "--high-elevation-threshold-m",
        type=float,
        default=10,
        help="Elevation at or above this is treated as a strong unlikely signal.",
    )
    parser.add_argument(
        "--low-elevation-threshold-m",
        type=float,
        default=5,
        help="Elevation at or below this is treated as a low-coast signal.",
    )
    parser.add_argument(
        "--slope-units",
        choices=["ratio", "degrees"],
        default="ratio",
        help="Units for the selected slope field.",
    )
    parser.add_argument(
        "--steep-slope-threshold",
        type=float,
        default=0.05,
        help="Slope at or above this is treated as steep. Default assumes rise/run ratio.",
    )
    parser.add_argument(
        "--gentle-slope-threshold",
        type=float,
        default=0.01,
        help="Slope at or below this is treated as gentle. Default assumes rise/run ratio.",
    )
    parser.add_argument(
        "--wave-exposure-threshold",
        type=float,
        default=0.5,
        help="Numeric wave/exposure values at or above this are treated as exposed.",
    )
    return parser.parse_args()


def read_vector(path: str, layer: str | None) -> gpd.GeoDataFrame:
    kwargs = {"layer": layer} if layer else {}
    gdf = gpd.read_file(path, **kwargs)
    if gdf.empty:
        raise SystemExit(f"No features found in {path}.")
    if gdf.crs is None:
        raise SystemExit(f"{path} has no CRS. Define one before running this utility.")
    gdf = gdf[gdf.geometry.notna()].copy()
    gdf["geometry"] = gdf.geometry.map(lambda geom: make_valid(geom) if not geom.is_valid else geom)
    return gdf


def resolve_coordinate_field(columns: list[str], explicit: str | None, candidates: list[str], axis: str) -> str:
    if explicit:
        if explicit not in columns:
            raise SystemExit(f"Requested GCC {axis} field '{explicit}' was not found.")
        return explicit

    lower_to_column = {column.lower(): column for column in columns}
    for candidate in candidates:
        if candidate.lower() in lower_to_column:
            return lower_to_column[candidate.lower()]
    raise SystemExit(
        f"Could not detect GCC {axis} coordinate field. "
        f"Pass --gcc-{axis}-field explicitly."
    )


def read_gcc(path: str, layer: str | None, args: argparse.Namespace) -> gpd.GeoDataFrame:
    if Path(path).suffix.lower() != ".csv":
        return read_vector(path, layer)

    table = pd.read_csv(path)
    if table.empty:
        raise SystemExit(f"No rows found in {path}.")
    columns = list(table.columns)
    x_field = resolve_coordinate_field(columns, args.gcc_x_field, ["lon", "longitude", "x"], "x")
    y_field = resolve_coordinate_field(columns, args.gcc_y_field, ["lat", "latitude", "y"], "y")
    table = table.dropna(subset=[x_field, y_field]).copy()
    if table.empty:
        raise SystemExit(f"No GCC rows in {path} have both {x_field} and {y_field}.")
    geometry = gpd.points_from_xy(table[x_field], table[y_field])
    return gpd.GeoDataFrame(table, geometry=geometry, crs=args.gcc_crs)


def require_geometry(gdf: gpd.GeoDataFrame, allowed: set[str], name: str) -> None:
    geom_types = set(gdf.geometry.geom_type.dropna().unique())
    invalid = geom_types - allowed
    if invalid:
        raise SystemExit(f"{name} must contain {sorted(allowed)} geometries; found {sorted(geom_types)}.")


def dissolve_aoi(aoi: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    polygons = {"Polygon", "MultiPolygon"}
    require_geometry(aoi, polygons, "AOI")
    if hasattr(aoi.geometry, "union_all"):
        dissolved = aoi.geometry.union_all()
    else:
        dissolved = aoi.geometry.unary_union
    return gpd.GeoDataFrame(geometry=[dissolved], crs=aoi.crs)


def choose_metric_crs(aoi: gpd.GeoDataFrame, explicit: str | None) -> CRS:
    if explicit:
        return CRS.from_user_input(explicit)

    aoi_wgs84 = aoi.to_crs("EPSG:4326")
    minx, miny, maxx, maxy = aoi_wgs84.total_bounds
    if (maxx - minx) <= 12 and (maxy - miny) <= 12:
        estimated = aoi_wgs84.estimate_utm_crs()
        if estimated is not None:
            return CRS.from_user_input(estimated)
    return CRS.from_epsg(8857)  # Equal Earth; works for broad AOIs and global runs.


def resolve_field(gdf: gpd.GeoDataFrame, explicit: str | None, candidates: list[str]) -> str | None:
    columns = [column for column in gdf.columns if column != gdf.geometry.name]
    if explicit:
        if explicit not in columns:
            raise SystemExit(f"Requested field '{explicit}' was not found.")
        return explicit

    lower_to_column = {column.lower(): column for column in columns}
    for candidate in candidates:
        if candidate.lower() in lower_to_column:
            return lower_to_column[candidate.lower()]
    for candidate in candidates:
        for column in columns:
            if candidate.lower() in column.lower():
                return column
    return None


def build_field_map(gcc: gpd.GeoDataFrame, args: argparse.Namespace) -> FieldMap:
    return FieldMap(
        coast_type=resolve_field(gcc, args.coast_type_field, FIELD_CANDIDATES["coast_type"]),
        slope=resolve_field(gcc, args.slope_field, FIELD_CANDIDATES["slope"]),
        elevation=resolve_field(gcc, args.elevation_field, FIELD_CANDIDATES["elevation"]),
        wave=resolve_field(gcc, args.wave_field, FIELD_CANDIDATES["wave"]),
        armoring=resolve_field(gcc, args.armoring_field, FIELD_CANDIDATES["armoring"]),
        vegetation=resolve_field(gcc, args.vegetation_field, FIELD_CANDIDATES["vegetation"]),
        sandy=resolve_field(gcc, args.sandy_field, FIELD_CANDIDATES["sandy"]),
        rocky=resolve_field(gcc, args.rocky_field, FIELD_CANDIDATES["rocky"]),
    )


def is_missing(value: Any) -> bool:
    try:
        return bool(pd.isna(value))
    except TypeError:
        return value is None


def as_float(value: Any) -> float | None:
    if is_missing(value):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def text_contains(value: Any, terms: set[str]) -> bool:
    if is_missing(value):
        return False
    text = str(value).lower().replace("_", " ").replace("-", " ")
    return any(term in text for term in terms)


def truthy(value: Any) -> bool:
    if is_missing(value):
        return False
    if isinstance(value, bool):
        return value
    numeric = as_float(value)
    if numeric is not None:
        return numeric > 0
    return str(value).strip().lower() in {"true", "yes", "y", "1", "present", "presence"}


def classify_gcc_row(row: pd.Series, fields: FieldMap, args: argparse.Namespace) -> Classification:
    evidence = 0
    reasons: list[str] = []
    unlikely = False
    compatible = False
    low_or_gentle = False
    exposed = False

    text_fields = [fields.coast_type, fields.armoring, fields.vegetation, fields.sandy, fields.rocky]
    for field in text_fields:
        if field is None or is_missing(row.get(field)):
            continue
        evidence += 1
        value = row[field]
        if text_contains(value, UNLIKELY_TERMS):
            unlikely = True
            reasons.append(f"{field} indicates steep/rocky/armored context")
        if text_contains(value, COMPATIBLE_TERMS):
            compatible = True
            reasons.append(f"{field} indicates habitat-compatible coastal setting")
        if text_contains(value, EXPOSURE_TERMS):
            exposed = True
            reasons.append(f"{field} indicates hazard exposure")

    if fields.vegetation and truthy(row.get(fields.vegetation)):
        compatible = True
        reasons.append(f"{fields.vegetation} indicates vegetation presence")
    if fields.sandy and truthy(row.get(fields.sandy)):
        compatible = True
        reasons.append(f"{fields.sandy} indicates sandy/beach setting")
    if fields.rocky and truthy(row.get(fields.rocky)):
        unlikely = True
        reasons.append(f"{fields.rocky} indicates rocky context")
    if fields.armoring and truthy(row.get(fields.armoring)):
        unlikely = True
        reasons.append(f"{fields.armoring} indicates hard defense or armoring")

    slope = as_float(row.get(fields.slope)) if fields.slope else None
    if slope is not None:
        evidence += 1
        steep_threshold = args.steep_slope_threshold
        gentle_threshold = args.gentle_slope_threshold
        if args.slope_units == "degrees" and steep_threshold == 0.05:
            steep_threshold = 3.0
        if args.slope_units == "degrees" and gentle_threshold == 0.01:
            gentle_threshold = 1.0
        if slope >= steep_threshold:
            unlikely = True
            reasons.append(f"{fields.slope}={slope:g} meets steep threshold")
        elif slope <= gentle_threshold:
            low_or_gentle = True
            reasons.append(f"{fields.slope}={slope:g} meets gentle threshold")

    elevation = as_float(row.get(fields.elevation)) if fields.elevation else None
    if elevation is not None:
        evidence += 1
        if elevation >= args.high_elevation_threshold_m:
            unlikely = True
            reasons.append(f"{fields.elevation}={elevation:g}m meets high-elevation threshold")
        elif elevation <= args.low_elevation_threshold_m:
            low_or_gentle = True
            reasons.append(f"{fields.elevation}={elevation:g}m meets low-elevation threshold")

    wave = as_float(row.get(fields.wave)) if fields.wave else None
    if wave is not None:
        evidence += 1
        if wave >= args.wave_exposure_threshold:
            exposed = True
            reasons.append(f"{fields.wave}={wave:g} meets exposure threshold")
    elif fields.wave and text_contains(row.get(fields.wave), EXPOSURE_TERMS):
        evidence += 1
        exposed = True
        reasons.append(f"{fields.wave} indicates exposure")

    if evidence == 0:
        return Classification(CLASS_INSUFFICIENT, "no usable classifier fields found")
    if unlikely:
        return Classification(CLASS_UNLIKELY, "; ".join(dict.fromkeys(reasons)))
    if exposed and low_or_gentle and compatible:
        return Classification(CLASS_LIKELY, "; ".join(dict.fromkeys(reasons)))
    if (exposed and (low_or_gentle or compatible)) or (low_or_gentle and compatible):
        return Classification(CLASS_POSSIBLE, "; ".join(dict.fromkeys(reasons)))
    return Classification(CLASS_POSSIBLE, "partial evidence only; conservative first-pass classification")


def line_parts(geometry: Any) -> list[LineString | MultiLineString]:
    if geometry is None or geometry.is_empty:
        return []
    if isinstance(geometry, (LineString, MultiLineString)):
        return [geometry]
    if isinstance(geometry, GeometryCollection):
        parts: list[LineString | MultiLineString] = []
        for part in geometry.geoms:
            parts.extend(line_parts(part))
        return parts
    return []


def build_voronoi_cells(gcc: gpd.GeoDataFrame, envelope: Any) -> gpd.GeoDataFrame:
    if len(gcc) == 1:
        return gpd.GeoDataFrame({"_gcc_internal_id": [gcc.iloc[0]["_gcc_internal_id"]]}, geometry=[envelope], crs=gcc.crs)

    diagram = voronoi_diagram(MultiPoint(list(gcc.geometry)), envelope=envelope, edges=False)
    polygons = list(diagram.geoms) if hasattr(diagram, "geoms") else [diagram]
    cells = gpd.GeoDataFrame({"_cell_id": range(len(polygons))}, geometry=polygons, crs=gcc.crs)
    cells = gpd.clip(cells, gpd.GeoDataFrame(geometry=[envelope], crs=gcc.crs))
    cells = gpd.sjoin_nearest(
        cells,
        gcc[["_gcc_internal_id", "geometry"]],
        how="left",
        distance_col="_cell_seed_distance_m",
    )
    return cells.drop(columns=["index_right"]).drop_duplicates("_cell_id")


def explode_to_linework(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    rows: list[dict[str, Any]] = []
    for _, row in gdf.iterrows():
        attrs = row.drop(labels=[gdf.geometry.name]).to_dict()
        for part in line_parts(row.geometry):
            rows.append({**attrs, "geometry": part})
    return gpd.GeoDataFrame(rows, geometry="geometry", crs=gdf.crs)


def main() -> None:
    args = parse_args()
    output_path = Path(args.output)
    if output_path.exists():
        if not args.overwrite:
            raise SystemExit(f"{output_path} already exists. Use --overwrite to replace it.")
        output_path.unlink()

    aoi = dissolve_aoi(read_vector(args.aoi, args.aoi_layer))
    coastline = read_vector(args.coastline, args.coastline_layer)
    gcc = read_gcc(args.gcc, args.gcc_layer, args)
    require_geometry(coastline, {"LineString", "MultiLineString"}, "Coastline")
    require_geometry(gcc, {"Point"}, "GCC")

    metric_crs = choose_metric_crs(aoi, args.metric_crs)
    aoi_m = aoi.to_crs(metric_crs)
    coastline_m = coastline.to_crs(metric_crs)
    gcc_m = gcc.to_crs(metric_crs).copy()

    coastline_m = coastline_m.reset_index(drop=True)
    coastline_m["_source_line_id"] = coastline_m.index + 1
    clipped_coast = gpd.clip(coastline_m[["_source_line_id", "geometry"]], aoi_m)
    clipped_coast = explode_to_linework(clipped_coast)
    if clipped_coast.empty:
        raise SystemExit("No coastline linework intersects the AOI.")

    aoi_buffer = aoi_m.geometry.iloc[0].buffer(args.gcc_buffer_m)
    gcc_m = gcc_m[gcc_m.geometry.intersects(aoi_buffer)].copy().reset_index(drop=True)
    if gcc_m.empty:
        raise SystemExit("No GCC points found within the AOI buffer. Increase --gcc-buffer-m or check inputs.")

    gcc_m["_gcc_internal_id"] = gcc_m.index + 1
    if args.gcc_id_field:
        if args.gcc_id_field not in gcc_m.columns:
            raise SystemExit(f"Requested GCC id field '{args.gcc_id_field}' was not found.")
        gcc_m["gcc_id"] = gcc_m[args.gcc_id_field].astype(str)
    else:
        gcc_m["gcc_id"] = gcc_m["_gcc_internal_id"].astype(str)

    fields = build_field_map(gcc_m, args)
    classifications = gcc_m.apply(lambda row: classify_gcc_row(row, fields, args), axis=1)
    gcc_m["classification"] = [item.label for item in classifications]
    gcc_m["classification_reason"] = [item.reason for item in classifications]

    cells = build_voronoi_cells(gcc_m[["_gcc_internal_id", "geometry"]], aoi_buffer)
    cell_attrs = cells.merge(
        gcc_m[
            [
                "_gcc_internal_id",
                "gcc_id",
                "classification",
                "classification_reason",
                *fields.selected_columns(),
            ]
        ],
        on="_gcc_internal_id",
        how="left",
    )

    segmented = gpd.overlay(
        clipped_coast,
        cell_attrs.drop(columns=["_cell_seed_distance_m"], errors="ignore"),
        how="intersection",
        keep_geom_type=False,
    )
    segmented = explode_to_linework(segmented)
    if segmented.empty:
        raise SystemExit("Voronoi segmentation produced no coastline linework.")

    gcc_points = dict(zip(gcc_m["_gcc_internal_id"], gcc_m.geometry))
    segmented["gcc_distance_m"] = segmented.apply(
        lambda row: row.geometry.distance(gcc_points[row["_gcc_internal_id"]]),
        axis=1,
    )
    if args.max_gcc_distance_m > 0:
        too_far = segmented["gcc_distance_m"] > args.max_gcc_distance_m
        segmented.loc[too_far, "classification"] = CLASS_INSUFFICIENT
        segmented.loc[too_far, "classification_reason"] = (
            "nearest GCC point is farther than "
            + str(args.max_gcc_distance_m)
            + " m from this coastline segment"
        )

    segmented["length_m"] = segmented.geometry.length
    segmented = segmented.reset_index(drop=True)
    segmented["segment_id"] = segmented.index + 1

    output_columns = [
        "segment_id",
        "_source_line_id",
        "gcc_id",
        "gcc_distance_m",
        "length_m",
        "classification",
        "classification_reason",
        *fields.selected_columns(),
        "geometry",
    ]
    output = segmented[[column for column in output_columns if column in segmented.columns]].to_crs(aoi.crs)

    driver = "GPKG" if output_path.suffix.lower() == ".gpkg" else None
    write_kwargs = {"layer": args.output_layer} if driver == "GPKG" else {}
    if driver:
        write_kwargs["driver"] = driver
    output.to_file(output_path, **write_kwargs)

    print(f"Wrote {len(output)} classified shoreline segments to {output_path}")
    print(f"Metric CRS used for processing: {metric_crs.to_string()}")
    print("Detected classifier fields:")
    for key, value in fields.__dict__.items():
        print(f"  {key}: {value or '(none)'}")


if __name__ == "__main__":
    main()
