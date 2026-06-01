# Coastal Natural Protection Model

This repository is a workspace for a global, coarse-screening model of where coastal habitats may reduce shoreline risk and where their loss or restoration could change exposure to coastal hazards.

The model is intended to improve on simple coastal vulnerability screening by asking a more specific question:

> Where could natural habitat plausibly reduce coastal risk, and how would risk change if that habitat were removed or restored?

The first version should be globally applicable, shoreline-segment based, and driven by public datasets.

## Model Phases

### 1. Need: Where Is Natural Coastal Protection Relevant?

The first phase screens shoreline segments for whether habitat-based protection could plausibly matter at all.

Some shorelines may not benefit much from habitat protection because they are already rocky, steep, high-elevation, armored, or otherwise not meaningfully exposed to the hazard pathway that habitat can reduce.

Potential screening inputs:

- [Global Coastal Characteristics, GCC](https://essd.copernicus.org/articles/16/3433/2024/index.html): global coastal transects with geophysical, hydrometeorological, and socioeconomic indicators.
- [GCC Zenodo dataset](https://zenodo.org/doi/10.5281/zenodo.8200199): downloadable GCC tables.
- [Coastal Hazard Wheel, CHW](https://www.coastalhazardwheel.org/): global coastal classification and coastal hazard screening framework.
- [CHW coastal classification](https://coastalhazardwheel.org/coastal-classification/): classification of coastal type, exposure, and hazard context.
- [GEBCO bathymetry](https://www.gebco.net/data-products/gridded-bathymetry-data): global bathymetry and elevation grid.
- [Copernicus DEM](https://dataspace.copernicus.eu/explore-data/data-collections/copernicus-contributing-missions/collections-description/COP-DEM): global elevation data.
- [CoDEC / GTSM extreme sea levels](https://research.vu.nl/en/datasets/codec-dataset-data-underlying-the-paper-a-high-resolution-global-/): tides, storm surge, and extreme sea-level estimates.
- [ERA5 wave and climate data](https://cds.climate.copernicus.eu/): global reanalysis data for wave and meteorological forcing.

Possible output:

- Shoreline segments classified as:
  - habitat protection likely relevant
  - habitat protection possibly relevant
  - habitat protection unlikely to matter
  - insufficient data

### 2. Habitat: What Protective Habitat Exists Or Could Exist There?

The second phase screens shoreline segments for current protective habitats and potential habitat opportunity.

Habitats in scope:

- mangroves
- saltmarshes
- coral reefs
- seagrasses
- dunes and sandy beaches

Current habitat datasets:

- [Global Mangrove Watch](https://www.wetlands.org/coasts-and-deltas/global-mangrove-watch/): mangrove extent and change.
- [Allen Coral Atlas](https://developers.google.com/earth-engine/datasets/catalog/ACA_reef_habitat_v2_0): global shallow coral reef geomorphic and benthic habitat maps.
- [UNEP-WCMC Global Distribution of Saltmarshes](https://wesr-search.unep.org/ckan/dataset/mapx-global-distribution-of-saltmarshes--2017-): global saltmarsh occurrence and extent.
- [UNEP-WCMC Global Distribution of Seagrasses](https://www.unep.org/resources/publication/global-distribution-seagrasses): global seagrass point and polygon occurrence data.
- [CoastSat](https://github.com/kvos/CoastSat): open-source shoreline extraction from Landsat and Sentinel-2, useful for sandy beaches and local shoreline change workflows.
- GCC and CHW sandy / vegetated coast indicators can serve as coarse global proxies where no better beach or dune layer is available.

Potential habitat screening could use simple ecological plausibility rules, such as:

- mangroves: tropical/subtropical, low-energy, intertidal, suitable temperature and tidal setting
- saltmarsh: temperate or subtropical sheltered intertidal settings
- coral reefs: warm shallow marine areas with existing or potential reef context
- seagrass: shallow nearshore marine areas with suitable substrate/light assumptions
- dunes/beaches: sandy coastlines with beach or sediment supply indicators

Possible output:

- shoreline segments with current protective habitat
- shoreline segments with plausible restoration or expansion opportunity
- habitat type, width, and position relative to shoreline

### 3. Mechanism: What Hazard Pathway Could Habitat Reduce?

The third phase links each habitat type to the hazard pathway it can plausibly reduce.

Examples:

| Habitat | Likely Protective Mechanism |
| --- | --- |
| Mangroves | Wave attenuation, flow resistance, erosion reduction, some surge/friction effects |
| Saltmarshes | Wave attenuation, surge/friction effects, erosion reduction |
| Coral reefs | Offshore wave breaking and wave-energy reduction |
| Seagrasses | Wave/current attenuation in shallow water, sediment stabilization |
| Dunes/beaches | Wave runup buffering, erosion buffering, overwash/flood pathway reduction |

This phase should also recognize that habitat benefits are event-dependent. A habitat may reduce wave energy or flood depth for frequent or moderate events, but may have less effect during extreme events that overwhelm the system.

Possible event classes:

- frequent event
- moderate storm event
- extreme storm event
- sea-level-rise scenario

### 4. Counterfactual: How Much Risk Changes If Habitat Is Removed Or Restored?

The fourth phase compares scenarios.

Core scenarios:

- **current habitat**
- **habitat removed or degraded**
- **habitat restored or expanded**

The first version should report screening metrics, not definitive avoided-damage estimates.

Possible outputs:

- change in wave height at shoreline
- change in wave energy at shoreline
- change in relative exposure score
- shoreline length where habitat removal would increase risk
- shoreline length where restoration could reduce risk
- later: people or assets potentially benefiting

Population and assets are not part of the first model concept, but can be added later using datasets such as [GHSL](https://data.jrc.ec.europa.eu/collection/ghsl/) or [WorldPop](https://www.worldpop.org/).

## Simplest Viable Global Model

The simplest useful global model is a shoreline-segment screening model with habitat attenuation factors.

## Coastline Classification Utility

The first utility creates classified shoreline linework for a polygon AOI using local source data:

- polygon AOI vector
- coastline line vector
- GCC point/transect vector

The script clips the coastline to the AOI, builds Voronoi cells from GCC points in a buffered AOI, intersects the clipped coastline with those cells, and assigns each resulting line segment to the nearest GCC point. The output is a line vector with a conservative first-pass classification of whether habitat-based coastal protection is physically relevant.

Install geospatial dependencies:

```bash
python -m pip install -r requirements.txt
```

Run the utility:

```bash
python scripts/classify_coastline.py \
  --aoi data/aoi/example_aoi.geojson \
  --coastline data/raw/coastline/coastline.gpkg \
  --gcc data/raw/gcc/gcc_points.gpkg \
  --output outputs/classified_coastline.gpkg \
  --overwrite
```

GCC can be supplied as a point vector file or as a CSV with coordinate columns:

```bash
--gcc data/raw/gcc/GCC_geophysical.csv
--gcc-x-field lon
--gcc-y-field lat
--gcc-crs EPSG:4326
```

Optional layer arguments are available for multi-layer files:

```bash
--aoi-layer AOI_LAYER
--coastline-layer COASTLINE_LAYER
--gcc-layer GCC_LAYER
--output-layer shoreline_segments
```

If the GCC field names are not detected correctly, pass explicit field mappings:

```bash
--gcc-id-field transect_id
--coast-type-field coastal_type
--slope-field slope
--elevation-field elevation
--wave-field wave_height
--armoring-field hard_defense
--vegetation-field vegetation
--sandy-field sandy
--rocky-field rocky
```

The output includes:

- `segment_id`
- source coastline line id
- nearest `gcc_id`
- `gcc_distance_m`
- `length_m`
- `classification`
- `classification_reason`
- detected or explicitly mapped GCC classifier fields

Classification values:

- `habitat protection likely relevant`
- `habitat protection possibly relevant`
- `habitat protection unlikely to matter`
- `insufficient data`

This is a screening classifier, not a hydrodynamic model. It classifies physical relevance of habitat-based protection, not confirmed habitat presence or avoided damages.

### Unit Of Analysis

Use shoreline segments that look good on a 2D map.

Each segment should have:

- geometry
- representative shoreline point or transect
- coastal type
- elevation/slope context
- hazard exposure context
- current habitat presence
- potential habitat presence
- scenario results

### Basic Workflow

1. Create or load global shoreline segments.
2. Join each segment to nearby GCC transects and CHW coastal classes.
3. Classify whether habitat protection is relevant.
4. Overlay current habitat datasets.
5. Estimate habitat width and position relative to the shoreline.
6. Apply simple habitat-specific attenuation relationships.
7. Compare current, removed, and restored habitat scenarios.
8. Export shoreline-segment results for mapping.

### Simple Attenuation Form

For a first global model, use an exponential attenuation form:

```text
H_out = H_in * exp(-k * W)

Where:

H_in = incoming wave height or exposure proxy
H_out = wave height or exposure proxy at the shoreline
W = effective habitat width
k = habitat-specific attenuation coefficient
```

This can be turned into a relative protection score:

```text
relative_reduction = 1 - exp(-k * W)
```

For scenario runs:

```text
risk_current = hazard * exposure * (1 - attenuation_current)
risk_removed = hazard * exposure
risk_restored = hazard * exposure * (1 - attenuation_restored)
```

Then:

```text
risk_increase_if_removed = risk_removed - risk_current
risk_reduction_if_restored = risk_current - risk_restored
```

This is not a full hydrodynamic model. It is a scalable screening model that can be run many times globally.

### Scenario Inputs

Scenario design should be simple and editable.

Examples:

```yaml
scenarios:
  current:
    use_observed_habitat: true

  removed:
    habitat_width_multiplier: 0

  restored:
    mangrove_width_added_m: 100
    saltmarsh_width_added_m: 100
    coral_reef_condition_multiplier: 1.2
    seagrass_width_added_m: 100
    dune_presence_added: true
```
