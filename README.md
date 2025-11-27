# Hex Bin Dwell Time Heatmap Analysis

A tool for visualizing camera dwell time patterns from high-altitude aircraft scientific missions. Generates hex bin heatmaps showing where the camera spent the most time, highlighting areas of crew interest.

## Table of Contents

- [Overview](#overview)
- [Quick Start](#quick-start)
- [Installation](#installation)
- [Configuration](#configuration)
- [Input Data Format](#input-data-format)
- [Usage Examples](#usage-examples)
- [Output](#output)
- [Customization Guide](#customization-guide)
- [Troubleshooting](#troubleshooting)
- [Technical Notes](#technical-notes)

---

## Overview

This script analyzes 15 Hz camera metadata logs to create heatmaps showing:
- **Camera frame center dwell time**: Where the camera was looking and for how long
- **Aircraft position dwell time**: Where the aircraft was flying and for how long
- **Crew interest**: Areas where the crew marked points of interest (highlighted with cyan borders)

Each run generates **two maps** by default: one for camera frame centers and one for aircraft position, allowing you to compare where the aircraft flew versus where the camera was pointed. The aircraft map can be disabled via `config$generate_aircraft_map <- FALSE`.

### Key Features

- **Dual heatmaps**: Generates both frame center and aircraft position maps (aircraft map toggleable)
- **Fast processing**: Uses H3 hexagonal indexing (~20 seconds for 35M points)
- **Flexible queries**: Filter by mission ID(s) or date range
- **Smart projections**: Auto-selects appropriate map projection based on data extent
- **Flexible extent control**: Six extent modes including data-based, clipmark-based, country, or manual
- **Overseas territory handling**: Automatically filters to polygon parts containing data (e.g., France without French Guiana)
- **Aspect ratio filling**: Map automatically expands to fill output dimensions (default 16:9)
- **Global map support**: Handles dateline-spanning countries and large extents cleanly
- **Configurable styling**: Colors, hex transparency, borders, time buckets all adjustable
- **DuckDB backend**: Efficient handling of large CSV datasets (20-30GB+)

### Frame Center vs Aircraft Position

| Map Type | What It Shows | Use When... |
|----------|---------------|-------------|
| **Frame Center** | Where the camera was pointed | Analyzing crew interest and target coverage |
| **Aircraft Position** | Where the aircraft flew | Analyzing flight patterns and route coverage |

The two maps often look quite different because the camera is typically pointed off to the side of the aircraft, not straight down. Comparing them helps identify areas where the aircraft made multiple passes versus areas that received sustained camera attention.

---

## Quick Start

```r
# 1. Load the script
source("hexbin_dwell_analysis.R")

# 2. Generate test data (choose one):
generate_synthetic_data(n_missions = 5)           # Random routes
# OR
generate_from_kml(kml_dir = "example_routes")     # From KML files

# 3. Run analysis (generates TWO maps)
result <- run_analysis(config, col_map)

# 4. View the plots
print(result$plots$frame_center)  # Where the camera was looking
print(result$plots$aircraft)       # Where the aircraft was flying

# 5. Check output file locations
result$output_files
# $frame_center: "output/all_missions_..._framecenter.png"
# $aircraft:     "output/all_missions_..._aircraft.png"
```

---

## Test Data Generation

### Option A: Random Missions (Quick Testing)

```r
generate_synthetic_data(output_dir = "data/missions", n_missions = 5)
```

Generates random flight paths with orbits and transits near a base location.

### Option B: KML Route Files (Controlled Testing)

For testing specific geographic scenarios (projection edge cases, dateline crossings, etc.), you can draw routes in Google Earth and use them to generate test data:

```r
# 1. Create .kml files with LineString paths
#    - Open Google Earth Pro
#    - Draw paths using Add > Path
#    - Save each path as a .kml file (not .kmz)

# 2. Place KML files in a folder
#    Default: data/routes/

# 3. Generate missions from routes
generate_from_kml(
  kml_dir = "data/routes",      # Folder with .kml files
  output_dir = "data/missions",  # Output folder
  flight_speed_kts = 300         # Simulated flight speed
)
```

Each `.kml` file becomes one mission. The generator:
- Parses LineString coordinates from the KML
- Interpolates waypoints to 15 Hz sample rate
- Generates camera frame centers (mixed nadir/side-looking)
- Creates synthetic clipmarks at waypoints

### Example KML Routes Included

The `example_routes/` folder contains test routes for common edge cases:

| File | Scenario | Tests |
|------|----------|-------|
| `01_central_asia.kml` | Mid-latitude baseline | Normal operation |
| `02_arctic_high_lat.kml` | Arctic (60-75°N) | High-latitude projection |
| `03_dateline_crossing.kml` | Pacific dateline | Coordinate wrapping |
| `04_long_transit_med_arctic.kml` | Mediterranean to Arctic | Large latitude span (50°+) |
| `05_equatorial_wide.kml` | Equatorial, wide longitude | Wide aspect ratio |
| `06_southern_hemisphere.kml` | South America | Negative latitudes |

To use them:
```r
generate_from_kml(kml_dir = "example_routes", output_dir = "data/missions")
config$use_existing_db <- FALSE  # Force reimport
result <- run_analysis(config, col_map)
```

---

## Installation

### Required R Packages

The script will attempt to install missing packages automatically, but you can install them manually:

```r
install.packages(c(
  "DBI", "duckdb",           # Database

"sf", "s2",                # Spatial operations
  "h3r",                     # H3 hexagonal indexing (via h3lib)
  "ggplot2", "viridis",      # Plotting
  "rnaturalearth",           # Basemaps
  "rnaturalearthdata",
  "dplyr", "lubridate",      # Data manipulation
  "glue", "cli"              # Utilities
))

# High-resolution country boundaries (from r-universe)
install.packages("rnaturalearthhires", 
                 repos = "https://ropensci.r-universe.dev")
```

**Note**: The `h3r` package provides access to Uber's H3 hexagonal indexing system. It requires `h3lib` which is installed automatically as a dependency.

### System Requirements

- R version 4.0 or higher
- ~4GB RAM for datasets with 10-20 million points
- ~8GB RAM for datasets with 30-40 million points
- Sufficient disk space for DuckDB database (~2x CSV file size)

---

## Configuration

All configuration is done in the `config` list at the top of the script.

### Data Source Settings

```r
config <- list(
  data_dir = "data/missions",    # Folder containing your CSV files
  db_path = "missions.duckdb",   # Database file (created automatically)
  use_existing_db = FALSE,       # Set TRUE after first run to skip import
  ...
)
```

**Tip**: After your first successful run, set `use_existing_db = TRUE` to skip the CSV import step on subsequent runs. This significantly speeds up iteration.

### Query Filters

You can filter data two ways:

**Option A: By Mission ID**
```r
config$mission_ids <- c("M2024001", "M2024002")  # Specific missions
config$mission_ids <- NULL                        # All missions
```

**Option B: By Date Range**
```r
config$date_start <- "2024-01-01"
config$date_end <- "2024-03-31"
```

### Hex Grid Size

```r
config$hex_diameter_km <- 25  # Target hex diameter in km
```

The script uses H3 resolution levels. Your target diameter is mapped to the closest H3 resolution:

| Target Diameter | H3 Resolution | Actual Edge Length |
|-----------------|---------------|-------------------|
| ~60 km | 3 | 59.8 km |
| ~25 km | 4 | 22.6 km |
| ~8 km | 5 | 8.5 km |
| ~3 km | 6 | 3.2 km |

### Time Categories

Adjust based on your analysis timeframe:

```r
# Single mission (1-2 days)
config$time_breaks_minutes <- c(0.5, 2, 5, 10)

# 30-90 days
config$time_breaks_minutes <- c(5, 15, 30, 60)

# 90+ days
config$time_breaks_minutes <- c(30, 60, 90, 120)
```

This creates `N+1` categories. For example, `c(0.5, 2, 5, 10)` creates:
- `<0.5 min`, `0.5-2 min`, `2-5 min`, `5-10 min`, `>10 min`

### Colors

```r
# Dwell time colors (low to high) - must have N+1 colors
config$dwell_colors <- c("#285DAB", "#5EA5DA", "#F1DCAA", "#F9A63F", "#CD5821")

# Basemap
config$water_color <- "#1a3a5c"   # Ocean background
config$land_color <- "#343434"    # Land fill

# Clipmark highlighting (hexes with crew-marked targets)
config$clipmark_border_color <- "cyan"
config$clipmark_border_width <- 0.8
config$clipmark_border_alpha <- 1.0    # Transparency (0-1)
```

### Hex Bin Styling

```r
config$hex_alpha <- 0.85           # Fill transparency (0-1)
config$hex_border_color <- "white" # Outline color
config$hex_border_width <- 0.1     # Outline width  
config$hex_border_alpha <- 0.3     # Outline transparency (0-1)
```

### Output Dimensions

```r
config$output_width <- 16    # Width in inches
config$output_height <- 9    # Height in inches (default 16:9)
config$output_dpi <- 300     # Resolution
```

**The map automatically expands to fill the aspect ratio.** If your data covers a tall region but you want 16:9 output, the map will show more area to the east/west to fill the frame.

### Map Generation

```r
config$generate_aircraft_map <- TRUE   # Generate both frame center + aircraft maps
config$generate_aircraft_map <- FALSE  # Generate only frame center map
```

### Geographic Extent

Control how the map bounds are determined with `extent_mode`:

```r
config$extent_mode <- "data_countries"     # Default: countries containing frame data
config$extent_mode <- "data_bbox"          # Raw bounding box of frame data
config$extent_mode <- "clipmark_countries" # Countries containing clipmarks
config$extent_mode <- "clipmark_bbox"      # Raw bounding box of clipmarks
config$extent_mode <- "country"            # Specific country (requires country_filter)
config$extent_mode <- "manual"             # Manual center + extent
```

**For `extent_mode = "country"`:**
```r
config$extent_mode <- "country"
config$country_filter <- "IRN"   # Iran (ISO3 code)
```

**For `extent_mode = "manual"`:**
```r
config$extent_mode <- "manual"
config$manual_center <- c(55, 80)     # c(lat, lon) - center of map
config$manual_extent_deg <- 60        # Width in degrees
```

**Buffer (all modes):**
```r
config$extent_buffer_km <- 50    # Buffer around extent
```

For countries that span the dateline (Russia, USA, Canada, etc.), the script uses sensible "continental" bounds automatically.

Find ISO3 codes: https://en.wikipedia.org/wiki/ISO_3166-1_alpha-3

### Projection Override

```r
config$projection_override <- NULL      # Auto-select (recommended)
config$projection_override <- "laea"    # Lambert Azimuthal Equal Area
config$projection_override <- "albers"  # Albers Equal Area Conic
config$projection_override <- "eqc"     # Equirectangular (Plate Carrée)
config$projection_override <- "robinson"# Robinson (world maps)
config$projection_override <- "mercator"# Mercator
```

Auto-selection logic:
- Small extent (<30°): LAEA
- Medium extent (30-60°): Albers or LAEA
- Large extent (>60° or >80° longitude): Equirectangular

---

## Input Data Format

### Required Files

Place your CSV files in the `data_dir` folder. Files are identified by naming pattern:

| File Pattern | Description | Required |
|--------------|-------------|----------|
| `<mission_id>_vuefast.csv` | 15 Hz sensor data | **Yes** |
| `<mission_id>_clipmarks.csv` | Crew-marked POIs | No |
| `<mission_id>_vueslow.csv` | 1 Hz housekeeping | No (not used) |

### Column Mapping

**IMPORTANT**: Update the `col_map` section to match your actual CSV column names.

```r
col_map <- list(
  # VUEFAST columns
  timestamp = "your_timestamp_column",
  frame_lat = "your_frame_center_lat_column",   # Camera frame center
  frame_lon = "your_frame_center_lon_column",
  aircraft_lat = "your_aircraft_lat_column",    # Aircraft position
  aircraft_lon = "your_aircraft_lon_column",
  
  # CLIPMARKS columns (if using)
  clip_start_time = "your_start_time_column",
  clip_lat = "your_poi_lat_column",
  clip_lon = "your_poi_lon_column",
  clip_description = "your_description_column"
)
```

**Note**: If aircraft position columns are not available in your data, the script will generate only the frame center heatmap and skip the aircraft position map.

### Expected Data

**vuefast.csv** should contain:
- Timestamp column
- Frame center latitude (decimal degrees, -90 to 90)
- Frame center longitude (decimal degrees, -180 to 180)
- Logged at 15 Hz (15 records per second)

**clipmarks.csv** should contain:
- Start/end timestamps
- POI coordinates
- Description text

---

## Usage Examples

### Analyze All Missions

```r
config$mission_ids <- NULL
config$date_start <- NULL
config$date_end <- NULL
result <- run_analysis(config, col_map)
```

### Analyze Specific Missions

```r
config$mission_ids <- c("M2024015", "M2024016", "M2024017")
result <- run_analysis(config, col_map)
```

### Analyze a Date Range

```r
config$mission_ids <- NULL
config$date_start <- "2024-06-01"
config$date_end <- "2024-06-30"
result <- run_analysis(config, col_map)
```

### Focus on a Country

```r
config$country_filter <- "IRN"  # Iran
config$extent_buffer_km <- 100   # Add buffer around country
result <- run_analysis(config, col_map)
```

### Manual Map Extent (Full Control)

```r
# Center on Central Asia with 80° span
config$manual_center <- c(45, 65)    # 45°N, 65°E
config$manual_extent_deg <- 80       # 80° wide (height auto from 16:9)
config$country_filter <- NULL        # Clear country filter
result <- run_analysis(config, col_map)

# Zoom in on a specific region
config$manual_center <- c(35, 50)    # 35°N, 50°E (Persian Gulf area)
config$manual_extent_deg <- 20       # Tighter zoom
result <- run_analysis(config, col_map)
```

### Quick Re-run (Skip Import)

```r
config$use_existing_db <- TRUE  # Reuse previously imported data
result <- run_analysis(config, col_map)
```

### Custom Aspect Ratio

```r
# Square output
config$output_width <- 10
config$output_height <- 10

# Ultra-wide
config$output_width <- 21
config$output_height <- 9

result <- run_analysis(config, col_map)
```

---

## Output

### Dual Map Output

Each run generates **two heatmaps**:

| Map Type | Filename Suffix | Shows |
|----------|-----------------|-------|
| Frame Center | `_framecenter.png` | Where the camera was looking |
| Aircraft Position | `_aircraft.png` | Where the aircraft was flying |

### File Naming

Output files are automatically named based on query parameters:

| Query Type | Example Filenames |
|------------|-------------------|
| Specific missions | `M2024001_M2024002_20251126_143022_framecenter.png`<br>`M2024001_M2024002_20251126_143022_aircraft.png` |
| Date range | `20251126_143022_90days_framecenter.png`<br>`20251126_143022_90days_aircraft.png` |
| All missions | `all_missions_20251126_143022_framecenter.png`<br>`all_missions_20251126_143022_aircraft.png` |
| With country filter | `..._AZE_framecenter.png` (appended) |

### Accessing Results

```r
result <- run_analysis(config, col_map)

# The ggplot objects (for further customization)
result$plots$frame_center    # Camera frame center heatmap
result$plots$aircraft        # Aircraft position heatmap

# Paths to saved files
result$output_files$frame_center
result$output_files$aircraft

# Display plots
print(result$plots$frame_center)
print(result$plots$aircraft)
```

---

## Customization Guide

### Changing Colors

Edit `config$dwell_colors`. Must have exactly `length(time_breaks_minutes) + 1` colors.

```r
# Example: Blue to Red gradient
config$dwell_colors <- c("#313695", "#4575b4", "#ffffbf", "#f46d43", "#a50026")

# Example: Viridis-style
config$dwell_colors <- c("#440154", "#3b528b", "#21918c", "#5ec962", "#fde725")
```

### Moving the Legend

In `create_heatmap()`, modify `legend.position`:

```r
# Inside plot (current default)
legend.position = c(0.98, 0.5),
legend.justification = c(1, 0.5),

# Outside on right
legend.position = "right",

# Bottom
legend.position = "bottom",

# Hide legend
legend.position = "none",
```

### Adding a Scale Bar or North Arrow

In the `create_heatmap()` function, add before `return(p)`:

```r
# Requires ggspatial package
library(ggspatial)
p <- p + 
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "tr", which_north = "true")
```

### Saving in Different Formats

```r
config$output_format <- "pdf"   # Vector format, good for reports
config$output_format <- "png"   # Raster, good for presentations
config$output_format <- "jpg"   # Smaller file size
```

---

## Troubleshooting

### "No vuefast CSV files found"

- Check that `config$data_dir` points to the correct folder
- Verify files are named `<something>_vuefast.csv`
- Run `generate_synthetic_data()` to create test data

### "Column not found" errors

- Update `col_map` to match your actual CSV column names
- Check for typos and case sensitivity

### Script hangs during binning

- This shouldn't happen with H3 (it's very fast)
- Check available RAM - may need to filter to fewer missions
- Try increasing `config$hex_diameter_km` to reduce hex count

### Map looks distorted or has cut-off edges

- Try `config$manual_center` and `config$manual_extent_deg` for full control
- The projection auto-selection should handle most cases
- For large extents spanning multiple continents, Equirectangular is used
- For high-latitude data, LAEA is automatically selected
- Try `config$projection_override <- "laea"` for problematic regions

### Countries like Russia cause artifacts

- For global views, problematic countries are automatically cropped to the visible extent
- Use `config$extent_mode <- "country"` with `config$country_filter <- "RUS"` to focus on Russia specifically
- Or use `config$extent_mode <- "manual"` with `config$manual_center` and `config$manual_extent_deg` for precise control

### Map zooms out to show distant territories (e.g., French Guiana when data is in France)

- The script automatically breaks country multipolygons into parts and only uses parts containing data
- This should be handled automatically in `data_countries` and `clipmark_countries` modes
- If still problematic, use `extent_mode = "data_bbox"` or `extent_mode = "manual"`

### Transit legs stretch the map too much

- Use `extent_mode = "clipmark_countries"` or `"clipmark_bbox"` to base extent on targets instead of all data
- Or use `extent_mode = "manual"` with specific center/extent values
- The data is still plotted, just the map extent is based on your chosen reference

### Colors don't match categories

- Ensure `length(config$dwell_colors) == length(config$time_breaks_minutes) + 1`

### Database errors on re-run

- Delete `missions.duckdb` file and run again
- Or set `config$use_existing_db <- FALSE`

### Legend overlaps data

- Move legend: `legend.position = "right"` or `legend.position = "bottom"`
- Or make output wider to give more space

---

## Technical Notes

### Why H3?

H3 (developed by Uber) is a hierarchical hexagonal grid system. Key advantages:

1. **Speed**: Assigns points to hexes using pure math (no spatial predicates)
2. **Equal area**: All hexes at a given resolution have the same area
3. **Consistent neighbors**: Each hex always has exactly 6 neighbors
4. **Hierarchical**: Can easily aggregate to coarser resolutions

The `h3r` package provides R bindings via `h3lib`.

### Data Flow

```
CSV Files → DuckDB Import → Query by Mission/Date → 
  → Coordinate Validation → H3 Binning (fast!) → 
  → Count Points per Hex → Classify Dwell Time → 
  → Generate Hex Polygons → Transform to Projection → 
  → Expand to Output Aspect Ratio → Render with ggplot2 → Save
```

### Extent Calculation

The map extent is determined by `extent_mode`:

1. **manual**: Uses `manual_center` and `manual_extent_deg`
2. **country**: Uses `country_filter` ISO3 code (special handling for dateline-spanning countries)
3. **data_countries** (default): Countries containing frame center data
4. **data_bbox**: Raw bounding box of frame center data
5. **clipmark_countries**: Countries containing clipmark locations
6. **clipmark_bbox**: Raw bounding box of clipmarks

After extent is determined:
- Buffer added (`extent_buffer_km`)
- Clamped to valid lat/lon ranges
- Expanded to match output aspect ratio (16:9 default)

**Overseas Territory Handling**

Countries with distant territories (France → French Guiana, UK → Falklands, Netherlands → Caribbean, etc.) are handled by breaking multipolygons into individual parts and only using the parts that actually contain data. This prevents the map from zooming out to show French Guiana when your data is only in metropolitan France.

**Special Country Handling**

Countries that span the dateline or have extreme extents use predefined "continental" bounds:

| Country | Extent Used |
|---------|-------------|
| Russia (RUS) | 20°E to 180°E, 41°N to 82°N |
| USA | 125°W to 66°W, 24°N to 50°N (CONUS) |
| Canada (CAN) | 141°W to 52°W, 42°N to 72°N |
| New Zealand (NZL) | 166°E to 179°E |
| Fiji (FJI) | 177°E to 180°E |

For global views, these countries are cropped to the visible extent to prevent rendering artifacts.

### Performance

Tested performance on Apple Silicon (M-series):

| Points | H3 Indexing | Polygon Generation | Total |
|--------|-------------|-------------------|-------|
| 3.8M | ~2s | ~0.6s | ~5s |
| 16M | ~6s | ~2s | ~15s |
| 35M | ~15s | ~4s | ~25s |

### Memory Usage

Approximate memory needs:
- 1 million points: ~500 MB
- 10 million points: ~2-3 GB
- 20 million points: ~4-5 GB
- 35 million points: ~6-8 GB

### Tips for Large Datasets

1. Set `use_existing_db = TRUE` after first import
2. Use larger hex sizes (`hex_diameter_km = 50`) for faster processing
3. Filter by mission or date to reduce data volume
4. Close other applications to free RAM
5. Consider using date ranges to process in chunks

---

## Support

For issues or questions:
1. Check [Troubleshooting](#troubleshooting) section
2. Review console output for error messages
3. Verify input data format matches expected structure

---

*Generated with Claude AI assistance - Version 1.4*
