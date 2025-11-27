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
- **Dwell time**: How long the camera looked at each area (measured in hex bins)
- **Crew interest**: Areas where the crew marked points of interest (highlighted with cyan borders)

The hex bins are colored by dwell time category, allowing you to quickly see which areas received the most attention during a mission or time period.

### Key Features

- **Fast processing**: Uses H3 hexagonal indexing (~20 seconds for 35M points)
- **Flexible queries**: Filter by mission ID(s) or date range
- **Smart projections**: Auto-selects appropriate map projection based on data extent
- **Aspect ratio filling**: Map automatically expands to fill output dimensions (default 16:9)
- **Outlier handling**: Uses trimmed extent (1st-99th percentile) to prevent transit legs from stretching the map
- **Configurable styling**: Colors, hex size, time buckets all adjustable
- **DuckDB backend**: Efficient handling of large CSV datasets (20-30GB+)

---

## Quick Start

```r
# 1. Load the script
source("hexbin_dwell_analysis.R")

# 2. Generate test data (if you don't have real data yet)
generate_synthetic_data(n_missions = 5)

# 3. Run analysis
result <- run_analysis(config, col_map)

# 4. View the plot
print(result$plot)

# 5. Check output file location
result$output_file
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

# Clipmark highlighting
config$clipmark_border_color <- "cyan"
config$clipmark_border_width <- 0.8
```

### Output Dimensions

```r
config$output_width <- 16    # Width in inches
config$output_height <- 9    # Height in inches (default 16:9)
config$output_dpi <- 300     # Resolution
```

**The map automatically expands to fill the aspect ratio.** If your data covers a tall region but you want 16:9 output, the map will show more area to the east/west to fill the frame.

### Geographic Focus

```r
# Focus on a specific country (use ISO3 codes)
config$country_filter <- "AZE"   # Azerbaijan
config$country_filter <- "TUR"   # Turkey
config$country_filter <- NULL    # Auto (based on data extent)
```

Find ISO3 codes: https://en.wikipedia.org/wiki/ISO_3166-1_alpha-3

### Projection Override

```r
config$projection_override <- NULL      # Auto-select (recommended)
config$projection_override <- "laea"    # Lambert Azimuthal Equal Area
config$projection_override <- "winkel3" # Winkel Tripel
config$projection_override <- "mercator"
config$projection_override <- "robinson"
config$projection_override <- "albers"
```

**Auto-selection logic:**
- Data touching high latitudes (>65°) or spanning >50° latitude → LAEA
- Small regional (<30° span) → LAEA
- Medium regional (30-60°) → Albers Equal Area Conic
- Continental (60-120°) → Winkel Tripel
- Global (>120°) → Robinson

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
  frame_lat = "your_frame_center_lat_column",
  frame_lon = "your_frame_center_lon_column",
  
  # CLIPMARKS columns (if using)
  clip_start_time = "your_start_time_column",
  clip_lat = "your_poi_lat_column",
  clip_lon = "your_poi_lon_column",
  clip_description = "your_description_column"
)
```

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

### File Naming

Output files are automatically named based on query parameters:

| Query Type | Example Filename |
|------------|------------------|
| Specific missions | `M2024001_M2024002_20251126_143022.png` |
| Date range | `20251126_143022_90days.png` |
| All missions | `all_missions_20251126_143022.png` |
| With country filter | `..._AZE.png` (appended) |

### Accessing Results

```r
result <- run_analysis(config, col_map)

# The ggplot object (for further customization)
result$plot

# Path to saved file
result$output_file

# Display the plot
print(result$plot)
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

- The projection auto-selection should handle most cases
- For high-latitude data, LAEA is automatically selected
- Try `config$projection_override <- "laea"` for problematic regions

### Transit legs stretch the map too much

- The script uses 1st-99th percentile of coordinates by default
- This excludes outlier points (like long transit legs) from the extent calculation
- The data is still plotted, just the map extent is based on the core coverage area

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

The map extent is calculated in three steps:

1. **Trimmed extent**: Uses 1st-99th percentile of hex centroids to exclude outliers
2. **Buffer**: Adds `extent_buffer_km` plus 15% context buffer
3. **Aspect ratio expansion**: Expands in X or Y direction to match output dimensions

This ensures transit legs don't stretch the map while still showing geographic context.

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

*Generated with Claude AI assistance - Version 1.1*
