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
- **Crew interest**: Areas where the crew marked points of interest (highlighted with borders)

The hex bins are colored by dwell time category, allowing you to quickly see which areas received the most attention during a mission or time period.

### Key Features

- **Fast processing**: Uses H3 hexagonal indexing (~10 seconds for 16M points)
- **Flexible queries**: Filter by mission ID(s) or date range
- **Smart projections**: Auto-selects appropriate map projection based on data extent
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
  "sf", "s2",                 # Spatial operations
  "h3r",                      # H3 hexagonal indexing
  "ggplot2", "viridis",       # Plotting
  "rnaturalearth",            # Basemaps
  "rnaturalearthdata",
  "dplyr", "lubridate",       # Data manipulation
  "glue", "cli"               # Utilities
))

# High-resolution country boundaries (from r-universe)
install.packages("rnaturalearthhires", 
                 repos = "https://ropensci.r-universe.dev")
```

### System Requirements

- R version 4.0 or higher
- ~4GB RAM for datasets with 10-20 million points
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

| Diameter | Best For | Detail Level |
|----------|----------|--------------|
| 8 km | Single mission analysis | High detail |
| 25 km | Regional/multi-mission | Medium detail |
| 60 km | Large area coverage | Low detail, fast |

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

### Colors

```r
# Dwell time colors (low to high) - must match number of categories
config$dwell_colors <- c("#285DAB", "#5EA5DA", "#F1DCAA", "#F9A63F", "#CD5821")

# Basemap
config$water_color <- "#1a3a5c"   # Ocean
config$land_color <- "#343434"    # Land

# Clipmark highlighting
config$clipmark_border_color <- "cyan"
config$clipmark_border_width <- 0.8
```

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

### Adjusting Legend Labels

Find the `break_labels` code in the `bin_dwell_times_fast()` function:

```r
break_labels <- c(
  paste0("<", time_breaks[1], " min"),
  paste0(time_breaks[-length(time_breaks)], "-", time_breaks[-1], " min"),
  paste0(">", time_breaks[length(time_breaks)], " min")
)
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

### Map looks distorted

- Try a different projection: `config$projection_override <- "winkel3"`
- For very wide areas, use `"robinson"` or `"mercator"`

### Colors don't match categories

- Ensure `length(config$dwell_colors) == length(config$time_breaks_minutes) + 1`

### Database errors on re-run

- Delete `missions.duckdb` file and run again
- Or set `config$use_existing_db <- FALSE`

---

## Technical Notes

### Why H3?

H3 (developed by Uber) is a hierarchical hexagonal grid system. Key advantages:

1. **Speed**: Assigns points to hexes using pure math (no spatial predicates)
2. **Equal area**: All hexes at a given resolution have the same area
3. **Consistent neighbors**: Each hex always has exactly 6 neighbors
4. **Hierarchical**: Can easily aggregate to coarser resolutions

### Data Flow

```
CSV Files → DuckDB → Query → Coordinate Validation → H3 Binning → 
  → Count Points per Hex → Classify Dwell Time → Generate Polygons → 
  → Project to Map CRS → Render with ggplot2 → Save Image
```

### Performance Tips

1. Set `use_existing_db = TRUE` after first import
2. Use larger hex sizes for faster processing
3. Filter by mission or date to reduce data volume
4. Close other applications to free RAM for large datasets

### Memory Usage

Approximate memory needs:
- 1 million points: ~500 MB
- 10 million points: ~2-3 GB
- 20 million points: ~4-5 GB

---

## Support

For issues or questions:
1. Check [Troubleshooting](#troubleshooting) section
2. Review console output for error messages
3. Verify input data format matches expected structure

---

*Generated with Claude AI assistance - Version 1.0*
