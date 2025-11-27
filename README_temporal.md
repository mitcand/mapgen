# Temporal Dwell Time Analysis

A companion tool to the Hex Bin Dwell Time Analysis for visualizing **when** camera attention occurred during high-altitude aircraft scientific missions. Generates faceted heatmaps showing dwell time patterns by day of week and time of day, with each country displayed in its local timezone including automatic daylight saving time handling.

## Table of Contents

- [Overview](#overview)
- [Quick Start](#quick-start)
- [Installation](#installation)
- [Configuration](#configuration)
- [Input Data Format](#input-data-format)
- [Usage Examples](#usage-examples)
- [Output](#output)
- [Timezone Handling](#timezone-handling)
- [Customization Guide](#customization-guide)
- [Troubleshooting](#troubleshooting)
- [Technical Notes](#technical-notes)

---

## Overview

This script analyzes the same 15 Hz camera metadata logs as the hexbin analysis, but focuses on **temporal patterns** rather than geographic patterns. It answers questions like:

- What days of the week see the most collection activity?
- What times of day (in local time) are most active?
- How do patterns differ between countries?

### Output Heatmaps

The script generates four faceted heatmaps:

| Plot | Data Source | X-Axis | Description |
|------|-------------|--------|-------------|
| `dow_vuefast` | Aircraft position | Day of week | Where the camera was pointing |
| `dow_clipmarks` | Crew markers | Day of week | What crew marked as interesting |
| `tod_vuefast` | Aircraft position | Time of day | 30-minute slots in LOCAL time |
| `tod_clipmarks` | Crew markers | Time of day | 30-minute slots in LOCAL time |

### Key Features

- **Automatic timezone lookup**: Uses `lutz` package to determine timezone from coordinates
- **Daylight saving time support**: Properly handles DST transitions via `lubridate`
- **Accurate country detection**: Uses Natural Earth boundaries via spatial join
- **Grid-based optimization**: Processes millions of points efficiently
- **Shared database**: Uses same DuckDB backend as hexbin analysis
- **Consistent styling**: Dark theme matching the hexbin heatmaps

---

## Quick Start

```r
# 1. Load the script
source("temporal_analysis.R")

# 2. Ensure data exists (use hexbin script's data generator if needed)
source("hexbin_dwell_analysis.R")
generate_synthetic_data(n_missions = 5)

# 3. Run analysis
result <- run_temporal_analysis(config, col_map)

# 4. View the plots
print(result$plots$dow_vuefast)    # Day of week
print(result$plots$tod_vuefast)    # Time of day

# 5. Check output files
result$output_files
```

---

## Installation

### Required R Packages

The script will attempt to install missing packages automatically, but you can install them manually:

```r
install.packages(c(
  "DBI", "duckdb",           # Database
  "sf",                       # Spatial operations
  "ggplot2",                  # Plotting
  "rnaturalearth",            # Country boundaries
  "rnaturalearthdata",
  "dplyr", "tidyr",           # Data manipulation
  "lubridate",                # Date/time handling with DST support
  "lutz",                     # Timezone lookup from coordinates
  "glue", "cli"               # Utilities
))

# High-resolution country boundaries (from r-universe)
install.packages("rnaturalearthhires", 
                 repos = "https://ropensci.r-universe.dev")
```

### System Requirements

- R version 4.0 or higher
- ~4GB RAM for datasets with 10-20 million points
- Sufficient disk space for DuckDB database

---

## Configuration

All configuration is at the top of the script in the `config` list.

### Data Source Settings

```r
config <- list(
  data_dir = "data/missions",    # Folder containing your CSV files
  db_path = "missions.duckdb",   # Database file (shared with hexbin script)
  use_existing_db = FALSE,       # Set TRUE after first run to skip import
  ...
)
```

**Tip**: If you've already run the hexbin analysis, set `use_existing_db = TRUE` to reuse the imported data.

### Query Filters

Filter data the same way as the hexbin script:

**Option A: By Mission ID**
```r
config$mission_ids <- c("M2024001", "M2024002")  # Specific missions
config$mission_ids <- NULL                        # All missions
```

**Option B: By Date Range**
```r
config$mission_ids <- NULL
config$date_start <- "2024-01-01"
config$date_end <- "2024-03-31"
```

### Timezone Lookup Method

```r
config$tz_lookup_method <- "fast"      # Fast Rcpp lookup (default)
config$tz_lookup_method <- "accurate"  # Slower sf spatial join
```

| Method | Speed | Accuracy | Best For |
|--------|-------|----------|----------|
| `fast` | ~10x faster | Good except near borders | Most use cases |
| `accurate` | Slower | Precise at borders | Critical timezone needs |

### Display Settings

```r
config$max_countries <- 12      # Maximum countries to show in facets
config$vertical_layout <- FALSE # Set TRUE for single-column stacked layout
```

**Layout options:**
- `vertical_layout = FALSE` (default): Grid layout with multiple columns
- `vertical_layout = TRUE`: Single column with facets stacked vertically. Each row has its own time/day axis labels for easy reading at any point in the image. Country labels appear on the left. Output height adjusts automatically based on the number of countries.

Facets are ordered by total dwell time (highest first), so you'll always see the most active countries.

### Colors

```r
# Heatmap gradient (low to high time)
config$heatmap_colors <- c("#285DAB", "#5EA5DA", "#F1DCAA", "#F9A63F", "#CD5821")

# Dark theme (matches hexbin maps)
config$background_color <- "#1a1a1a"
config$text_color <- "#ffffff"
```

---

## Input Data Format

### Required Files

Uses the same files as the hexbin analysis:

| File Pattern | Description | Required |
|--------------|-------------|----------|
| `<mission_id>_vuefast.csv` | 15 Hz sensor data | **Yes** |
| `<mission_id>_clipmarks.csv` | Crew-marked POIs | No (but recommended) |

### Column Mapping

Update `col_map` to match your CSV column names:

```r
col_map <- list(
  # VUEFAST columns
  timestamp = "your_timestamp_column",
  frame_lat = "your_frame_lat_column",
  frame_lon = "your_frame_lon_column",
  
  # CLIPMARKS columns
  clip_start_time = "your_start_time_column",
  clip_end_time = "your_end_time_column",
  clip_lat = "your_poi_lat_column",
  clip_lon = "your_poi_lon_column"
)
```

### Timestamp Format

Timestamps should be in a format R can parse (ISO 8601 recommended):
- `2024-03-15 14:30:00`
- `2024-03-15T14:30:00Z`

The script assumes timestamps are in **UTC**. They are converted to local time for each point based on its coordinates.

---

## Usage Examples

### Analyze All Missions

```r
config$mission_ids <- NULL
config$date_start <- NULL
config$date_end <- NULL
result <- run_temporal_analysis(config, col_map)
```

### Analyze Specific Missions

```r
config$mission_ids <- c("M2024015", "M2024016")
result <- run_temporal_analysis(config, col_map)
```

### Analyze a Date Range

```r
config$mission_ids <- NULL
config$date_start <- "2024-06-01"
config$date_end <- "2024-06-30"
result <- run_temporal_analysis(config, col_map)
```

### Quick Re-run (Skip Import)

```r
config$use_existing_db <- TRUE
result <- run_temporal_analysis(config, col_map)
```

### Use Accurate Timezone Lookup

```r
config$tz_lookup_method <- "accurate"
result <- run_temporal_analysis(config, col_map)
```

---

## Output

### Generated Files

Four image files are saved to `config$output_dir`:

| Analysis | Filename Pattern |
|----------|------------------|
| Day of Week - Aircraft | `{prefix}_dow_vuefast_{timestamp}.png` |
| Day of Week - Clipmarks | `{prefix}_dow_clipmarks_{timestamp}.png` |
| Time of Day - Aircraft | `{prefix}_tod_vuefast_{timestamp}.png` |
| Time of Day - Clipmarks | `{prefix}_tod_clipmarks_{timestamp}.png` |

Where `{prefix}` is:
- Mission IDs: `M2024001_M2024002`
- Date range: `90days`
- All missions: `all_missions`

### Accessing Results

```r
result <- run_temporal_analysis(config, col_map)

# Individual plots (ggplot objects)
result$plots$dow_vuefast
result$plots$tod_vuefast
result$plots$dow_clipmarks
result$plots$tod_clipmarks

# File paths
result$output_files$dow_vuefast
# etc.

# Display a plot
print(result$plots$tod_vuefast)
```

### Reading the Plots

**Day of Week Heatmaps:**
- Facet labels show full country name and total time (e.g., "Iran (73.1 hrs)")
- X-axis: Monday through Sunday (in local time)
- Each cell shows minutes of dwell time
- Values printed in cells for easy reading

**Time of Day Heatmaps:**
- Facet labels show full country name, total time, and UTC offset
- X-axis: 00:00 to 24:00 in 30-minute slots
- Times are LOCAL to each country with DST applied

---

## Timezone Handling

### How It Works

1. **Timezone Lookup**: The `lutz` package determines the IANA timezone (e.g., "Asia/Tehran") from each point's coordinates

2. **Grid Optimization**: To avoid looking up 20M+ points individually:
   - Coordinates are rounded to a 0.1° grid (~10km)
   - Timezone is looked up once per unique grid cell
   - Results are applied to all points in that cell

3. **Automatic DST Handling**: `lubridate::with_tz()` converts each UTC timestamp to local time, automatically applying the correct offset based on that specific timestamp's date. Data spanning an entire year will have DST transitions handled correctly - no manual configuration needed.

### Example

A point at coordinates in Iran on March 20, 2024:
- `lutz` returns: `"Asia/Tehran"`
- UTC timestamp: `2024-03-20 10:00:00 UTC`
- Iran observes DST starting March 21, so on March 20:
  - Standard time offset: UTC+3:30
  - Local time: `2024-03-20 13:30:00` (1:30 PM)

The same point on March 22 (after DST starts):
- UTC timestamp: `2024-03-22 10:00:00 UTC`
- DST offset: UTC+4:30
- Local time: `2024-03-22 14:30:00` (2:30 PM)

**Note**: The UTC offset shown in facet labels reflects the *current* offset (at time of running the script), but the actual time conversion for each data point uses the historically correct offset for that timestamp's date.

### Timezone Display

Each facet label shows the full country name and UTC offset:
```
Iran (73.1 hrs)
UTC+3:30
```

Fractional offsets are displayed with minutes (e.g., "UTC+3:30" for Iran, "UTC+5:45" for Nepal).

### Countries with Multiple Timezones

For countries spanning multiple timezones (Russia, USA, Kazakhstan, etc.), the script:
1. Looks up the actual timezone for each point's coordinates
2. Displays the most common timezone in that country's facet label
3. Still uses the correct local time for each individual point

---

## Customization Guide

### Changing Colors

```r
# Blue to red gradient
config$heatmap_colors <- c("#313695", "#4575b4", "#ffffbf", "#f46d43", "#a50026")

# Viridis-style
config$heatmap_colors <- c("#440154", "#3b528b", "#21918c", "#5ec962", "#fde725")
```

### Adjusting Number of Countries

```r
config$max_countries <- 8    # Show fewer countries (larger facets)
config$max_countries <- 16   # Show more countries (smaller facets)
```

### Changing Output Format

```r
config$output_format <- "pdf"   # Vector format for reports
config$output_format <- "png"   # Raster for presentations
config$output_format <- "jpg"   # Smaller file size
```

### Adjusting Plot Dimensions

```r
config$output_width <- 16    # Wider for more countries
config$output_height <- 12   # Taller for more detail
config$output_dpi <- 150     # Lower DPI for faster saves
```

---

## Troubleshooting

### "No vuefast CSV files found"

- Check that `config$data_dir` points to the correct folder
- Verify files are named `<something>_vuefast.csv`
- Run `generate_synthetic_data()` from hexbin script to create test data

### "Column not found" errors

- Update `col_map` to match your actual CSV column names
- Check for typos and case sensitivity

### Script hangs during "Assigning countries and timezones"

This step does spatial joins, which can take a few minutes for 20M+ points. The grid-based optimization helps, but very large datasets may still take time.

If it takes more than 5-10 minutes:
- Check available RAM
- Filter to fewer missions
- Consider using `tz_lookup_method = "fast"`

### Unexpected local times

- Verify your input timestamps are in UTC
- Check that coordinates are valid (lat: -90 to 90, lon: -180 to 180)
- Try `tz_lookup_method = "accurate"` for border regions

### "Error in with_tz" or timezone errors

- Ensure `lubridate` package is up to date
- Some obscure timezone names may not be recognized by R
- Check R's timezone database: `OlsonNames()`

### Database errors

- Delete `missions.duckdb` file and run again
- Or set `config$use_existing_db <- FALSE`

---

## Technical Notes

### Country Assignment

Uses Natural Earth boundaries via `sf` spatial join:

1. Round all coordinates to 0.1° grid (~10km cells)
2. Get unique grid cells (typically 10-50K instead of millions)
3. Spatial join grid cells to country polygons
4. Lookup timezone for each grid cell via `lutz`
5. Apply country and timezone to all points

### Timezone Lookup with lutz

The `lutz` package provides two methods:

| Method | Implementation | Pros | Cons |
|--------|----------------|------|------|
| `fast` | Rcpp + precomputed grid | Very fast, works in oceans | Less accurate near borders |
| `accurate` | sf spatial join | Precise boundaries | Slower, requires sf |

### DST Handling

`lubridate::with_tz()` handles DST correctly:
- Uses the system's timezone database (Olson/IANA)
- Applies DST rules based on the actual date of each timestamp
- Handles historical DST changes

### Data Flow

```
CSV Files → DuckDB → Query → Coordinate Validation → 
  → Grid-based Country + Timezone Assignment → 
  → DST-aware Local Time Conversion → 
  → Aggregate by Country + Time Dimension → Plot → Save
```

### Performance

| Dataset Size | Timezone Lookup | Total Runtime |
|--------------|-----------------|---------------|
| 1M points | ~10 seconds | ~30 seconds |
| 10M points | ~30 seconds | ~2 minutes |
| 20M points | ~60 seconds | ~3-4 minutes |

### Memory Usage

Approximate memory needs (similar to hexbin script):
- 1 million points: ~500 MB
- 10 million points: ~2-3 GB
- 20 million points: ~4-5 GB

---

## Relationship to Hexbin Analysis

This script is designed to work alongside `hexbin_dwell_analysis.R`:

| Aspect | Hexbin | Temporal |
|--------|--------|----------|
| Focus | Where (geography) | When (time patterns) |
| Output | Map with hex bins | Faceted heatmaps |
| Database | Creates DuckDB | Can reuse same DB |
| Countries | Optional filter | Automatic faceting |
| Timezones | N/A | Automatic via lutz |

**Workflow tip**: Run hexbin first to identify countries of interest, then run temporal to understand when activity occurred in those countries.

---

## Support

For issues or questions:
1. Check [Troubleshooting](#troubleshooting) section
2. Review console output for error messages
3. Verify input data format matches expected structure

---

*Generated with Claude AI assistance - Version 1.2*
