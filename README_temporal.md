# Temporal Dwell Time Analysis

Analyzes **when** camera attention occurred during high-altitude aircraft scientific missions. Generates faceted heatmaps showing dwell time patterns by day of week and time of day, with automatic local timezone conversion including daylight saving time handling.

## Table of Contents

- [Overview](#overview)
- [Quick Start](#quick-start)
- [Installation](#installation)
- [Configuration Reference](#configuration-reference)
- [Column Mapping](#column-mapping)
- [Input Data Format](#input-data-format)
- [Usage Examples](#usage-examples)
- [Output](#output)
- [Timezone Handling](#timezone-handling)
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

- **Timestamp-based dwell time** - Uses actual timestamp differences, not sample counting
- **Discrete color buckets** - Matches hexbin_dwell_analysis.R style for visual consistency
- **Automatic timezone lookup** from coordinates using `lutz` package
- **Automatic DST handling** - correctly handles data spanning entire years
- **Full country names** in facet labels with UTC offset
- **Two layout options** - grid or vertical single-column
- **Accurate country detection** using Natural Earth boundaries

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

The script will attempt to install missing packages automatically:

```r
install.packages(c(
  "DBI", "duckdb",           # Database
  "sf",                       # Spatial operations
  "ggplot2", "scales",        # Plotting
  "rnaturalearth",            # Country boundaries
  "rnaturalearthdata",
  "dplyr", "tidyr",           # Data manipulation
  "lubridate",                # Date/time handling
  "lutz",                     # Timezone lookup from coordinates
  "patchwork", "cowplot",     # Plot composition
  "glue", "cli"               # Utilities
))

# High-resolution country boundaries (from r-universe)
install.packages("rnaturalearthhires", 
                 repos = "https://ropensci.r-universe.dev")
```

### System Requirements

- R version 4.0 or higher
- ~4GB RAM for datasets with 10-20 million points

---

## Configuration Reference

All configuration is at the top of the script in the `config` list.

### Data Source

| Parameter | Default | Description |
|-----------|---------|-------------|
| `data_dir` | `"data/missions"` | Folder containing CSV files |
| `db_path` | `"missions.duckdb"` | DuckDB database path (created automatically) |
| `use_existing_db` | `FALSE` | `TRUE` = reuse existing DB, `FALSE` = reimport CSVs |

**Tip**: Set `use_existing_db = TRUE` after the first run to skip CSV import.

### Query Filters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `mission_ids` | `NULL` | Vector of mission IDs, e.g., `c("M2024001", "M2024002")`. `NULL` = all missions |
| `date_start` | `NULL` | Start date filter, e.g., `"2024-01-01"` |
| `date_end` | `NULL` | End date filter, e.g., `"2024-03-31"` |

### Dwell Time Buckets

| Parameter | Default | Description |
|-----------|---------|-------------|
| `time_breaks_minutes` | `c(0.5, 2, 5, 10)` | Break points for categorizing dwell time. Creates N+1 buckets. |

**Example**: `c(0.5, 2, 5, 10)` creates these buckets:
- `<0.5` (under 30 seconds)
- `0.5-2` (30 sec to 2 min)
- `2-5` (2 to 5 min)
- `5-10` (5 to 10 min)
- `>10` (over 10 min)

**Important**: Match these to your `hexbin_dwell_analysis.R` config for visual consistency between the geographic and temporal analyses.

### Timezone Lookup

| Parameter | Default | Description |
|-----------|---------|-------------|
| `tz_lookup_method` | `"fast"` | Method for timezone lookup: `"fast"` (Rcpp, may be inaccurate near borders) or `"accurate"` (sf spatial join, slower but precise) |

### Display

| Parameter | Default | Description |
|-----------|---------|-------------|
| `max_countries` | `12` | Maximum countries to display in faceted plots |
| `vertical_layout` | `FALSE` | `TRUE` = single column with axis labels on each row; `FALSE` = grid layout |
| `facet_cols_dow` | `4` | Number of columns for day-of-week grid layout |
| `facet_cols_tod` | `3` | Number of columns for time-of-day grid layout |

### Colors

| Parameter | Default | Description |
|-----------|---------|-------------|
| `heatmap_colors` | `c("#285DAB", "#5EA5DA", "#F1DCAA", "#F9A63F", "#CD5821")` | Colors for each dwell time bucket (low to high). **Must have exactly `length(time_breaks_minutes) + 1` colors.** |
| `background_color` | `"#1a1a1a"` | Plot background color |
| `text_color` | `"#ffffff"` | Text and label color |
| `grid_color` | `"#333333"` | Grid lines between tiles |

**Note**: The default 5 colors match the default 4 time breaks (creating 5 buckets). If you change `time_breaks_minutes`, update `heatmap_colors` to have the matching count.

### Output

| Parameter | Default | Description |
|-----------|---------|-------------|
| `output_dir` | `"output"` | Directory for saved images |
| `output_format` | `"png"` | Format: `"png"`, `"pdf"`, or `"jpg"` |
| `output_width` | `14` | Image width in inches (grid layout) |
| `output_height` | `10` | Image height in inches (grid layout) |
| `output_dpi` | `300` | Resolution (300 = print quality) |

**Note**: Vertical layout automatically calculates dimensions based on number of countries.

---

## Column Mapping

Update `col_map` to match your CSV column names:

```r
col_map <- list(
  # VUEFAST columns (15 Hz sensor data)
  timestamp = "timestamp",           # Timestamp column
  frame_lat = "frame_center_lat",    # Camera frame center latitude
  frame_lon = "frame_center_lon",    # Camera frame center longitude
  
  # CLIPMARKS columns (crew-marked points of interest)
  clip_start_time = "start_time",    # Clip start timestamp
  clip_end_time = "end_time",        # Clip end timestamp
  clip_lat = "poi_lat",              # Point of interest latitude
  clip_lon = "poi_lon"               # Point of interest longitude
)
```

### Finding Your Column Names

1. Open one of your vuefast CSV files
2. Look at the header row
3. Replace the placeholder names with your actual column names

---

## Input Data Format

### Required Files

| File Pattern | Description | Required |
|--------------|-------------|----------|
| `<mission_id>_vuefast.csv` | 15 Hz sensor data | **Yes** |
| `<mission_id>_clipmarks.csv` | Crew-marked POIs | No (but recommended) |

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

### Use Vertical Layout

```r
config$vertical_layout <- TRUE
result <- run_temporal_analysis(config, col_map)
```

### Match Hexbin Color Buckets

```r
# Use same breaks as hexbin_dwell_analysis.R
config$time_breaks_minutes <- c(0.5, 2, 5, 10)
config$heatmap_colors <- c("#285DAB", "#5EA5DA", "#F1DCAA", "#F9A63F", "#CD5821")
result <- run_temporal_analysis(config, col_map)
```

### Custom Bucket Sizes

```r
# Finer granularity for shorter dwell times
config$time_breaks_minutes <- c(0.25, 0.5, 1, 2, 5)
config$heatmap_colors <- c("#440154", "#3b528b", "#21918c", "#5ec962", "#fde725", "#ff0000")
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
```

### Reading the Plots

**Day of Week Heatmaps:**
- Facet labels show full country name and total time
- X-axis: Monday through Sunday
- Each cell shows minutes of dwell time
- Colors indicate dwell time bucket (legend at bottom)

**Time of Day Heatmaps:**
- Facet labels show full country name, total time, and UTC offset
- X-axis: 00:00 to 24:00 in 30-minute slots
- Times are LOCAL to each country
- Colors indicate dwell time bucket (legend at bottom)

**Vertical Layout:**
- Each row has its own time/day axis labels for easy reading
- Country labels on the left side
- Useful for tall displays or printing

---

## Timezone Handling

### How It Works

1. **Timezone Lookup**: The `lutz` package determines the IANA timezone (e.g., "Asia/Tehran") from each point's coordinates

2. **Grid Optimization**: Coordinates are rounded to a 0.1° grid (~10km). Timezone is looked up once per unique grid cell, then applied to all points in that cell.

3. **Automatic DST**: `lubridate::with_tz()` converts each UTC timestamp to local time, automatically applying the correct offset based on that specific timestamp's date. Data spanning an entire year will have DST transitions handled correctly.

### Timezone Display

Each time-of-day facet label shows the full country name and UTC offset:
```
Iran (73.1 hrs)
UTC+3:30
```

**Note**: The displayed UTC offset reflects the current offset (when the script runs), but the actual time conversion for each data point uses the historically correct offset for that timestamp's date.

### Countries with Multiple Timezones

For countries spanning multiple timezones (Russia, Kazakhstan, etc.):
1. Each point gets the correct timezone for its coordinates
2. The facet label shows the most common timezone in that country's data
3. Individual points still use their correct local time

### Disputed Areas

Some regions return multiple overlapping timezones (e.g., "Asia/Hebron; Asia/Jerusalem"). The script automatically uses the first timezone listed.

---

## Troubleshooting

### "No vuefast CSV files found"

- Check that `config$data_dir` points to the correct folder
- Verify files are named `<something>_vuefast.csv`

### "Column not found" errors

- Update `col_map` to match your actual CSV column names
- Check for typos and case sensitivity

### Colors don't match bucket count

If you see a warning about color count mismatch:
- Count your buckets: `length(time_breaks_minutes) + 1`
- Ensure `heatmap_colors` has exactly that many colors

### Script runs slowly

- Use `config$tz_lookup_method = "fast"` instead of "accurate"
- Filter to fewer missions or a shorter date range
- Set `config$use_existing_db = TRUE` after first run

### Database errors

- Delete `missions.duckdb` file and run again
- Or set `config$use_existing_db <- FALSE`

---

## Technical Notes

### Dwell Time Calculation

Dwell time is calculated from actual timestamp differences:
1. Data is sorted by timestamp within each country/day or country/time-slot group
2. Time difference between consecutive points is calculated
3. Gaps are capped at 1 minute max (to avoid counting gaps between missions)
4. This approach works regardless of actual data sample rate

### Country Assignment

Uses Natural Earth boundaries (scale 50) via `sf` spatial join:
1. Round coordinates to 0.1° grid (~10km cells)
2. Spatial join grid cells to country polygons
3. Apply country code to all points in each cell

### Timezone Lookup

Uses `lutz` package with two methods:
- `"fast"`: Rcpp-based, works in international waters, ~10x faster
- `"accurate"`: sf spatial join with detailed timezone boundaries

### Performance

| Dataset Size | Timezone Lookup | Total Runtime |
|--------------|-----------------|---------------|
| 1M points | ~10 seconds | ~30 seconds |
| 10M points | ~30 seconds | ~2 minutes |
| 20M points | ~60 seconds | ~3-4 minutes |

### Memory Usage

- 1 million points: ~500 MB
- 10 million points: ~2-3 GB
- 20 million points: ~4-5 GB

---

## Relationship to Hexbin Analysis

This script works alongside `hexbin_dwell_analysis.R`:

| Aspect | Hexbin | Temporal |
|--------|--------|----------|
| Focus | Where (geography) | When (time patterns) |
| Output | Map with hex bins | Faceted heatmaps |
| Color Scale | Discrete buckets | Discrete buckets (same) |
| Database | Creates DuckDB | Can reuse same DB |

**Workflow tip**: 
1. Set matching `time_breaks_minutes` and `heatmap_colors` in both scripts
2. Run hexbin first to identify countries of interest
3. Run temporal to understand when activity occurred
4. Both analyses use consistent color coding for easy comparison

---

*Generated with Claude AI assistance*
