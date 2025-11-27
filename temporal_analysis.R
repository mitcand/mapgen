# ==============================================================================
# Temporal Dwell Time Analysis
# ==============================================================================
# 
# PURPOSE:
#   Analyzes camera frame center positions and clipmarks to visualize temporal
#   patterns from high-altitude aircraft scientific missions. Generates four
#   faceted heatmaps showing:
#     1a. Dwell time by day of week (aircraft position from vuefast)
#     1b. Collection time by day of week (clipmarks)
#     2a. Dwell time by time of day in 30-min segments (aircraft position)
#     2b. Collection time by time of day (clipmarks)
#
#   Each facet represents a country, ordered by total time (highest first).
#   Time-of-day plots show each country in its LOCAL timezone with automatic
#   daylight saving time handling.
#
# AUTHOR: Generated with Claude AI assistance
# VERSION: 1.2
# LAST UPDATED: 2025-11-27
#
# USAGE:
#   1. Update the CONFIG section with your parameters
#   2. Update the COLUMN MAPPING section with your CSV column names
#   3. Run: result <- run_temporal_analysis(config, col_map)
#   4. View: print(result$plots$dow_vuefast)
#
# INPUT FILES:
#   - <mission_id>_vuefast.csv  : 15 Hz sensor data with frame center coordinates
#   - <mission_id>_clipmarks.csv: Crew-marked points of interest
#
# OUTPUT:
#   - Four PNG/PDF images saved to output_dir
#   - Returns list with $plots (ggplot objects) and $output_files (paths)
#
# DEPENDENCIES:
#   DBI, duckdb, sf, ggplot2, rnaturalearth, rnaturalearthdata,
#   dplyr, lubridate, glue, cli, tidyr, lutz
#
# ==============================================================================

# ------------------------------------------------------------------------------
# CONFIGURATION - Adjust these parameters as needed
# ------------------------------------------------------------------------------
# 
# This section contains all user-configurable parameters. Modify these values
# to customize the analysis for your specific needs.
#
# TIP: For iterative analysis, set use_existing_db = TRUE after the first run
#      to skip the CSV import step and speed up subsequent runs.
# ------------------------------------------------------------------------------

config <- list(
  
  # === DATA SOURCE ===
  data_dir = "data/missions",           # Folder containing CSV files
  db_path = "missions.duckdb",          # DuckDB database path (created automatically)
  use_existing_db = FALSE,              # TRUE = reuse existing DB, FALSE = reimport CSVs
  
  # === QUERY FILTERS ===
  # Use ONE of these approaches to filter data:
  # Option A: Specific missions
  mission_ids = NULL,                   # e.g., c("M2024001", "M2024002") or NULL for all
  # Option B: Date range
  date_start = NULL,                    # e.g., "2024-01-01" or NULL
  date_end = NULL,                      # e.g., "2024-03-31" or NULL
  
  # === TIMEZONE LOOKUP ===
  # Method for looking up timezones from coordinates (via lutz package)
  # "fast"     - Very fast, uses Rcpp. May be inaccurate near timezone borders.
  # "accurate" - Slower, uses sf spatial join. More accurate near borders.
  tz_lookup_method = "fast",
  
  # === DISPLAY ===
  max_countries = 12,                   # Maximum countries to display in faceted plots
  
  # === COLORS ===
  # Heatmap colors (low to high dwell time) - gradient for time intensity
  heatmap_colors = c("#285DAB", "#5EA5DA", "#F1DCAA", "#F9A63F", "#CD5821"),
  
  # Background styling (dark theme to match hexbin maps)
  background_color = "#1a1a1a",         # Plot background
  text_color = "#ffffff",               # Text and labels
  grid_color = "#333333",               # Grid lines between tiles
  
  # === OUTPUT ===
  output_dir = "output",                # Directory for saved images
  output_format = "png",                # Format: "png", "pdf", or "jpg"
  output_width = 14,                    # Image width in inches
  output_height = 10,                   # Image height in inches
  output_dpi = 300                      # Resolution (300 = print quality)
)

# ------------------------------------------------------------------------------
# COLUMN MAPPING - Update these to match your actual CSV column names
# ------------------------------------------------------------------------------
#
# IMPORTANT: These are placeholder names. You MUST update them to match your
# actual CSV column names before running on real data.
#
# How to find your column names:
#   1. Open one of your vuefast CSV files
#   2. Look at the header row
#   3. Replace the placeholder names below with your actual column names
#
# Example: If your CSV has "gps_timestamp" instead of "timestamp",
#          change: timestamp = "timestamp"
#          to:     timestamp = "gps_timestamp"
# ------------------------------------------------------------------------------

col_map <- list(
  
  # --- VUEFAST COLUMNS (15 Hz sensor data) ---
  timestamp = "timestamp",              # UPDATE: Timestamp column (e.g., "gps_time", "utc")
  frame_lat = "frame_center_lat",       # UPDATE: Camera frame center latitude
  frame_lon = "frame_center_lon",       # UPDATE: Camera frame center longitude
  
  # NOTE: If your timestamp is split across date and time columns, you'll need
  # to modify the query functions to combine them.
  
  # --- CLIPMARKS COLUMNS (crew-marked points of interest) ---
  clip_start_time = "start_time",       # UPDATE: Clip start timestamp
  clip_end_time = "end_time",           # UPDATE: Clip end timestamp
  clip_lat = "poi_lat",                 # UPDATE: Point of interest latitude
  clip_lon = "poi_lon"                  # UPDATE: Point of interest longitude
)

# ==============================================================================
# END OF USER CONFIGURATION
# ==============================================================================
# 
# The code below handles all processing. You should not need to modify it
# unless you want to customize the analysis behavior.
#
# ==============================================================================

# ------------------------------------------------------------------------------
# PACKAGES
# ------------------------------------------------------------------------------

required_packages <- c(
  "DBI", "duckdb",                      # Database
  "sf",                                 # Spatial operations
  "ggplot2",                            # Plotting
  "rnaturalearth", "rnaturalearthdata", # Basemaps/country boundaries
  "dplyr", "tidyr", "lubridate",        # Data manipulation
  "glue", "cli",                        # Utilities
  "lutz"                                # Timezone lookup from coordinates
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# High-resolution boundaries for better country detection
if (!requireNamespace("rnaturalearthhires", quietly = TRUE)) {
  install.packages("rnaturalearthhires", repos = "https://ropensci.r-universe.dev")
}

# ------------------------------------------------------------------------------
# DATABASE FUNCTIONS
# ------------------------------------------------------------------------------

init_temporal_database <- function(config, col_map) {
  cli_alert_info("Initializing DuckDB database: {config$db_path}")
  
  con <- dbConnect(duckdb(), config$db_path)
  
  # Check if we should use existing DB
  if (config$use_existing_db && dbExistsTable(con, "vuefast")) {
    cli_alert_success("Using existing database")
    return(con)
  }
  
  # Find and import CSV files
  csv_files <- list.files(config$data_dir, pattern = "_vuefast\\.csv$", full.names = TRUE)
  
  if (length(csv_files) == 0) {
    cli_alert_warning("No vuefast CSV files found in {config$data_dir}")
    cli_alert_info("Run generate_synthetic_data() from hexbin_dwell_analysis.R to create test data")
    dbDisconnect(con)
    return(NULL)
  }
  
  cli_alert_info("Found {length(csv_files)} mission files")
  
  # Create tables
  dbExecute(con, "DROP TABLE IF EXISTS vuefast")
  dbExecute(con, "DROP TABLE IF EXISTS clipmarks")
  
  # Import vuefast files
  cli_alert("Importing vuefast data...")
  
  for (i in seq_along(csv_files)) {
    file_path <- csv_files[i]
    mission_id <- gsub("_vuefast\\.csv$", "", basename(file_path))
    
    cli_progress_step("Mission {mission_id} ({i}/{length(csv_files)})")
    
    if (i == 1) {
      query <- glue("
        CREATE TABLE vuefast AS 
        SELECT '{mission_id}' as mission_id, * 
        FROM read_csv_auto('{file_path}')
      ")
    } else {
      query <- glue("
        INSERT INTO vuefast 
        SELECT '{mission_id}' as mission_id, * 
        FROM read_csv_auto('{file_path}')
      ")
    }
    dbExecute(con, query)
  }
  
  # Import clipmarks
  clip_files <- list.files(config$data_dir, pattern = "_clipmarks\\.csv$", full.names = TRUE)
  
  if (length(clip_files) > 0) {
    cli_alert("Importing clipmarks data...")
    
    for (i in seq_along(clip_files)) {
      file_path <- clip_files[i]
      mission_id <- gsub("_clipmarks\\.csv$", "", basename(file_path))
      
      if (i == 1) {
        query <- glue("
          CREATE TABLE clipmarks AS 
          SELECT '{mission_id}' as mission_id, * 
          FROM read_csv_auto('{file_path}')
        ")
      } else {
        query <- glue("
          INSERT INTO clipmarks 
          SELECT '{mission_id}' as mission_id, * 
          FROM read_csv_auto('{file_path}')
        ")
      }
      dbExecute(con, query)
    }
  }
  
  # Create indexes
  cli_alert("Creating indexes...")
  dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_vuefast_mission ON vuefast(mission_id)")
  dbExecute(con, glue("CREATE INDEX IF NOT EXISTS idx_vuefast_time ON vuefast({col_map$timestamp})"))
  
  # Report stats
  n_rows <- dbGetQuery(con, "SELECT COUNT(*) as n FROM vuefast")$n
  cli_alert_success("Database ready: {format(n_rows, big.mark=',')} vuefast records")
  
  if (dbExistsTable(con, "clipmarks")) {
    n_clips <- dbGetQuery(con, "SELECT COUNT(*) as n FROM clipmarks")$n
    cli_alert_success("Found {format(n_clips, big.mark=',')} clipmark records")
  }
  
  return(con)
}

# ------------------------------------------------------------------------------
# COUNTRY AND TIMEZONE ASSIGNMENT FUNCTIONS
# ------------------------------------------------------------------------------

assign_countries_and_timezones <- function(coords_df, config, lat_col = "lat", lon_col = "lon") {
  #' Assign ISO3 country codes and IANA timezones to coordinates
  #'
  #' Uses a grid-based lookup for speed:
  #' 1. Round coordinates to 0.1 degree grid (~10km)
  #' 2. Get unique grid cells
  #' 3. Spatial join grid cells to country polygons

  #' 4. Lookup timezones using lutz package
  #' 5. Look up country/timezone for each point from its grid cell
  #'
  #' @param coords_df Data frame with lat/lon columns
  #' @param config Configuration list (for tz_lookup_method)
  #' @param lat_col Name of latitude column
  #' @param lon_col Name of longitude column
  #' @return Data frame with country and tz columns added
  
  n_points <- nrow(coords_df)
  cli_alert("Assigning countries and timezones to {format(n_points, big.mark=',')} coordinates...")
  
  # Disable s2 spherical geometry - Natural Earth polygons sometimes have issues
  sf_use_s2(FALSE)
  
  # Get country boundaries
  countries <- ne_countries(scale = 50, returnclass = "sf")
  
  # Use grid-based lookup for speed
  cli_alert("Creating grid cells...")
  grid_resolution <- 0.1  # ~10km at equator
  
  # Round coordinates to grid
  grid_lat <- round(coords_df[[lat_col]] / grid_resolution) * grid_resolution
  grid_lon <- round(coords_df[[lon_col]] / grid_resolution) * grid_resolution
  grid_id <- paste(grid_lat, grid_lon, sep = "_")
  
  # Get unique grid cells
  unique_grid_ids <- unique(grid_id)
  unique_cells <- data.frame(
    grid_id = unique_grid_ids,
    grid_lat = as.numeric(sapply(strsplit(unique_grid_ids, "_"), `[`, 1)),
    grid_lon = as.numeric(sapply(strsplit(unique_grid_ids, "_"), `[`, 2)),
    stringsAsFactors = FALSE
  )
  
  n_cells <- nrow(unique_cells)
  cli_alert("Processing {format(n_cells, big.mark=',')} unique grid cells...")
  
  # Create sf points for grid cell centers
  cells_sf <- st_as_sf(
    unique_cells,
    coords = c("grid_lon", "grid_lat"),
    crs = 4326,
    remove = FALSE
  )
  
  # Spatial join grid cells to countries
  cli_alert("Assigning countries via spatial join...")
  cells_joined <- st_join(cells_sf, countries[, c("iso_a3")], left = TRUE)
  
  # Lookup timezones using lutz
  cli_alert("Looking up timezones (method: {config$tz_lookup_method})...")
  tz_names <- tz_lookup_coords(
    lat = unique_cells$grid_lat,
    lon = unique_cells$grid_lon,
    method = config$tz_lookup_method
  )
  
  # Create lookup table
  cell_lookup <- data.frame(
    grid_id = cells_joined$grid_id,
    country = cells_joined$iso_a3,
    tz = tz_names,
    stringsAsFactors = FALSE
  )
  cell_lookup$country[cell_lookup$country == "-99"] <- NA
  
  # Look up countries and timezones for all points
  cli_alert("Applying to all points...")
  match_idx <- match(grid_id, cell_lookup$grid_id)
  coords_df$country <- cell_lookup$country[match_idx]
  coords_df$tz <- cell_lookup$tz[match_idx]
  
  # Summary
  n_assigned <- sum(!is.na(coords_df$country))
  n_unassigned <- sum(is.na(coords_df$country))
  cli_alert_success("Assigned {format(n_assigned, big.mark=',')} points to countries ({format(n_unassigned, big.mark=',')} unassigned)")
  
  # Show top countries
  if (n_assigned > 0) {
    country_counts <- sort(table(coords_df$country[!is.na(coords_df$country)]), decreasing = TRUE)
    top_countries <- head(country_counts, 5)
    cli_alert_info("Top countries: {paste(names(top_countries), collapse=', ')}")
  }
  
  # Show timezone info
  tz_counts <- sort(table(coords_df$tz[!is.na(coords_df$tz)]), decreasing = TRUE)
  n_tz <- length(tz_counts)
  cli_alert_info("Found {n_tz} unique timezones")
  
  # Re-enable s2
  sf_use_s2(TRUE)
  
  return(coords_df)
}

# ------------------------------------------------------------------------------
# DATA QUERY FUNCTIONS
# ------------------------------------------------------------------------------

build_where_clause <- function(config, col_map, timestamp_col = NULL) {
  if (is.null(timestamp_col)) {
    timestamp_col <- col_map$timestamp
  }
  
  where_clauses <- c()
  
  if (!is.null(config$mission_ids)) {
    missions_str <- paste0("'", config$mission_ids, "'", collapse = ", ")
    where_clauses <- c(where_clauses, glue("mission_id IN ({missions_str})"))
  }
  
  if (!is.null(config$date_start)) {
    where_clauses <- c(where_clauses, glue("{timestamp_col} >= '{config$date_start}'"))
  }
  
  if (!is.null(config$date_end)) {
    where_clauses <- c(where_clauses, glue("{timestamp_col} <= '{config$date_end}'"))
  }
  
  if (length(where_clauses) > 0) {
    return(paste("WHERE", paste(where_clauses, collapse = " AND ")))
  }
  
  return("")
}

query_vuefast_temporal <- function(con, config, col_map) {
  #' Query vuefast data with temporal components extracted (in UTC)
  #' Country/timezone assignment and local time conversion happens in R
  
  where_sql <- build_where_clause(config, col_map)
  
  query <- glue("
    SELECT 
      {col_map$frame_lat} as lat,
      {col_map$frame_lon} as lon,
      {col_map$timestamp} as timestamp_utc
    FROM vuefast
    {where_sql}
  ")
  
  cli_alert_info("Querying vuefast temporal data...")
  data <- dbGetQuery(con, query)
  cli_alert_success("Retrieved {format(nrow(data), big.mark=',')} points")
  
  # Remove invalid coordinates
  valid <- !is.na(data$lat) & !is.na(data$lon) &
           data$lat >= -90 & data$lat <= 90 &
           data$lon >= -180 & data$lon <= 180
  
  n_invalid <- sum(!valid)
  if (n_invalid > 0) {
    cli_alert_warning("Removed {format(n_invalid, big.mark=',')} rows with invalid coordinates")
    data <- data[valid, ]
  }
  
  # Convert timestamp to POSIXct if needed
  if (!inherits(data$timestamp_utc, "POSIXct")) {
    data$timestamp_utc <- as.POSIXct(data$timestamp_utc, tz = "UTC")
  }
  
  return(data)
}

query_clipmarks_temporal <- function(con, config, col_map) {
  #' Query clipmarks data with temporal components and duration
  
  if (!dbExistsTable(con, "clipmarks")) {
    cli_alert_warning("No clipmarks table found")
    return(NULL)
  }
  
  where_sql <- build_where_clause(config, col_map, col_map$clip_start_time)
  
  query <- glue("
    SELECT 
      {col_map$clip_lat} as lat,
      {col_map$clip_lon} as lon,
      {col_map$clip_start_time} as start_time_utc,
      {col_map$clip_end_time} as end_time_utc
    FROM clipmarks
    {where_sql}
  ")
  
  cli_alert_info("Querying clipmarks temporal data...")
  data <- tryCatch({
    dbGetQuery(con, query)
  }, error = function(e) {
    cli_alert_warning("Error querying clipmarks: {e$message}")
    return(NULL)
  })
  
  if (is.null(data) || nrow(data) == 0) {
    cli_alert_warning("No clipmarks found")
    return(NULL)
  }
  
  cli_alert_success("Retrieved {nrow(data)} clipmarks")
  
  # Remove invalid
  valid <- !is.na(data$lat) & !is.na(data$lon) &
           data$lat >= -90 & data$lat <= 90 &
           data$lon >= -180 & data$lon <= 180
  
  data <- data[valid, ]
  
  # Convert timestamps
  if (!inherits(data$start_time_utc, "POSIXct")) {
    data$start_time_utc <- as.POSIXct(data$start_time_utc, tz = "UTC")
    data$end_time_utc <- as.POSIXct(data$end_time_utc, tz = "UTC")
  }
  
  # Calculate duration in minutes
  data$duration_minutes <- as.numeric(difftime(data$end_time_utc, data$start_time_utc, units = "mins"))
  data <- data[data$duration_minutes > 0, ]
  
  return(data)
}

# ------------------------------------------------------------------------------
# AGGREGATION FUNCTIONS WITH LOCAL TIMEZONE + DST SUPPORT
# ------------------------------------------------------------------------------

add_local_time_components <- function(data, timestamp_col = "timestamp_utc") {
  #' Add day_of_week and time_slot columns based on each point's local timezone
  #' Properly handles daylight saving time via lubridate::with_tz()
  #'
  #' @param data Data frame with tz and timestamp columns
  #' @param timestamp_col Name of the UTC timestamp column
  #' @return Data frame with day_of_week and time_slot columns added
  
  cli_alert("Converting to local timezones (with DST handling)...")
  
  # Get unique timezones
  unique_tz <- unique(data$tz[!is.na(data$tz)])
  cli_alert_info("Processing {length(unique_tz)} unique timezones")
  
  # Initialize columns
  data$local_hour <- NA_real_
  data$local_minute <- NA_real_
  data$day_of_week <- NA_integer_
  
  # Process each timezone separately (with_tz requires single timezone)
  # For efficiency, we vectorize within each timezone group
  for (tz_name in unique_tz) {
    idx <- which(data$tz == tz_name)
    if (length(idx) == 0) next
    
    # Convert UTC to local time for this timezone
    # with_tz properly handles DST based on the actual date
    local_times <- with_tz(data[[timestamp_col]][idx], tzone = tz_name)
    
    # Extract components
    data$local_hour[idx] <- hour(local_times)
    data$local_minute[idx] <- minute(local_times)
    data$day_of_week[idx] <- wday(local_times, week_start = 1)  # 1=Mon, 7=Sun
  }
  
  # Calculate time slot (0-47, each slot = 30 minutes)
  data$time_slot <- floor((data$local_hour * 60 + data$local_minute) / 30)
  
  # Handle any remaining NAs (ocean points, etc.) - use UTC
  na_idx <- which(is.na(data$day_of_week))
  if (length(na_idx) > 0) {
    data$day_of_week[na_idx] <- wday(data[[timestamp_col]][na_idx], week_start = 1)
    data$time_slot[na_idx] <- floor((hour(data[[timestamp_col]][na_idx]) * 60 + 
                                      minute(data[[timestamp_col]][na_idx])) / 30)
  }
  
  cli_alert_success("Local time conversion complete")
  
  return(data)
}

aggregate_by_country_dow <- function(data, sample_rate_hz = 15, is_clipmarks = FALSE) {
  #' Aggregate data by country and day of week
  
  if (is_clipmarks) {
    agg <- data %>%
      filter(!is.na(country)) %>%
      group_by(country, day_of_week) %>%
      summarise(minutes = sum(duration_minutes, na.rm = TRUE), .groups = "drop")
  } else {
    agg <- data %>%
      filter(!is.na(country)) %>%
      group_by(country, day_of_week) %>%
      summarise(n_samples = n(), .groups = "drop") %>%
      mutate(minutes = n_samples / sample_rate_hz / 60) %>%
      select(country, day_of_week, minutes)
  }
  
  return(agg)
}

aggregate_by_country_tod <- function(data, sample_rate_hz = 15, is_clipmarks = FALSE) {
  #' Aggregate data by country and time of day slot
  
  if (is_clipmarks) {
    agg <- data %>%
      filter(!is.na(country)) %>%
      group_by(country, time_slot) %>%
      summarise(minutes = sum(duration_minutes, na.rm = TRUE), .groups = "drop")
  } else {
    agg <- data %>%
      filter(!is.na(country)) %>%
      group_by(country, time_slot) %>%
      summarise(n_samples = n(), .groups = "drop") %>%
      mutate(minutes = n_samples / sample_rate_hz / 60) %>%
      select(country, time_slot, minutes)
  }
  
  return(agg)
}

get_country_totals <- function(data, is_clipmarks = FALSE, sample_rate_hz = 15) {
  #' Get total time per country for ordering facets
  
  if (is_clipmarks) {
    totals <- data %>%
      filter(!is.na(country)) %>%
      group_by(country) %>%
      summarise(total_minutes = sum(duration_minutes, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(total_minutes))
  } else {
    totals <- data %>%
      filter(!is.na(country)) %>%
      group_by(country) %>%
      summarise(n_samples = n(), .groups = "drop") %>%
      mutate(total_minutes = n_samples / sample_rate_hz / 60) %>%
      arrange(desc(total_minutes)) %>%
      select(country, total_minutes)
  }
  
  return(totals)
}

get_country_primary_timezone <- function(data) {
  #' Get the most common timezone for each country (for display purposes)
  #'
  #' @param data Data frame with country and tz columns
  #' @return Named vector: country -> primary timezone name
  
  tz_by_country <- data %>%
    filter(!is.na(country) & !is.na(tz)) %>%
    group_by(country, tz) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(country) %>%
    slice_max(n, n = 1) %>%
    ungroup()
  
  setNames(tz_by_country$tz, tz_by_country$country)
}

# ------------------------------------------------------------------------------
# PLOTTING FUNCTIONS
# ------------------------------------------------------------------------------

create_dow_heatmap <- function(agg_data, country_order, config, title, subtitle = NULL) {
  #' Create faceted day-of-week heatmap
  
  day_labels <- c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
  
  # Limit to max_countries
  countries_to_plot <- head(country_order, config$max_countries)
  
  # Filter and prepare data
  plot_data <- agg_data %>%
    filter(country %in% countries_to_plot) %>%
    mutate(
      country = factor(country, levels = countries_to_plot),
      day_of_week = factor(day_of_week, levels = 1:7, labels = day_labels)
    )
  
  # Complete the grid (fill missing day/country combinations with 0)
  plot_data <- plot_data %>%
    complete(country, day_of_week, fill = list(minutes = 0))
  
  # Calculate totals for facet labels
  country_totals <- plot_data %>%
    group_by(country) %>%
    summarise(total = sum(minutes), .groups = "drop")
  
  # Create facet labels with totals
  facet_labels <- setNames(
    paste0(country_totals$country, "\n(", 
           ifelse(country_totals$total >= 60, 
                  paste0(round(country_totals$total / 60, 1), " hrs"),
                  paste0(round(country_totals$total, 0), " min")),
           ")"),
    country_totals$country
  )
  
  max_minutes <- max(plot_data$minutes, na.rm = TRUE)
  
  p <- ggplot(plot_data, aes(x = day_of_week, y = 1, fill = minutes)) +
    geom_tile(color = config$grid_color, linewidth = 0.5) +
    geom_text(aes(label = ifelse(minutes > 0, 
                                  ifelse(minutes >= 1, round(minutes), 
                                         sprintf("%.1f", minutes)), 
                                  "")),
              color = config$text_color, size = 3) +
    facet_wrap(~ country, ncol = 4, labeller = labeller(country = facet_labels)) +
    scale_fill_gradientn(
      colors = config$heatmap_colors,
      limits = c(0, max_minutes),
      name = "Minutes"
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = NULL
    ) +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = config$background_color, color = NA),
      panel.background = element_rect(fill = config$background_color, color = NA),
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "#2a2a2a", color = NA),
      strip.text = element_text(color = config$text_color, face = "bold", size = 10),
      axis.text.x = element_text(color = config$text_color, size = 9),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      legend.background = element_rect(fill = config$background_color, color = NA),
      legend.text = element_text(color = config$text_color),
      legend.title = element_text(color = config$text_color),
      plot.title = element_text(color = config$text_color, size = 14, face = "bold"),
      plot.subtitle = element_text(color = "#aaaaaa", size = 10)
    )
  
  return(p)
}

create_tod_heatmap <- function(agg_data, country_order, country_tz_map, config, title, subtitle = NULL) {
  #' Create faceted time-of-day heatmap
  #' Each country is shown in its LOCAL timezone
  #'
  #' @param agg_data Aggregated data by country and time_slot
  #' @param country_order Vector of country codes in display order
  #' @param country_tz_map Named vector mapping country -> timezone name
  #' @param config Configuration list
  #' @param title Plot title
  #' @param subtitle Plot subtitle
  
  time_labels <- sprintf("%02d:00", seq(0, 24, by = 4))
  time_breaks <- seq(0, 48, by = 8)
  
  countries_to_plot <- head(country_order, config$max_countries)
  
  plot_data <- agg_data %>%
    filter(country %in% countries_to_plot) %>%
    mutate(country = factor(country, levels = countries_to_plot))
  
  # Complete the grid
  all_slots <- expand.grid(
    country = factor(countries_to_plot, levels = countries_to_plot),
    time_slot = 0:47
  )
  
  plot_data <- all_slots %>%
    left_join(plot_data, by = c("country", "time_slot")) %>%
    replace_na(list(minutes = 0))
  
  # Calculate totals and get timezone info for facet labels
  country_info <- plot_data %>%
    group_by(country) %>%
    summarise(total = sum(minutes), .groups = "drop") %>%
    mutate(
      tz_name = country_tz_map[as.character(country)],
      # Simplify timezone name for display (e.g., "Asia/Tehran" -> "Tehran")
      tz_short = sapply(tz_name, function(tz) {
        if (is.na(tz)) return("UTC")
        parts <- strsplit(tz, "/")[[1]]
        tail(parts, 1)
      })
    )
  
  # Facet labels with total time AND timezone
  facet_labels <- setNames(
    paste0(country_info$country, " (", 
           ifelse(country_info$total >= 60, 
                  paste0(round(country_info$total / 60, 1), " hrs"),
                  paste0(round(country_info$total, 0), " min")),
           ")\n", country_info$tz_short),
    country_info$country
  )
  
  max_minutes <- max(plot_data$minutes, na.rm = TRUE)
  
  p <- ggplot(plot_data, aes(x = time_slot, y = 1, fill = minutes)) +
    geom_tile(color = NA) +
    facet_wrap(~ country, ncol = 3, labeller = labeller(country = facet_labels)) +
    scale_fill_gradientn(
      colors = config$heatmap_colors,
      limits = c(0, max_minutes),
      name = "Minutes"
    ) +
    scale_x_continuous(
      breaks = time_breaks,
      labels = time_labels,
      expand = c(0, 0)
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Time of Day (Local)",
      y = NULL
    ) +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = config$background_color, color = NA),
      panel.background = element_rect(fill = config$background_color, color = NA),
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "#2a2a2a", color = NA),
      strip.text = element_text(color = config$text_color, face = "bold", size = 10),
      axis.text.x = element_text(color = config$text_color, size = 8),
      axis.text.y = element_blank(),
      axis.title.x = element_text(color = config$text_color, size = 9),
      axis.ticks = element_blank(),
      legend.background = element_rect(fill = config$background_color, color = NA),
      legend.text = element_text(color = config$text_color),
      legend.title = element_text(color = config$text_color),
      plot.title = element_text(color = config$text_color, size = 14, face = "bold"),
      plot.subtitle = element_text(color = "#aaaaaa", size = 10)
    )
  
  return(p)
}

# ------------------------------------------------------------------------------
# FILENAME GENERATION
# ------------------------------------------------------------------------------

generate_temporal_filename <- function(config, analysis_type) {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  if (!is.null(config$mission_ids)) {
    prefix <- paste(config$mission_ids, collapse = "_")
  } else if (!is.null(config$date_start) || !is.null(config$date_end)) {
    start_date <- if (!is.null(config$date_start)) as.Date(config$date_start) else as.Date("1970-01-01")
    end_date <- if (!is.null(config$date_end)) as.Date(config$date_end) else Sys.Date()
    n_days <- as.integer(end_date - start_date)
    prefix <- paste0(n_days, "days")
  } else {
    prefix <- "all_missions"
  }
  
  filename <- paste0(prefix, "_", analysis_type, "_", timestamp, ".", config$output_format)
  
  if (!dir.exists(config$output_dir)) {
    dir.create(config$output_dir, recursive = TRUE)
  }
  
  return(file.path(config$output_dir, filename))
}

# ------------------------------------------------------------------------------
# MAIN EXECUTION
# ------------------------------------------------------------------------------

run_temporal_analysis <- function(config, col_map) {
  #' Run the complete temporal analysis
  #'
  #' @param config Configuration list
  #' @param col_map Column mapping list
  #' @return List with $plots (named list of ggplot objects) and $output_files (named list of paths)
  
  cli_h1("Temporal Dwell Time Analysis")
  
  # Initialize database
  con <- init_temporal_database(config, col_map)
  if (is.null(con)) return(NULL)
  
  results <- list(plots = list(), output_files = list())
  
  # Build subtitle from filter description
  if (!is.null(config$mission_ids)) {
    filter_desc <- paste("Missions:", paste(config$mission_ids, collapse = ", "))
  } else if (!is.null(config$date_start) && !is.null(config$date_end)) {
    filter_desc <- paste(config$date_start, "to", config$date_end)
  } else if (!is.null(config$date_start)) {
    filter_desc <- paste("From", config$date_start)
  } else if (!is.null(config$date_end)) {
    filter_desc <- paste("Until", config$date_end)
  } else {
    filter_desc <- "All missions"
  }
  
  # ============================================================================
  # 1. VUEFAST DATA (Aircraft Position)
  # ============================================================================
  
  cli_h2("Processing Vuefast Data (Aircraft Position)")
  
  vuefast_data <- query_vuefast_temporal(con, config, col_map)
  
  if (nrow(vuefast_data) > 0) {
    
    # Assign countries and timezones
    vuefast_data <- assign_countries_and_timezones(vuefast_data, config)
    
    # Get primary timezone per country for display
    country_tz_map <- get_country_primary_timezone(vuefast_data)
    
    # Add local time components (day of week, time slot with DST handling)
    vuefast_data <- add_local_time_components(vuefast_data, "timestamp_utc")
    
    # Get country order by total time
    country_totals <- get_country_totals(vuefast_data, is_clipmarks = FALSE)
    country_order <- country_totals$country
    
    if (length(country_order) > 0) {
      
      # --- 1a: Day of Week ---
      cli_h3("1a. Day of Week - Aircraft Position")
      
      dow_agg <- aggregate_by_country_dow(vuefast_data, is_clipmarks = FALSE)
      
      if (nrow(dow_agg) > 0) {
        p_dow <- create_dow_heatmap(
          dow_agg, 
          country_order, 
          config,
          title = "Dwell Time by Day of Week - Aircraft Position",
          subtitle = filter_desc
        )
        
        output_file <- generate_temporal_filename(config, "dow_vuefast")
        ggsave(output_file, plot = p_dow, 
               width = config$output_width, height = config$output_height, 
               dpi = config$output_dpi)
        cli_alert_success("Saved: {output_file}")
        
        results$plots$dow_vuefast <- p_dow
        results$output_files$dow_vuefast <- output_file
      }
      
      # --- 2a: Time of Day ---
      cli_h3("2a. Time of Day - Aircraft Position (Local Timezones with DST)")
      
      tod_agg <- aggregate_by_country_tod(vuefast_data, is_clipmarks = FALSE)
      
      if (nrow(tod_agg) > 0) {
        p_tod <- create_tod_heatmap(
          tod_agg,
          country_order,
          country_tz_map,
          config,
          title = "Dwell Time by Time of Day (30-min slots) - Aircraft Position",
          subtitle = paste0(filter_desc, " | Local time with DST")
        )
        
        output_file <- generate_temporal_filename(config, "tod_vuefast")
        ggsave(output_file, plot = p_tod,
               width = config$output_width, height = config$output_height,
               dpi = config$output_dpi)
        cli_alert_success("Saved: {output_file}")
        
        results$plots$tod_vuefast <- p_tod
        results$output_files$tod_vuefast <- output_file
      }
    }
  } else {
    cli_alert_warning("No vuefast data found for query parameters")
  }
  
  # ============================================================================
  # 2. CLIPMARKS DATA (Collections)
  # ============================================================================
  
  cli_h2("Processing Clipmarks Data (Collections)")
  
  clipmarks_data <- query_clipmarks_temporal(con, config, col_map)
  
  if (!is.null(clipmarks_data) && nrow(clipmarks_data) > 0) {
    
    # Assign countries and timezones
    clipmarks_data <- assign_countries_and_timezones(clipmarks_data, config)
    
    # Get primary timezone per country for display
    country_tz_map <- get_country_primary_timezone(clipmarks_data)
    
    # Add local time components
    clipmarks_data <- add_local_time_components(clipmarks_data, "start_time_utc")
    
    # Get country order
    country_totals <- get_country_totals(clipmarks_data, is_clipmarks = TRUE)
    country_order <- country_totals$country
    
    if (length(country_order) > 0) {
      
      # --- 1b: Day of Week ---
      cli_h3("1b. Day of Week - Clipmarks")
      
      dow_agg <- aggregate_by_country_dow(clipmarks_data, is_clipmarks = TRUE)
      
      if (nrow(dow_agg) > 0) {
        p_dow <- create_dow_heatmap(
          dow_agg,
          country_order,
          config,
          title = "Collection Time by Day of Week - Clipmarks",
          subtitle = filter_desc
        )
        
        output_file <- generate_temporal_filename(config, "dow_clipmarks")
        ggsave(output_file, plot = p_dow,
               width = config$output_width, height = config$output_height,
               dpi = config$output_dpi)
        cli_alert_success("Saved: {output_file}")
        
        results$plots$dow_clipmarks <- p_dow
        results$output_files$dow_clipmarks <- output_file
      }
      
      # --- 2b: Time of Day ---
      cli_h3("2b. Time of Day - Clipmarks (Local Timezones with DST)")
      
      tod_agg <- aggregate_by_country_tod(clipmarks_data, is_clipmarks = TRUE)
      
      if (nrow(tod_agg) > 0) {
        p_tod <- create_tod_heatmap(
          tod_agg,
          country_order,
          country_tz_map,
          config,
          title = "Collection Time by Time of Day (30-min slots) - Clipmarks",
          subtitle = paste0(filter_desc, " | Local time with DST")
        )
        
        output_file <- generate_temporal_filename(config, "tod_clipmarks")
        ggsave(output_file, plot = p_tod,
               width = config$output_width, height = config$output_height,
               dpi = config$output_dpi)
        cli_alert_success("Saved: {output_file}")
        
        results$plots$tod_clipmarks <- p_tod
        results$output_files$tod_clipmarks <- output_file
      }
    }
  } else {
    cli_alert_warning("No clipmarks data found for query parameters")
  }
  
  # Cleanup
  dbDisconnect(con)
  
  cli_h1("Analysis Complete")
  
  if (length(results$output_files) > 0) {
    cli_alert_success("Generated {length(results$output_files)} images:")
    for (name in names(results$output_files)) {
      cat(paste0("  - ", name, ": ", results$output_files[[name]], "\n"))
    }
  }
  
  return(results)
}

# ==============================================================================
# RUN
# ==============================================================================

cat("\n")
cli_alert_info("Temporal Analysis script loaded. To run:")
cat("
1. Ensure test data exists (or run generate_synthetic_data() from hexbin_dwell_analysis.R)

2. Adjust 'config' and 'col_map' settings at top of script

3. Run analysis:
   result <- run_temporal_analysis(config, col_map)

4. Display plots:
   print(result$plots$dow_vuefast)    # Day of week - aircraft position
   print(result$plots$tod_vuefast)    # Time of day - aircraft position
   print(result$plots$dow_clipmarks)  # Day of week - collections
   print(result$plots$tod_clipmarks)  # Time of day - collections

5. Check output files:
   result$output_files

TIMEZONE HANDLING:
  - Uses 'lutz' package to look up timezone from coordinates
  - Uses 'lubridate::with_tz()' for DST-aware local time conversion
  - Timezone lookup method configurable: 'fast' or 'accurate'
  - Each facet shows the primary timezone name for that country
")
