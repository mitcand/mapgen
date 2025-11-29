# ==============================================================================
# GIMBAL SLEW RATE ANALYSIS
# ==============================================================================
#
# Analyzes gimbal pointing behavior from aircraft camera metadata logs.
# Generates:
#   1. Overall slew rate charts (azimuth and elevation rate of change)
#   2. On-target slew rate charts (filtered by clipmark windows)
#   3. On-target pointing hotspot (where gimbal spends most time while on target)
#
# ==============================================================================

# ------------------------------------------------------------------------------
# CONFIGURATION
# ------------------------------------------------------------------------------

config <- list(
  # === DATA SOURCE ===
  data_dir = "data/missions",
  db_path = "missions.duckdb",
  use_existing_db = FALSE,
  
  # === QUERY FILTERS ===
  mission_ids = NULL,
  date_start = NULL,
  date_end = NULL,
  
  # === ANALYSIS PARAMETERS ===
  sample_rate_hz = 15,                    # Data sample rate for rate calculation
  slew_rate_units = "deg_per_sec",        # "deg_per_sec" or "deg_per_frame"
  
  # Hotspot settings
 hotspot_az_bins = 72,                   # Number of azimuth bins (72 = 5° bins)
  hotspot_el_bins = 36,                   # Number of elevation bins (36 = 5° bins)
  
  # === OUTPUT ===
  output_dir = "output",
  output_format = "png",
  output_width = 12,
  output_height = 8,
  output_dpi = 300
)

# ------------------------------------------------------------------------------
# COLUMN MAPPING - Update these to match your CSV column names
# ------------------------------------------------------------------------------

col_map <- list(
  # Timestamp
  timestamp = "timestamp",
  
  # Mission identifier
  mission_id = "mission_id",
  
  # Gimbal pointing (sensor-relative or absolute)
  gimbal_azimuth = "sensor_az",           # Azimuth in degrees
  gimbal_elevation = "sensor_el",         # Elevation in degrees (negative = looking down)
  
  # Aircraft position (for context)
  aircraft_lat = "ac_lat",
  aircraft_lon = "ac_lon",
  aircraft_alt = "ac_alt",
  aircraft_heading = "ac_heading"
)

# Clipmark column mapping
clip_col_map <- list(
  mission_id = "mission_id",
  start_time = "start_time",
  stop_time = "stop_time",
  target_lat = "target_lat",
  target_lon = "target_lon"
)

# ------------------------------------------------------------------------------
# PACKAGES
# ------------------------------------------------------------------------------

required_packages <- c(
  "DBI", "duckdb", "dplyr", "lubridate", "ggplot2", "glue", "cli", "viridis", "patchwork", "zoo", "scales"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# ------------------------------------------------------------------------------
# DATABASE FUNCTIONS
# ------------------------------------------------------------------------------

init_database <- function(config, col_map) {
  cli_h2("Initializing Database")
  
  if (config$use_existing_db && file.exists(config$db_path)) {
    cli_alert_success("Using existing database: {config$db_path}")
    con <- dbConnect(duckdb(), config$db_path)
    return(con)
  }
  
  # Create new database
  if (file.exists(config$db_path)) {
    file.remove(config$db_path)
  }
  
  con <- dbConnect(duckdb(), config$db_path)
  
  # Find and import CSV files
  csv_files <- list.files(config$data_dir, pattern = "_vuefast\\.csv$", 
                          full.names = TRUE, recursive = TRUE)
  
  if (length(csv_files) == 0) {
    cli_alert_danger("No vuefast CSV files found in {config$data_dir}")
    dbDisconnect(con)
    return(NULL)
  }
  
  cli_alert_info("Found {length(csv_files)} vuefast files")
  
  # Import each file
  for (i in seq_along(csv_files)) {
    file_path <- csv_files[i]
    mission_id <- gsub("_vuefast\\.csv$", "", basename(file_path))
    
    cli_alert("Importing {basename(file_path)}...")
    
    tryCatch({
      df <- read.csv(file_path, stringsAsFactors = FALSE)
      df$mission_id <- mission_id
      
      if (i == 1) {
        dbWriteTable(con, "vuefast", df, overwrite = TRUE)
      } else {
        dbWriteTable(con, "vuefast", df, append = TRUE)
      }
    }, error = function(e) {
      cli_alert_warning("Failed to import {basename(file_path)}: {e$message}")
    })
  }
  
  # Import clipmarks
  clip_files <- list.files(config$data_dir, pattern = "_clipmarks\\.csv$",
                           full.names = TRUE, recursive = TRUE)
  
  if (length(clip_files) > 0) {
    cli_alert_info("Found {length(clip_files)} clipmark files")
    
    for (i in seq_along(clip_files)) {
      file_path <- clip_files[i]
      mission_id <- gsub("_clipmarks\\.csv$", "", basename(file_path))
      
      tryCatch({
        df <- read.csv(file_path, stringsAsFactors = FALSE)
        df$mission_id <- mission_id
        
        if (i == 1) {
          dbWriteTable(con, "clipmarks", df, overwrite = TRUE)
        } else {
          dbWriteTable(con, "clipmarks", df, append = TRUE)
        }
      }, error = function(e) {
        cli_alert_warning("Failed to import clipmarks: {e$message}")
      })
    }
  }
  
  cli_alert_success("Database initialized")
  return(con)
}

# ------------------------------------------------------------------------------
# QUERY FUNCTIONS
# ------------------------------------------------------------------------------

query_gimbal_data <- function(con, config, col_map) {
  # Build column selection
  cols <- glue("{col_map$timestamp} as timestamp,
                {col_map$mission_id} as mission_id,
                {col_map$gimbal_azimuth} as azimuth,
                {col_map$gimbal_elevation} as elevation")
  
  # Build WHERE clause
  where_parts <- c()
  
  if (!is.null(config$mission_ids)) {
    ids <- paste(sprintf("'%s'", config$mission_ids), collapse = ", ")
    where_parts <- c(where_parts, glue("{col_map$mission_id} IN ({ids})"))
  }
  
  where_clause <- if (length(where_parts) > 0) {
    paste("WHERE", paste(where_parts, collapse = " AND "))
  } else {
    ""
  }
  
  query <- glue("SELECT {cols} FROM vuefast {where_clause} ORDER BY timestamp")
  
  cli_alert("Querying gimbal data...")
  data <- dbGetQuery(con, query)
  cli_alert_success("Retrieved {format(nrow(data), big.mark=',')} rows")
  
  return(data)
}

query_clipmarks <- function(con, config, clip_col_map) {
  # Check if clipmarks table exists
  tables <- dbListTables(con)
  if (!"clipmarks" %in% tables) {
    cli_alert_warning("No clipmarks table found")
    return(NULL)
  }
  
  cols <- glue("{clip_col_map$mission_id} as mission_id,
                {clip_col_map$start_time} as start_time,
                {clip_col_map$stop_time} as stop_time")
  
  # Build WHERE clause
  where_parts <- c()
  
  if (!is.null(config$mission_ids)) {
    ids <- paste(sprintf("'%s'", config$mission_ids), collapse = ", ")
    where_parts <- c(where_parts, glue("{clip_col_map$mission_id} IN ({ids})"))
  }
  
  where_clause <- if (length(where_parts) > 0) {
    paste("WHERE", paste(where_parts, collapse = " AND "))
  } else {
    ""
  }
  
  query <- glue("SELECT {cols} FROM clipmarks {where_clause}")
  
  cli_alert("Querying clipmarks...")
  data <- dbGetQuery(con, query)
  cli_alert_success("Retrieved {nrow(data)} clipmarks")
  
  return(data)
}

# ------------------------------------------------------------------------------
# ANALYSIS FUNCTIONS
# ------------------------------------------------------------------------------
calculate_slew_rates <- function(gimbal_data, config) {
  cli_alert("Calculating slew rates...")
  
  # Sort by timestamp
  gimbal_data <- gimbal_data %>%
    arrange(mission_id, timestamp)
  
  # Calculate differences
  gimbal_data <- gimbal_data %>%
    group_by(mission_id) %>%
    mutate(
      # Time delta
      dt = as.numeric(difftime(timestamp, lag(timestamp), units = "secs")),
      
      # Angle deltas (handle wrap-around for azimuth)
      d_az = azimuth - lag(azimuth),
      d_el = elevation - lag(elevation),
      
      # Fix azimuth wrap-around (e.g., 359 -> 1 should be +2, not -358)
      d_az = case_when(
        d_az > 180 ~ d_az - 360,
        d_az < -180 ~ d_az + 360,
        TRUE ~ d_az
      ),
      
      # Slew rates
      az_rate = case_when(
        dt > 0 & dt < 1 ~ d_az / dt,  # deg/sec (filter out large gaps)
        TRUE ~ NA_real_
      ),
      el_rate = case_when(
        dt > 0 & dt < 1 ~ d_el / dt,
        TRUE ~ NA_real_
      ),
      
      # Absolute slew rate (total angular velocity)
      total_rate = sqrt(az_rate^2 + el_rate^2)
    ) %>%
    ungroup()
  
  # Remove first row of each mission (no rate calculable)
  gimbal_data <- gimbal_data %>%
    filter(!is.na(az_rate))
  
  n_valid <- sum(!is.na(gimbal_data$total_rate))
  cli_alert_success("Calculated {format(n_valid, big.mark=',')} valid slew rate samples")
  
  return(gimbal_data)
}

filter_on_target <- function(gimbal_data, clipmarks) {
  if (is.null(clipmarks) || nrow(clipmarks) == 0) {
    cli_alert_warning("No clipmarks available for on-target filtering")
    return(NULL)
  }
  
  cli_alert("Filtering to on-target periods...")
  
  # Convert timestamps
  gimbal_data$timestamp <- as.POSIXct(gimbal_data$timestamp)
  clipmarks$start_time <- as.POSIXct(clipmarks$start_time)
  clipmarks$stop_time <- as.POSIXct(clipmarks$stop_time)
  
  # Flag rows that fall within any clipmark window
  gimbal_data$on_target <- FALSE
  
  for (i in seq_len(nrow(clipmarks))) {
    clip <- clipmarks[i, ]
    in_window <- gimbal_data$mission_id == clip$mission_id &
                 gimbal_data$timestamp >= clip$start_time &
                 gimbal_data$timestamp <= clip$stop_time
    gimbal_data$on_target[in_window] <- TRUE
  }
  
  on_target_data <- gimbal_data %>% filter(on_target == TRUE)
  
  n_on_target <- nrow(on_target_data)
  pct_on_target <- round(100 * n_on_target / nrow(gimbal_data), 1)
  
  cli_alert_success("Found {format(n_on_target, big.mark=',')} on-target samples ({pct_on_target}% of data)")
  
  return(on_target_data)
}

# ------------------------------------------------------------------------------
# PLOTTING FUNCTIONS
# ------------------------------------------------------------------------------

plot_slew_rate_histogram <- function(data, title_suffix = "", config) {
  # Filter extreme values for visualization
  data_filtered <- data %>%
    filter(abs(az_rate) < 100, abs(el_rate) < 100)
  
  # Azimuth rate histogram
  p_az <- ggplot(data_filtered, aes(x = az_rate)) +
    geom_histogram(bins = 100, fill = "#4575b4", color = NA, alpha = 0.8) +
    scale_y_continuous(labels = scales::comma) +
    labs(
      title = paste("Azimuth Slew Rate", title_suffix),
      x = "Rate (°/sec)",
      y = "Count"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  # Elevation rate histogram
  p_el <- ggplot(data_filtered, aes(x = el_rate)) +
    geom_histogram(bins = 100, fill = "#d73027", color = NA, alpha = 0.8) +
    scale_y_continuous(labels = scales::comma) +
    labs(
      title = paste("Elevation Slew Rate", title_suffix),
      x = "Rate (°/sec)",
      y = "Count"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  # Total slew rate
  p_total <- ggplot(data_filtered, aes(x = total_rate)) +
    geom_histogram(bins = 100, fill = "#1a9850", color = NA, alpha = 0.8) +
    scale_y_continuous(labels = scales::comma) +
    labs(
      title = paste("Total Slew Rate", title_suffix),
      x = "Rate (°/sec)",
      y = "Count"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  # Combine
  combined <- p_az / p_el / p_total
  
  return(combined)
}

plot_slew_rate_time_series <- function(data, title_suffix = "", config) {
  # Downsample for plotting if too many points
  if (nrow(data) > 50000) {
    data <- data %>% sample_n(50000) %>% arrange(timestamp)
  }
  
  data$timestamp <- as.POSIXct(data$timestamp)
  
  # Clamp extreme values for visualization
  data <- data %>%
    mutate(
      az_rate_clamped = pmax(pmin(az_rate, 50), -50),
      el_rate_clamped = pmax(pmin(el_rate, 50), -50)
    )
  
  p <- ggplot(data, aes(x = timestamp)) +
    geom_line(aes(y = az_rate_clamped, color = "Azimuth"), alpha = 0.5, linewidth = 0.3) +
    geom_line(aes(y = el_rate_clamped, color = "Elevation"), alpha = 0.5, linewidth = 0.3) +
    scale_color_manual(values = c("Azimuth" = "#4575b4", "Elevation" = "#d73027")) +
    labs(
      title = paste("Slew Rate Over Time", title_suffix),
      x = "Time",
      y = "Rate (°/sec)",
      color = "Axis"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "bottom"
    )
  
  return(p)
}

plot_pointing_hotspot <- function(data, title_suffix = "", config) {
  cli_alert("Generating pointing hotspot...")
  
  # Bin the azimuth and elevation data
  data <- data %>%
    mutate(
      az_bin = cut(azimuth, breaks = seq(0, 360, length.out = config$hotspot_az_bins + 1), 
                   include.lowest = TRUE, labels = FALSE),
      el_bin = cut(elevation, breaks = seq(-90, 90, length.out = config$hotspot_el_bins + 1),
                   include.lowest = TRUE, labels = FALSE)
    ) %>%
    filter(!is.na(az_bin), !is.na(el_bin))
  
  # Count time in each bin
  hotspot <- data %>%
    group_by(az_bin, el_bin) %>%
    summarise(
      count = n(),
      dwell_sec = n() / config$sample_rate_hz,
      .groups = "drop"
    ) %>%
    mutate(
      # Convert bin indices back to center values
      azimuth = (az_bin - 0.5) * (360 / config$hotspot_az_bins),
      elevation = (el_bin - 0.5) * (180 / config$hotspot_el_bins) - 90
    )
  
  # Create plot
  p <- ggplot(hotspot, aes(x = azimuth, y = elevation, fill = dwell_sec)) +
    geom_tile() +
    scale_fill_viridis(
      option = "inferno",
      name = "Dwell (sec)",
      labels = scales::comma
    ) +
    scale_x_continuous(
      breaks = seq(0, 360, 45),
      labels = c("0°\n(Fwd)", "45°", "90°\n(Right)", "135°", "180°\n(Aft)", 
                 "225°", "270°\n(Left)", "315°", "360°")
    ) +
    scale_y_continuous(
      breaks = seq(-90, 90, 30),
      labels = c("-90°\n(Down)", "-60°", "-30°", "0°\n(Horizon)", 
                 "30°", "60°", "90°\n(Up)")
    ) +
    labs(
      title = paste("Gimbal Pointing Hotspot", title_suffix),
      subtitle = "Time spent at each azimuth/elevation combination",
      x = "Azimuth (relative to aircraft heading)",
      y = "Elevation"
    ) +
    coord_fixed(ratio = 2) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
  
  return(p)
}

plot_slew_vs_pointing <- function(data, title_suffix = "", config) {
  # Show slew rate by pointing direction
  cli_alert("Generating slew rate by pointing direction...")
  
  # Bin the data
  data <- data %>%
    mutate(
      az_bin = cut(azimuth, breaks = seq(0, 360, length.out = config$hotspot_az_bins + 1), 
                   include.lowest = TRUE, labels = FALSE),
      el_bin = cut(elevation, breaks = seq(-90, 90, length.out = config$hotspot_el_bins + 1),
                   include.lowest = TRUE, labels = FALSE)
    ) %>%
    filter(!is.na(az_bin), !is.na(el_bin), !is.na(total_rate))
  
  # Mean slew rate in each bin
  slew_by_pointing <- data %>%
    group_by(az_bin, el_bin) %>%
    summarise(
      mean_rate = mean(total_rate, na.rm = TRUE),
      median_rate = median(total_rate, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      azimuth = (az_bin - 0.5) * (360 / config$hotspot_az_bins),
      elevation = (el_bin - 0.5) * (180 / config$hotspot_el_bins) - 90
    )
  
  p <- ggplot(slew_by_pointing, aes(x = azimuth, y = elevation, fill = median_rate)) +
    geom_tile() +
    scale_fill_viridis(
      option = "plasma",
      name = "Median Rate\n(°/sec)",
      limits = c(0, quantile(slew_by_pointing$median_rate, 0.95, na.rm = TRUE))
    ) +
    scale_x_continuous(breaks = seq(0, 360, 45)) +
    scale_y_continuous(breaks = seq(-90, 90, 30)) +
    labs(
      title = paste("Slew Rate by Pointing Direction", title_suffix),
      subtitle = "Median slew rate at each azimuth/elevation",
      x = "Azimuth (°)",
      y = "Elevation (°)"
    ) +
    coord_fixed(ratio = 2) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold")
    )
  
  return(p)
}

# ------------------------------------------------------------------------------
# MAIN ANALYSIS
# ------------------------------------------------------------------------------

run_gimbal_analysis <- function(config, col_map, clip_col_map) {
  cli_h1("Gimbal Slew Rate Analysis")
  
  # Initialize database
  con <- init_database(config, col_map)
  if (is.null(con)) return(NULL)
  
  # Query data
  gimbal_data <- query_gimbal_data(con, config, col_map)
  clipmarks <- query_clipmarks(con, config, clip_col_map)
  
  dbDisconnect(con)
  
  if (nrow(gimbal_data) == 0) {
    cli_alert_danger("No gimbal data found")
    return(NULL)
  }
  
  # Parse timestamps
  gimbal_data$timestamp <- as.POSIXct(gimbal_data$timestamp)
  
  # Calculate slew rates
  gimbal_data <- calculate_slew_rates(gimbal_data, config)
  
  # Filter to on-target data
  on_target_data <- filter_on_target(gimbal_data, clipmarks)
  
  # Create output directory
  dir.create(config$output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Generate plots
  results <- list()
  
  # 1. Overall slew rate histogram
  cli_h2("Generating Overall Slew Rate Charts")
  p_hist_all <- plot_slew_rate_histogram(gimbal_data, " - All Data", config)
  
  hist_file <- file.path(config$output_dir, paste0("slew_rate_histogram_all.", config$output_format))
  ggsave(hist_file, p_hist_all, width = config$output_width, height = config$output_height, dpi = config$output_dpi)
  cli_alert_success("Saved: {hist_file}")
  results$histogram_all <- p_hist_all
  
  # 2. Overall time series
  p_ts_all <- plot_slew_rate_time_series(gimbal_data, " - All Data", config)
  
  ts_file <- file.path(config$output_dir, paste0("slew_rate_timeseries_all.", config$output_format))
  ggsave(ts_file, p_ts_all, width = config$output_width, height = config$output_height / 2, dpi = config$output_dpi)
  cli_alert_success("Saved: {ts_file}")
  results$timeseries_all <- p_ts_all
  
  # 3. On-target analysis (if clipmarks available)
  if (!is.null(on_target_data) && nrow(on_target_data) > 0) {
    cli_h2("Generating On-Target Charts")
    
    # On-target histogram
    p_hist_target <- plot_slew_rate_histogram(on_target_data, " - On Target", config)
    
    hist_target_file <- file.path(config$output_dir, paste0("slew_rate_histogram_ontarget.", config$output_format))
    ggsave(hist_target_file, p_hist_target, width = config$output_width, height = config$output_height, dpi = config$output_dpi)
    cli_alert_success("Saved: {hist_target_file}")
    results$histogram_ontarget <- p_hist_target
    
    # On-target pointing hotspot
    p_hotspot <- plot_pointing_hotspot(on_target_data, " - On Target", config)
    
    hotspot_file <- file.path(config$output_dir, paste0("pointing_hotspot_ontarget.", config$output_format))
    ggsave(hotspot_file, p_hotspot, width = config$output_width, height = config$output_height, dpi = config$output_dpi)
    cli_alert_success("Saved: {hotspot_file}")
    results$hotspot_ontarget <- p_hotspot
    
    # Slew rate by pointing direction (on-target)
    p_slew_pointing <- plot_slew_vs_pointing(on_target_data, " - On Target", config)
    
    slew_pointing_file <- file.path(config$output_dir, paste0("slew_by_pointing_ontarget.", config$output_format))
    ggsave(slew_pointing_file, p_slew_pointing, width = config$output_width, height = config$output_height, dpi = config$output_dpi)
    cli_alert_success("Saved: {slew_pointing_file}")
    results$slew_by_pointing <- p_slew_pointing
  }
  
  # 4. Overall pointing hotspot (for comparison)
  cli_h2("Generating Overall Pointing Hotspot")
  p_hotspot_all <- plot_pointing_hotspot(gimbal_data, " - All Data", config)
  
  hotspot_all_file <- file.path(config$output_dir, paste0("pointing_hotspot_all.", config$output_format))
  ggsave(hotspot_all_file, p_hotspot_all, width = config$output_width, height = config$output_height, dpi = config$output_dpi)
  cli_alert_success("Saved: {hotspot_all_file}")
  results$hotspot_all <- p_hotspot_all
  
  # Summary statistics
  cli_h2("Summary Statistics")
  
  cli_alert_info("Overall:")
  cli_alert_info("  Median azimuth rate: {round(median(gimbal_data$az_rate, na.rm=TRUE), 2)}°/sec")
  cli_alert_info("  Median elevation rate: {round(median(gimbal_data$el_rate, na.rm=TRUE), 2)}°/sec")
  cli_alert_info("  Median total rate: {round(median(gimbal_data$total_rate, na.rm=TRUE), 2)}°/sec")
  
  if (!is.null(on_target_data) && nrow(on_target_data) > 0) {
    cli_alert_info("On-Target:")
    cli_alert_info("  Median azimuth rate: {round(median(on_target_data$az_rate, na.rm=TRUE), 2)}°/sec")
    cli_alert_info("  Median elevation rate: {round(median(on_target_data$el_rate, na.rm=TRUE), 2)}°/sec")
    cli_alert_info("  Median total rate: {round(median(on_target_data$total_rate, na.rm=TRUE), 2)}°/sec")
  }
  
  cli_h1("Analysis Complete")
  
  return(results)
}

# ==============================================================================
# RUN
# ==============================================================================

# Uncomment to run:
# results <- run_gimbal_analysis(config, col_map, clip_col_map)
