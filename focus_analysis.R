# ==============================================================================
# FOCUS HUNTING DETECTION - PHASE 1: ANALYSIS
# ==============================================================================
#
# Analyzes camera focus behavior to detect:
#   - Autofocus events (full 0-100% sweeps)
#   - Manual focus adjustments
#   - "Locked in" periods (stable focus while on target)
#
# Outputs:
#   - Focus event timeline charts (faceted by camera/zoom)
#   - focus_events.csv for Phase 2 video frame extraction
#   - Summary statistics
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
  sample_rate_hz = 15,
  
 # Autofocus detection
  autofocus_sweep_threshold = 70,        # Minimum range (%) to consider autofocus
  autofocus_time_window_sec = 3,         # Max time for a sweep to complete
  
  # Focus change detection
  focus_change_threshold = 5,            # Minimum change (%) to flag as event
  focus_stable_threshold = 1,            # Max variation (%) to consider "locked"
  focus_stable_duration_sec = 2,         # Min duration to consider "locked"
  
  # Event filtering
  min_event_gap_sec = 1,                 # Merge events closer than this
  
  # === OUTPUT ===
  output_dir = "output",
  output_format = "png",
  output_width = 14,
  output_height = 10,
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
  
  # Camera 1 focus values (0-100%)
  cam1_focus = "cam1_focus",
  
  # Camera 2 focus values (0-100%)
  cam2_focus = "cam2_focus",
  
  # Zoom levels (0%, 50%, 100%)
  cam1_zoom = "cam1_zoom",
  cam2_zoom = "cam2_zoom"
)

# Clipmark column mapping
clip_col_map <- list(
  mission_id = "mission_id",
  start_time = "start_time",
  stop_time = "stop_time"
)

# ------------------------------------------------------------------------------
# PACKAGES
# ------------------------------------------------------------------------------

required_packages <- c(
  "DBI", "duckdb", "dplyr", "tidyr", "lubridate", 
  "ggplot2", "glue", "cli", "viridis", "patchwork", "zoo", "scales"
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
  
  if (file.exists(config$db_path)) {
    file.remove(config$db_path)
  }
  
  con <- dbConnect(duckdb(), config$db_path)
  
  csv_files <- list.files(config$data_dir, pattern = "_vuefast\\.csv$", 
                          full.names = TRUE, recursive = TRUE)
  
  if (length(csv_files) == 0) {
    cli_alert_danger("No vuefast CSV files found in {config$data_dir}")
    dbDisconnect(con)
    return(NULL)
  }
  
  cli_alert_info("Found {length(csv_files)} vuefast files")
  
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

query_focus_data <- function(con, config, col_map) {
  cols <- glue("{col_map$timestamp} as timestamp,
                {col_map$mission_id} as mission_id,
                {col_map$cam1_focus} as cam1_focus,
                {col_map$cam2_focus} as cam2_focus,
                {col_map$cam1_zoom} as cam1_zoom,
                {col_map$cam2_zoom} as cam2_zoom")
  
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
  
  cli_alert("Querying focus data...")
  data <- dbGetQuery(con, query)
  cli_alert_success("Retrieved {format(nrow(data), big.mark=',')} rows")
  
  return(data)
}

query_clipmarks <- function(con, config, clip_col_map) {
  tables <- dbListTables(con)
  if (!"clipmarks" %in% tables) {
    cli_alert_warning("No clipmarks table found")
    return(NULL)
  }
  
  cols <- glue("{clip_col_map$mission_id} as mission_id,
                {clip_col_map$start_time} as start_time,
                {clip_col_map$stop_time} as stop_time")
  
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
  
  data <- dbGetQuery(con, query)
  return(data)
}

# ------------------------------------------------------------------------------
# ANALYSIS FUNCTIONS
# ------------------------------------------------------------------------------

reshape_focus_data <- function(focus_data) {
  # Convert wide format (cam1_focus, cam2_focus) to long format
  # with camera and zoom as factors
  
  cli_alert("Reshaping focus data to long format...")
  
  focus_data$timestamp <- as.POSIXct(focus_data$timestamp)
  
  # Create long format with camera, zoom, and focus
  cam1 <- focus_data %>%
    select(timestamp, mission_id, focus = cam1_focus, zoom = cam1_zoom) %>%
    mutate(camera = "Camera 1")
  
  cam2 <- focus_data %>%
    select(timestamp, mission_id, focus = cam2_focus, zoom = cam2_zoom) %>%
    mutate(camera = "Camera 2")
  
  long_data <- bind_rows(cam1, cam2) %>%
    arrange(camera, timestamp) %>%
    mutate(
      zoom_level = case_when(
        zoom == 0 ~ "Wide (0%)",
        zoom == 50 ~ "Mid (50%)",
        zoom == 100 ~ "Tele (100%)",
        TRUE ~ paste0(zoom, "%")
      ),
      zoom_level = factor(zoom_level, levels = c("Wide (0%)", "Mid (50%)", "Tele (100%)"))
    )
  
  cli_alert_success("Reshaped to {format(nrow(long_data), big.mark=',')} rows")
  
  return(long_data)
}

detect_focus_events <- function(focus_data, config) {
  cli_alert("Detecting focus events...")
  
  events <- list()
  
  # Process each camera separately
  for (cam in unique(focus_data$camera)) {
    cam_data <- focus_data %>%
      filter(camera == cam) %>%
      arrange(timestamp)
    
    # Calculate focus changes
    cam_data <- cam_data %>%
      group_by(mission_id) %>%
      mutate(
        focus_diff = focus - lag(focus),
        focus_diff = ifelse(is.na(focus_diff), 0, focus_diff),
        dt = as.numeric(difftime(timestamp, lag(timestamp), units = "secs")),
        dt = ifelse(is.na(dt), 0, dt),
        focus_rate = ifelse(dt > 0 & dt < 1, focus_diff / dt, 0)
      ) %>%
      ungroup()
    
    # Detect autofocus sweeps (large range in short time)
    cam_data <- detect_autofocus_sweeps(cam_data, config)
    
    # Detect manual focus adjustments
    cam_data <- detect_manual_adjustments(cam_data, config)
    
    # Detect stable "locked" periods
    cam_data <- detect_locked_periods(cam_data, config)
    
    events[[cam]] <- cam_data
  }
  
  combined <- bind_rows(events)
  
  n_autofocus <- sum(combined$is_autofocus_start, na.rm = TRUE)
  n_manual <- sum(combined$is_manual_adjust, na.rm = TRUE)
  n_locked <- sum(combined$locked_start, na.rm = TRUE)
  
  cli_alert_success("Detected: {n_autofocus} autofocus, {n_manual} manual adjustments, {n_locked} locked periods")
  
  return(combined)
}

detect_autofocus_sweeps <- function(data, config) {
  # Look for rapid sweeps covering most of the focus range
  
  window_samples <- config$autofocus_time_window_sec * config$sample_rate_hz
  
  data <- data %>%
    mutate(
      # Rolling max and min over window
      focus_roll_max = zoo::rollmax(focus, window_samples, fill = NA, align = "center"),
      focus_roll_min = zoo::rollmin(focus, window_samples, fill = NA, align = "center"),
      focus_roll_range = focus_roll_max - focus_roll_min,
      
      # Flag autofocus events
      is_autofocus = focus_roll_range >= config$autofocus_sweep_threshold,
      
      # Find start of autofocus sequences
      is_autofocus_start = is_autofocus & !lag(is_autofocus, default = FALSE)
    )
  
  return(data)
}

detect_manual_adjustments <- function(data, config) {
  # Detect significant focus changes that aren't autofocus
  
  data <- data %>%
    mutate(
      # Significant change in focus
      significant_change = abs(focus_diff) >= config$focus_change_threshold,
      
      # Manual adjust = significant change but not during autofocus
      is_manual_adjust = significant_change & !is_autofocus
    )
  
  return(data)
}

detect_locked_periods <- function(data, config) {
  # Detect periods where focus is stable
  
  window_samples <- config$focus_stable_duration_sec * config$sample_rate_hz
  
  data <- data %>%
    mutate(
      # Rolling standard deviation
      focus_roll_sd = zoo::rollapply(focus, window_samples, sd, fill = NA, align = "center"),
      
      # Locked = low variation
      is_locked = focus_roll_sd <= config$focus_stable_threshold,
      
      # Find start of locked periods
      locked_start = is_locked & !lag(is_locked, default = FALSE),
      locked_end = !is_locked & lag(is_locked, default = FALSE)
    )
  
  return(data)
}

flag_on_target <- function(focus_data, clipmarks) {
  if (is.null(clipmarks) || nrow(clipmarks) == 0) {
    focus_data$on_target <- FALSE
    return(focus_data)
  }
  
  cli_alert("Flagging on-target periods...")
  
  clipmarks$start_time <- as.POSIXct(clipmarks$start_time)
  clipmarks$stop_time <- as.POSIXct(clipmarks$stop_time)
  
  focus_data$on_target <- FALSE
  
  for (i in seq_len(nrow(clipmarks))) {
    clip <- clipmarks[i, ]
    in_window <- focus_data$mission_id == clip$mission_id &
                 focus_data$timestamp >= clip$start_time &
                 focus_data$timestamp <= clip$stop_time
    focus_data$on_target[in_window] <- TRUE
  }
  
  n_on_target <- sum(focus_data$on_target)
  cli_alert_success("Flagged {format(n_on_target, big.mark=',')} on-target samples")
  
  return(focus_data)
}

extract_event_summary <- function(focus_data, config) {
  # Create summary table of events for video extraction
  
  cli_alert("Extracting event summary for video processing...")
  
  # Autofocus events
  autofocus_events <- focus_data %>%
    filter(is_autofocus_start) %>%
    mutate(
      event_type = "autofocus",
      before_focus = lag(focus, n = config$sample_rate_hz),
      after_focus = lead(focus, n = config$sample_rate_hz * 2)
    ) %>%
    select(mission_id, camera, zoom_level, timestamp, event_type,
           before_focus, after_focus, on_target)
  
  # Manual adjustment events
  manual_events <- focus_data %>%
    filter(is_manual_adjust) %>%
    mutate(
      event_type = "manual",
      before_focus = lag(focus),
      after_focus = focus
    ) %>%
    select(mission_id, camera, zoom_level, timestamp, event_type,
           before_focus, after_focus, on_target)
  
  # Locked period starts (interesting for "in focus" examples)
  locked_events <- focus_data %>%
    filter(locked_start & on_target) %>%
    mutate(
      event_type = "locked_on_target",
      before_focus = focus,
      after_focus = focus
    ) %>%
    select(mission_id, camera, zoom_level, timestamp, event_type,
           before_focus, after_focus, on_target)
  
  all_events <- bind_rows(autofocus_events, manual_events, locked_events) %>%
    arrange(timestamp) %>%
    mutate(
      event_id = row_number(),
      time_hhmmss = format(timestamp, "%H%M%S"),
      date = as.Date(timestamp)
    )
  
  cli_alert_success("Extracted {nrow(all_events)} events")
  
  return(all_events)
}

# ------------------------------------------------------------------------------
# PLOTTING FUNCTIONS
# ------------------------------------------------------------------------------

plot_focus_timeline <- function(focus_data, title = "Focus Timeline", config) {
  cli_alert("Generating focus timeline plot...")
  
  # Downsample if needed
  if (nrow(focus_data) > 100000) {
    focus_data <- focus_data %>%
      group_by(camera, zoom_level) %>%
      sample_frac(0.1) %>%
      ungroup() %>%
      arrange(timestamp)
  }
  
  # Create zoom level background rectangles
  zoom_changes <- focus_data %>%
    group_by(camera, mission_id) %>%
    mutate(
      zoom_changed = zoom_level != lag(zoom_level),
      zoom_changed = ifelse(is.na(zoom_changed), TRUE, zoom_changed),
      zoom_group = cumsum(zoom_changed)
    ) %>%
    group_by(camera, mission_id, zoom_group, zoom_level) %>%
    summarise(
      xmin = min(timestamp),
      xmax = max(timestamp),
      .groups = "drop"
    )
  
  p <- ggplot() +
    # Zoom level backgrounds
    geom_rect(
      data = zoom_changes,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = zoom_level),
      alpha = 0.2
    ) +
    # Focus line
    geom_line(
      data = focus_data,
      aes(x = timestamp, y = focus),
      color = "black",
      linewidth = 0.3,
      alpha = 0.7
    ) +
    # Autofocus markers
    geom_point(
      data = focus_data %>% filter(is_autofocus_start),
      aes(x = timestamp, y = focus),
      color = "red",
      size = 2,
      shape = 17
    ) +
    # Manual adjustment markers
    geom_point(
      data = focus_data %>% filter(is_manual_adjust),
      aes(x = timestamp, y = focus),
      color = "orange",
      size = 1.5,
      shape = 16
    ) +
    # Facet by camera
    facet_wrap(~ camera, ncol = 1, scales = "free_x") +
    scale_fill_manual(
      values = c("Wide (0%)" = "#a6cee3", "Mid (50%)" = "#b2df8a", "Tele (100%)" = "#fb9a99"),
      name = "Zoom Level"
    ) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
    labs(
      title = title,
      subtitle = "Red triangles = autofocus | Orange dots = manual adjustments | Background = zoom level",
      x = "Time",
      y = "Focus (%)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      strip.text = element_text(size = 12, face = "bold"),
      legend.position = "bottom"
    )
  
  return(p)
}

plot_focus_histogram <- function(focus_data, title = "Focus Value Distribution", config) {
  p <- ggplot(focus_data, aes(x = focus, fill = zoom_level)) +
    geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
    facet_grid(camera ~ zoom_level) +
    scale_fill_manual(
      values = c("Wide (0%)" = "#a6cee3", "Mid (50%)" = "#b2df8a", "Tele (100%)" = "#fb9a99")
    ) +
    labs(
      title = title,
      x = "Focus (%)",
      y = "Count"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "none"
    )
  
  return(p)
}

plot_focus_rate_histogram <- function(focus_data, title = "Focus Change Rate", config) {
  # Filter to reasonable rates
  focus_data_filtered <- focus_data %>%
    filter(abs(focus_rate) < 100)
  
  p <- ggplot(focus_data_filtered, aes(x = focus_rate, fill = zoom_level)) +
    geom_histogram(bins = 100, alpha = 0.7) +
    facet_grid(camera ~ zoom_level) +
    scale_fill_manual(
      values = c("Wide (0%)" = "#a6cee3", "Mid (50%)" = "#b2df8a", "Tele (100%)" = "#fb9a99")
    ) +
    labs(
      title = title,
      subtitle = "Higher rates may indicate autofocus or rapid manual changes",
      x = "Focus Rate (%/sec)",
      y = "Count"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "none"
    )
  
  return(p)
}

plot_event_summary <- function(events, title = "Focus Events Summary", config) {
  if (nrow(events) == 0) {
    cli_alert_warning("No events to plot")
    return(NULL)
  }
  
  # Events by type and camera
  summary_data <- events %>%
    group_by(camera, zoom_level, event_type) %>%
    summarise(count = n(), .groups = "drop")
  
  p <- ggplot(summary_data, aes(x = zoom_level, y = count, fill = event_type)) +
    geom_col(position = "dodge") +
    facet_wrap(~ camera) +
    scale_fill_manual(
      values = c("autofocus" = "#e41a1c", "manual" = "#ff7f00", "locked_on_target" = "#4daf4a"),
      name = "Event Type"
    ) +
    labs(
      title = title,
      x = "Zoom Level",
      y = "Event Count"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  return(p)
}

# ------------------------------------------------------------------------------
# MAIN ANALYSIS
# ------------------------------------------------------------------------------

run_focus_analysis <- function(config, col_map, clip_col_map) {
  cli_h1("Focus Hunting Detection - Phase 1")
  
  # Initialize database
  con <- init_database(config, col_map)
  if (is.null(con)) return(NULL)
  
  # Query data
  focus_data <- query_focus_data(con, config, col_map)
  clipmarks <- query_clipmarks(con, config, clip_col_map)
  
  dbDisconnect(con)
  
  if (nrow(focus_data) == 0) {
    cli_alert_danger("No focus data found")
    return(NULL)
  }
  
  # Reshape to long format
  focus_long <- reshape_focus_data(focus_data)
  
  # Detect focus events
  focus_analyzed <- detect_focus_events(focus_long, config)
  
  # Flag on-target periods
  focus_analyzed <- flag_on_target(focus_analyzed, clipmarks)
  
  # Extract event summary
  events <- extract_event_summary(focus_analyzed, config)
  
  # Create output directory
  dir.create(config$output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Save events CSV
  events_file <- file.path(config$output_dir, "focus_events.csv")
  write.csv(events, events_file, row.names = FALSE)
  cli_alert_success("Saved: {events_file}")
  
  # Generate plots
  results <- list(
    data = focus_analyzed,
    events = events
  )
  
  # 1. Focus timeline
  cli_h2("Generating Focus Timeline")
  p_timeline <- plot_focus_timeline(focus_analyzed, "Focus Timeline - All Data", config)
  
  timeline_file <- file.path(config$output_dir, paste0("focus_timeline.", config$output_format))
  ggsave(timeline_file, p_timeline, width = config$output_width, height = config$output_height, dpi = config$output_dpi)
  cli_alert_success("Saved: {timeline_file}")
  results$timeline <- p_timeline
  
  # 2. Focus histogram
  p_hist <- plot_focus_histogram(focus_analyzed, "Focus Value Distribution", config)
  
  hist_file <- file.path(config$output_dir, paste0("focus_histogram.", config$output_format))
  ggsave(hist_file, p_hist, width = config$output_width, height = config$output_height, dpi = config$output_dpi)
  cli_alert_success("Saved: {hist_file}")
  results$histogram <- p_hist
  
  # 3. Focus rate histogram
  p_rate <- plot_focus_rate_histogram(focus_analyzed, "Focus Change Rate Distribution", config)
  
  rate_file <- file.path(config$output_dir, paste0("focus_rate_histogram.", config$output_format))
  ggsave(rate_file, p_rate, width = config$output_width, height = config$output_height, dpi = config$output_dpi)
  cli_alert_success("Saved: {rate_file}")
  results$rate_histogram <- p_rate
  
  # 4. Event summary
  if (nrow(events) > 0) {
    p_events <- plot_event_summary(events, "Focus Events by Camera and Zoom", config)
    
    if (!is.null(p_events)) {
      events_plot_file <- file.path(config$output_dir, paste0("focus_events_summary.", config$output_format))
      ggsave(events_plot_file, p_events, width = config$output_width, height = config$output_height / 2, dpi = config$output_dpi)
      cli_alert_success("Saved: {events_plot_file}")
      results$events_plot <- p_events
    }
  }
  
  # Summary statistics
  cli_h2("Summary Statistics")
  
  for (cam in unique(focus_analyzed$camera)) {
    cli_alert_info("{cam}:")
    cam_events <- events %>% filter(camera == cam)
    cli_alert_info("  Autofocus events: {sum(cam_events$event_type == 'autofocus')}")
    cli_alert_info("  Manual adjustments: {sum(cam_events$event_type == 'manual')}")
    cli_alert_info("  Locked on-target: {sum(cam_events$event_type == 'locked_on_target')}")
  }
  
  cli_h1("Phase 1 Complete")
  cli_alert_info("Run Phase 2 (focus_extract_frames.R) after placing videos in the video folder")
  
  return(results)
}

# ==============================================================================
# RUN
# ==============================================================================

# Uncomment to run:
# results <- run_focus_analysis(config, col_map, clip_col_map)
