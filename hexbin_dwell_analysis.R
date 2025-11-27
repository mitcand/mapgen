# ==============================================================================
# Hex Bin Dwell Time Heatmap Analysis
# ==============================================================================
# 
# PURPOSE:
#   Analyzes camera frame center positions to visualize dwell time patterns
#   from high-altitude aircraft scientific missions. Generates hex bin heatmaps
#   showing where the camera spent the most time, highlighting areas of crew
#   interest.
#
# AUTHOR: Generated with Claude AI assistance
# VERSION: 1.0
# LAST UPDATED: 2025-11-26
#
# USAGE:
#   1. Update the CONFIG section with your parameters
#   2. Update the COLUMN MAPPING section with your CSV column names
#   3. Run: result <- run_analysis(config, col_map)
#   4. View: print(result$plot)
#
# INPUT FILES:
#   - <mission_id>_vuefast.csv  : 15 Hz sensor data with frame center coordinates
#   - <mission_id>_clipmarks.csv: Crew-marked points of interest (optional)
#   - <mission_id>_vueslow.csv  : 1 Hz housekeeping data (not used by this script)
#
# OUTPUT:
#   - PNG/PDF/JPG heatmap image saved to output_dir
#   - Returns list with $plot (ggplot object) and $output_file (path)
#
# DEPENDENCIES:
#   DBI, duckdb, sf, s2, h3r, ggplot2, viridis, rnaturalearth, 
#   rnaturalearthdata, dplyr, lubridate, glue, cli
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
  
  # === MAP GENERATION ===
  generate_aircraft_map = TRUE,         # TRUE = generate both maps, FALSE = frame center only
  
  # === GEOGRAPHIC EXTENT ===
  # extent_mode: How to determine map bounds
  #   "data_countries"     - Countries containing frame center data (default)
  #   "data_bbox"          - Raw bounding box of frame center data
  #   "clipmark_countries" - Countries containing clipmark locations
  #   "clipmark_bbox"      - Raw bounding box of clipmark locations
  #   "country"            - Specific country (set country_filter below)
  #   "manual"             - Manual center point + extent (set manual_center, manual_extent_deg)
  extent_mode = "data_countries",
  
  # For extent_mode = "country":
  #   ISO3 codes: "AZE" (Azerbaijan), "TUR" (Turkey), "USA", "GBR", etc.
  #   Find codes at: https://en.wikipedia.org/wiki/ISO_3166-1_alpha-3
  country_filter = NULL,
  

  # For extent_mode = "manual":
  #   manual_center = c(lat, lon)  - Center of the map
  #   manual_extent_deg = 50       - Width in degrees (height auto from aspect ratio)
  manual_center = NULL,
  manual_extent_deg = NULL,
  
  # Buffer around extent (applies to all modes)
  extent_buffer_km = 50,
  
  # === HEX GRID ===
  hex_diameter_km = 25,                 # Target hex diameter in km (H3 picks closest size)
                                        # Smaller = more detail, larger = faster processing
                                        # Reference: 8km (fine), 25km (medium), 60km (coarse)
  
  # === DWELL TIME CLASSIFICATION ===
  # Defines the buckets for categorizing dwell time. Creates N+1 categories.
  # Adjust based on your analysis timeframe:
  #   - Single mission (1-2 days): c(0.5, 2, 5, 10)  
  #   - 30-90 days: c(5, 15, 30, 60)
  #   - 90+ days: c(30, 60, 90, 120)
  time_breaks_minutes = c(0.5, 2, 5, 10),
  
  # === COLORS ===
  # Dwell time colors (low to high). Must have length(time_breaks_minutes) + 1 colors.
  dwell_colors = c("#285DAB", "#5EA5DA", "#F1DCAA", "#F9A63F", "#CD5821"),
  
  # Basemap styling
  water_color = "#1a3a5c",              # Ocean/water background
  land_color = "#343434",               # Land fill color
  border_color = "#555555",             # Country border color
  border_width = 0.3,                   # Country border line width
  
  # Clipmark highlighting (hexes containing crew-marked POIs)
  clipmark_border_color = "cyan",       # Border color for highlighted hexes
  clipmark_border_width = 0.8,          # Border width for highlighted hexes
  clipmark_border_alpha = 1.0,          # Border transparency (0-1)
  
  # Hex bin styling
  hex_alpha = 0.85,                     # Transparency of hex fills (0-1)
  hex_border_color = "white",           # Border color for hex outlines
  hex_border_width = 0.1,               # Border width for hex outlines
  hex_border_alpha = 0.3,               # Transparency of hex borders (0-1)
  
  # === PROJECTION ===
  # projection_override: Force a specific projection, or NULL for auto-select
  #   Options: "laea" (Lambert Azimuthal), "albers", "robinson", "eqc" (Equirectangular),
  #            "winkel3", "mercator", NULL (auto based on extent)
  projection_override = NULL,
  
  # === OUTPUT ===
  output_dir = "output",                # Directory for saved images
  output_format = "png",                # Format: "png", "pdf", or "jpg"
  output_width = 16,                    # Image width in inches (map expands to fill)
  output_height = 9,                    # Image height in inches (16:9 default)
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
#          to:     timestamp = "gps_time"
# ------------------------------------------------------------------------------

col_map <- list(
  
  # --- VUEFAST COLUMNS (15 Hz sensor data) ---
  timestamp = "timestamp",              # UPDATE: Timestamp column (e.g., "gps_time", "utc")
  frame_lat = "frame_center_lat",       # UPDATE: Camera frame center latitude
  frame_lon = "frame_center_lon",       # UPDATE: Camera frame center longitude
  aircraft_lat = "aircraft_lat",        # UPDATE: Aircraft position latitude (optional)
  aircraft_lon = "aircraft_lon",        # UPDATE: Aircraft position longitude (optional)
  
  # NOTE: If your timestamp is split across date and time columns, you'll need
  # to modify the query_frame_centers() function to combine them.
  
  # --- CLIPMARKS COLUMNS (crew-marked points of interest) ---
  clip_start_time = "start_time",       # UPDATE: Clip start timestamp
  clip_end_time = "end_time",           # UPDATE: Clip end timestamp
  clip_lat = "poi_lat",                 # UPDATE: Point of interest latitude
  clip_lon = "poi_lon",                 # UPDATE: Point of interest longitude
  clip_description = "description"      # UPDATE: Crew description text
)

# ------------------------------------------------------------------------------
# PACKAGES
# ------------------------------------------------------------------------------

required_packages <- c(
  "DBI", "duckdb",                      # Database

  "sf", "s2",                           # Spatial operations
  "h3r",                                # Fast H3 hexagonal indexing
  "ggplot2", "viridis",                 # Plotting
  "rnaturalearth", "rnaturalearthdata", # Basemaps
  "dplyr", "lubridate", "glue",         # Data manipulation
  "cli"                                 # Progress/messages
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# Also need rnaturalearthhires for country boundaries
if (!requireNamespace("rnaturalearthhires", quietly = TRUE)) {
  install.packages("rnaturalearthhires", repos = "https://ropensci.r-universe.dev")
}

# ------------------------------------------------------------------------------
# SYNTHETIC DATA GENERATOR (for development without real data)
# ------------------------------------------------------------------------------

# KML ROUTE-BASED GENERATOR
# Place .kml files with LineString elements in a folder, and this will generate
# synthetic mission data that follows those routes.

parse_kml_route <- function(kml_file) {
  #' Parse a KML file and extract LineString coordinates

  #' 

  #' @param kml_file Path to KML file

  #' @return Data frame with lon, lat columns (in that order, as KML stores them)
  

  if (!file.exists(kml_file)) {
    cli_alert_danger("KML file not found: {kml_file}")
    return(NULL)
  }
  
  # Read KML as text
  kml_text <- paste(readLines(kml_file, warn = FALSE), collapse = "\n")
  
  # Extract coordinates from LineString elements
  # KML format: <coordinates>lon,lat,alt lon,lat,alt ...</coordinates>
  coord_pattern <- "<coordinates>\\s*([^<]+)\\s*</coordinates>"
  coord_matches <- regmatches(kml_text, gregexpr(coord_pattern, kml_text, perl = TRUE))[[1]]
  
  if (length(coord_matches) == 0) {
    cli_alert_warning("No LineString coordinates found in {basename(kml_file)}")
    return(NULL)
  }
  
  # Parse all coordinate strings
  all_coords <- data.frame(lon = numeric(), lat = numeric())
  
  for (match in coord_matches) {
    # Extract just the coordinate text
    coord_text <- gsub("</?coordinates>", "", match)
    coord_text <- trimws(coord_text)
    
    # Split by whitespace (each point is lon,lat,alt)
    points <- strsplit(coord_text, "\\s+")[[1]]
    points <- points[nchar(points) > 0]
    
    for (point in points) {
      parts <- as.numeric(strsplit(point, ",")[[1]])
      if (length(parts) >= 2 && !any(is.na(parts[1:2]))) {
        all_coords <- rbind(all_coords, data.frame(lon = parts[1], lat = parts[2]))
      }
    }
  }
  
  if (nrow(all_coords) == 0) {
    cli_alert_warning("Could not parse coordinates from {basename(kml_file)}")
    return(NULL)
  }
  
  cli_alert_success("Parsed {nrow(all_coords)} waypoints from {basename(kml_file)}")
  return(all_coords)
}

interpolate_route <- function(waypoints, target_samples, speed_variation = 0.1) {
  #' Interpolate waypoints to create smooth flight path at desired sample rate
  #' 
  #' @param waypoints Data frame with lon, lat columns
  #' @param target_samples Approximate number of output samples
  #' @param speed_variation Random speed variation factor (0-1)
  #' @return Data frame with interpolated lon, lat
  
  n_waypoints <- nrow(waypoints)
  if (n_waypoints < 2) return(waypoints)
  
  # Calculate segment distances (approximate degrees)
  segment_dists <- sqrt(
    diff(waypoints$lon)^2 + diff(waypoints$lat)^2
  )
  total_dist <- sum(segment_dists)
  
  # Allocate samples proportional to segment length
  samples_per_segment <- pmax(1, round((segment_dists / total_dist) * target_samples))
  
  # Interpolate each segment
  result <- data.frame(lon = numeric(), lat = numeric())
  
  for (i in 1:(n_waypoints - 1)) {
    n <- samples_per_segment[i]
    t <- seq(0, 1, length.out = n + 1)[1:n]  # Don't include endpoint (next segment starts there)
    
    # Add some speed variation (slow down/speed up randomly)
    if (speed_variation > 0) {
      t <- t + cumsum(rnorm(n, 0, speed_variation / n))
      t <- pmax(0, pmin(1, t))  # Clamp to [0,1]
      t <- sort(t)
    }
    
    seg_lon <- waypoints$lon[i] + t * (waypoints$lon[i+1] - waypoints$lon[i])
    seg_lat <- waypoints$lat[i] + t * (waypoints$lat[i+1] - waypoints$lat[i])
    
    result <- rbind(result, data.frame(lon = seg_lon, lat = seg_lat))
  }
  
  # Add final waypoint
  result <- rbind(result, waypoints[n_waypoints, ])
  
  return(result)
}

generate_camera_from_aircraft <- function(aircraft_lon, aircraft_lat, 
                                           look_mode = "mixed",
                                           offset_range_deg = c(0.1, 0.8)) {

  #' Generate camera frame center positions based on aircraft position
  #' 
  #' @param aircraft_lon Vector of aircraft longitudes
  #' @param aircraft_lat Vector of aircraft latitudes
  #' @param look_mode One of: "mixed", "side_look", "nadir", "forward"
  #' @param offset_range_deg Range of offset distance in degrees
  #' @return Data frame with frame_lon, frame_lat
  
  n <- length(aircraft_lon)
  frame_lon <- numeric(n)
  frame_lat <- numeric(n)
  
  # Calculate heading at each point
  heading <- numeric(n)
  for (i in 2:n) {
    heading[i] <- atan2(
      aircraft_lon[i] - aircraft_lon[i-1],
      aircraft_lat[i] - aircraft_lat[i-1]
    )
  }
  heading[1] <- heading[2]
  
  # Smooth heading
  for (i in 2:(n-1)) {
    heading[i] <- mean(heading[max(1,i-5):min(n,i+5)])
  }
  
  for (i in 1:n) {
    mode <- if (look_mode == "mixed") {
      sample(c("side_look", "nadir", "forward"), 1, prob = c(0.6, 0.25, 0.15))
    } else {
      look_mode
    }
    
    offset_dist <- runif(1, offset_range_deg[1], offset_range_deg[2])
    
    if (mode == "nadir") {
      # Looking straight down
      frame_lon[i] <- aircraft_lon[i] + rnorm(1, 0, 0.01)
      frame_lat[i] <- aircraft_lat[i] + rnorm(1, 0, 0.01)
      
    } else if (mode == "side_look") {
      # Looking perpendicular to track
      side <- sample(c(-1, 1), 1)
      offset_angle <- heading[i] + side * pi/2 + rnorm(1, 0, 0.1)
      frame_lon[i] <- aircraft_lon[i] + offset_dist * sin(offset_angle)
      frame_lat[i] <- aircraft_lat[i] + offset_dist * cos(offset_angle)
      
    } else if (mode == "forward") {
      # Looking ahead
      offset_angle <- heading[i] + rnorm(1, 0, 0.2)
      frame_lon[i] <- aircraft_lon[i] + offset_dist * 0.5 * sin(offset_angle)
      frame_lat[i] <- aircraft_lat[i] + offset_dist * 0.5 * cos(offset_angle)
    }
  }
  
  return(data.frame(frame_lon = frame_lon, frame_lat = frame_lat))
}

generate_from_kml <- function(kml_dir = "data/routes",
                               output_dir = "data/missions",
                               sample_rate_hz = 15,
                               flight_speed_kts = 300,
                               seed = 42) {
  #' Generate synthetic mission data from KML route files
  #' 
  #' @param kml_dir Directory containing .kml files with LineString routes
  #' @param output_dir Output directory for generated CSV files
  #' @param sample_rate_hz Sample rate (default 15 Hz)
  #' @param flight_speed_kts Approximate flight speed in knots
  #' @param seed Random seed for reproducibility
  #' @return List with generated mission info
  
  set.seed(seed)
  
  # Find KML files
  kml_files <- list.files(kml_dir, pattern = "\\.kml$", full.names = TRUE, 
                          ignore.case = TRUE)
  
  if (length(kml_files) == 0) {
    cli_alert_danger("No KML files found in {kml_dir}")
    cli_alert_info("Create a folder with .kml files containing LineString routes")
    cli_alert_info("You can create these in Google Earth by drawing paths and saving as KML")
    return(NULL)
  }
  
  cli_alert_info("Found {length(kml_files)} KML route file(s)")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Speed in degrees per second (approximate at mid-latitudes)
  # 1 knot ≈ 1.852 km/h, 1 degree ≈ 111 km at equator
  speed_deg_per_sec <- (flight_speed_kts * 1.852) / (111 * 3600)
  
  missions <- list()
  
  for (i in seq_along(kml_files)) {
    kml_file <- kml_files[i]
    mission_id <- sprintf("M2024%03d", i)
    route_name <- tools::file_path_sans_ext(basename(kml_file))
    
    cli_h2("Processing: {route_name}")
    
    # Parse route
    waypoints <- parse_kml_route(kml_file)
    if (is.null(waypoints)) next
    
    # Calculate route length and duration
    total_dist_deg <- sum(sqrt(diff(waypoints$lon)^2 + diff(waypoints$lat)^2))
    flight_duration_sec <- total_dist_deg / speed_deg_per_sec
    n_samples <- as.integer(flight_duration_sec * sample_rate_hz)
    
    cli_alert_info("Route length: ~{round(total_dist_deg * 111)} km")
    cli_alert_info("Flight duration: ~{round(flight_duration_sec / 3600, 1)} hours")
    cli_alert_info("Generating {format(n_samples, big.mark=',')} samples at {sample_rate_hz} Hz")
    
    # Interpolate route to sample rate
    path <- interpolate_route(waypoints, n_samples, speed_variation = 0.05)
    n_samples <- nrow(path)  # May differ slightly
    
    # Generate timestamps
    timestamps <- seq(
      from = as.POSIXct(paste0("2024-", sprintf("%02d", (i %% 12) + 1), "-15 06:00:00"), tz = "UTC"),
      by = 1/sample_rate_hz,
      length.out = n_samples
    )
    
    # Generate camera positions
    camera <- generate_camera_from_aircraft(path$lon, path$lat, look_mode = "mixed")
    
    # Create vuefast dataframe
    vuefast <- data.frame(
      timestamp = timestamps,
      frame_center_lat = camera$frame_lat,
      frame_center_lon = camera$frame_lon,
      aircraft_lat = path$lat,
      aircraft_lon = path$lon,
      altitude_ft = 60000 + rnorm(n_samples, 0, 500),
      camera_azimuth = runif(n_samples, 0, 360),
      camera_elevation = runif(n_samples, -90, 0)
    )
    
    # Create vueslow (1 Hz subset)
    vueslow_idx <- seq(1, n_samples, by = sample_rate_hz)
    vueslow <- data.frame(
      timestamp = timestamps[vueslow_idx],
      aircraft_lat = path$lat[vueslow_idx],
      aircraft_lon = path$lon[vueslow_idx],
      outside_air_temp_c = -56 + rnorm(length(vueslow_idx), 0, 2),
      sensor_temp_c = 20 + rnorm(length(vueslow_idx), 0, 1)
    )
    
    # Create clipmarks at interesting points along route
    n_clips <- min(15, max(3, nrow(waypoints)))
    clip_indices <- round(seq(1, n_samples, length.out = n_clips + 2)[2:(n_clips + 1)])
    
    clipmarks <- data.frame(
      start_time = timestamps[clip_indices],
      end_time = timestamps[pmin(clip_indices + sample_rate_hz * 60, n_samples)],
      poi_lat = camera$frame_lat[clip_indices],
      poi_lon = camera$frame_lon[clip_indices],
      description = paste("Waypoint", 1:n_clips)
    )
    
    # Save files
    vuefast_file <- file.path(output_dir, paste0(mission_id, "_vuefast.csv"))
    vueslow_file <- file.path(output_dir, paste0(mission_id, "_vueslow.csv"))
    clips_file <- file.path(output_dir, paste0(mission_id, "_clipmarks.csv"))
    
    write.csv(vuefast, vuefast_file, row.names = FALSE)
    write.csv(vueslow, vueslow_file, row.names = FALSE)
    write.csv(clipmarks, clips_file, row.names = FALSE)
    
    cli_alert_success("Generated: {mission_id} from {route_name}")
    cli_alert_info("  {format(n_samples, big.mark=',')} vuefast records")
    cli_alert_info("  {length(vueslow_idx)} vueslow records")
    cli_alert_info("  {nrow(clipmarks)} clipmarks")
    
    missions[[mission_id]] <- list(
      route_name = route_name,
      kml_file = kml_file,
      n_samples = n_samples,
      duration_hours = flight_duration_sec / 3600,
      extent = list(
        lat_range = range(path$lat),
        lon_range = range(path$lon)
      )
    )
  }
  
  cli_h2("KML Route Generation Complete")
  cli_alert_success("Generated {length(missions)} missions in {output_dir}")
  
  return(missions)
}

# RANDOM MISSION GENERATOR (original)

generate_synthetic_data <- function(output_dir = "data/missions",
                                    n_missions = 5,
                                    base_location = c(40.4093, 49.8671), # Baku, Azerbaijan
                                    seed = 42) {
  set.seed(seed)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  cli_alert_info("Generating {n_missions} synthetic missions...")
  
  for (m in 1:n_missions) {
    mission_id <- sprintf("M2024%03d", m)
    cli_alert("Mission {mission_id}")
    
    # Mission duration: 8-20 hours
    mission_hours <- runif(1, 8, 20)
    mission_seconds <- mission_hours * 3600
    n_samples <- as.integer(mission_seconds * 15)  # 15 Hz
    
    # Generate flight path with orbits and transits
    # Start near base
    ac_lat <- base_location[1]
    ac_lon <- base_location[2]
    
    timestamps <- seq(
      from = as.POSIXct(paste0("2024-", sprintf("%02d", (m %% 12) + 1), "-15 06:00:00"), tz = "UTC"),
      by = 1/15,
      length.out = n_samples
    )
    
    # Build path with segments: transit, orbit, transit, orbit, return
    # Initialize with base location (not zeros) to avoid 0,0 coordinates
    path_lat <- rep(base_location[1], n_samples)
    path_lon <- rep(base_location[2], n_samples)
    frame_lat <- rep(base_location[1], n_samples)
    frame_lon <- rep(base_location[2], n_samples)
    
    # Aircraft speed: ~300 knots = 0.00154 deg/sec at equator
    speed_deg_per_sample <- 0.00154 / 15
    
    # Define waypoints (random destinations within range)
    n_waypoints <- sample(3:6, 1)
    waypoints <- data.frame(
      lat = base_location[1] + runif(n_waypoints, -15, 15),
      lon = base_location[2] + runif(n_waypoints, -20, 20),
      orbit_time = sample(c(0, runif(n_waypoints - 1, 1800, 7200)), n_waypoints)  # 0 = transit only
    )
    # Add return to base
    waypoints <- rbind(waypoints, data.frame(lat = base_location[1], lon = base_location[2], orbit_time = 0))
    
    current_lat <- base_location[1]
    current_lon <- base_location[2]
    current_heading <- runif(1, 0, 2 * pi)
    idx <- 1
    wp_idx <- 1
    
    while (idx <= n_samples && wp_idx <= nrow(waypoints)) {
      target_lat <- waypoints$lat[wp_idx]
      target_lon <- waypoints$lon[wp_idx]
      orbit_samples <- as.integer(waypoints$orbit_time[wp_idx] * 15)
      
      # Transit to waypoint
      while (idx <= n_samples) {
        dist_to_target <- sqrt((target_lat - current_lat)^2 + (target_lon - current_lon)^2)
        
        if (dist_to_target < 0.05) break  # Close enough, start orbit or next waypoint
        
        # Head toward target with some wandering
        target_heading <- atan2(target_lon - current_lon, target_lat - current_lat)
        current_heading <- current_heading + 0.3 * (target_heading - current_heading) + rnorm(1, 0, 0.02)
        
        current_lat <- current_lat + speed_deg_per_sample * cos(current_heading)
        current_lon <- current_lon + speed_deg_per_sample * sin(current_heading)
        
        path_lat[idx] <- current_lat
        path_lon[idx] <- current_lon
        
        # Frame center: usually offset from aircraft, sometimes tracking a point
        if (runif(1) < 0.7) {
          # Looking off to the side
          offset_angle <- current_heading + sample(c(-1, 1), 1) * pi/2 + rnorm(1, 0, 0.1)
          offset_dist <- runif(1, 0.1, 0.8)  # 0.1 to 0.8 degrees ~ 10-80 km
          frame_lat[idx] <- current_lat + offset_dist * cos(offset_angle)
          frame_lon[idx] <- current_lon + offset_dist * sin(offset_angle)
        } else {
          # Looking down
          frame_lat[idx] <- current_lat + rnorm(1, 0, 0.01)
          frame_lon[idx] <- current_lon + rnorm(1, 0, 0.01)
        }
        
        idx <- idx + 1
      }
      
      # Orbit waypoint if specified
      if (orbit_samples > 0 && idx <= n_samples) {
        orbit_center_lat <- target_lat + rnorm(1, 0, 0.05)
        orbit_center_lon <- target_lon + rnorm(1, 0, 0.05)
        orbit_radius <- runif(1, 0.1, 0.3)  # degrees
        orbit_start_angle <- runif(1, 0, 2 * pi)
        
        for (o in 1:min(orbit_samples, n_samples - idx + 1)) {
          angle <- orbit_start_angle + (o / orbit_samples) * 2 * pi * sample(1:3, 1)  # 1-3 orbits
          
          current_lat <- orbit_center_lat + orbit_radius * cos(angle)
          current_lon <- orbit_center_lon + orbit_radius * sin(angle)
          
          path_lat[idx] <- current_lat
          path_lon[idx] <- current_lon
          
          # During orbit, camera often locks to center
          if (runif(1) < 0.8) {
            frame_lat[idx] <- orbit_center_lat + rnorm(1, 0, 0.02)
            frame_lon[idx] <- orbit_center_lon + rnorm(1, 0, 0.02)
          } else {
            frame_lat[idx] <- current_lat + rnorm(1, 0, 0.05)
            frame_lon[idx] <- current_lon + rnorm(1, 0, 0.05)
          }
          
          idx <- idx + 1
        }
      }
      
      wp_idx <- wp_idx + 1
    }
    
    # Fill any remaining samples (shouldn't happen often)
    if (idx <= n_samples) {
      path_lat[idx:n_samples] <- current_lat
      path_lon[idx:n_samples] <- current_lon
      frame_lat[idx:n_samples] <- current_lat
      frame_lon[idx:n_samples] <- current_lon
    }
    
    # Create vuefast dataframe
    vuefast <- data.frame(
      timestamp = timestamps,
      frame_center_lat = frame_lat,
      frame_center_lon = frame_lon,
      aircraft_lat = path_lat,
      aircraft_lon = path_lon,
      altitude_ft = 60000 + rnorm(n_samples, 0, 500),
      camera_azimuth = runif(n_samples, 0, 360),
      camera_elevation = runif(n_samples, -90, 0)
    )
    
    # Create vueslow (1 Hz subset with additional columns)
    vueslow_idx <- seq(1, n_samples, by = 15)
    vueslow <- data.frame(
      timestamp = timestamps[vueslow_idx],
      aircraft_lat = path_lat[vueslow_idx],
      aircraft_lon = path_lon[vueslow_idx],
      outside_air_temp_c = -56 + rnorm(length(vueslow_idx), 0, 2),
      sensor_temp_c = 20 + rnorm(length(vueslow_idx), 0, 1),
      fault_code = sample(c(0, 0, 0, 0, 1, 2), length(vueslow_idx), replace = TRUE)
    )
    
    # Create clipmarks (crew-marked POIs)
    n_clips <- sample(5:15, 1)
    clip_starts <- sort(sample(1:(n_samples - 1000), n_clips))
    clip_durations <- sample(300:3000, n_clips, replace = TRUE)  # 20 sec to 3 min in samples
    
    clipmarks <- data.frame(
      start_time = timestamps[clip_starts],
      end_time = timestamps[pmin(clip_starts + clip_durations, n_samples)],
      poi_lat = frame_lat[clip_starts] + rnorm(n_clips, 0, 0.01),
      poi_lon = frame_lon[clip_starts] + rnorm(n_clips, 0, 0.01),
      description = sample(c(
        "Industrial facility", "Port activity", "Urban area", 
        "Vehicle convoy", "Agricultural pattern", "Coastal feature",
        "Mountain terrain", "River crossing", "Airport",
        "Unknown structure"
      ), n_clips, replace = TRUE)
    )
    
    # Write CSVs
    write.csv(vuefast, file.path(output_dir, paste0(mission_id, "_vuefast.csv")), row.names = FALSE)
    write.csv(vueslow, file.path(output_dir, paste0(mission_id, "_vueslow.csv")), row.names = FALSE)
    write.csv(clipmarks, file.path(output_dir, paste0(mission_id, "_clipmarks.csv")), row.names = FALSE)
  }
  
  cli_alert_success("Generated {n_missions} missions in {output_dir}")
}

# ------------------------------------------------------------------------------
# DATABASE FUNCTIONS
# ------------------------------------------------------------------------------

init_database <- function(config, col_map) {
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
    cli_alert_info("Run generate_synthetic_data() to create test data")
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
      # Create table from first file
      query <- glue("
        CREATE TABLE vuefast AS 
        SELECT '{mission_id}' as mission_id, * 
        FROM read_csv_auto('{file_path}')
      ")
    } else {
      # Append to existing table
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
  
  # Create indexes for faster queries
  cli_alert("Creating indexes...")
  dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_vuefast_mission ON vuefast(mission_id)")
  dbExecute(con, glue("CREATE INDEX IF NOT EXISTS idx_vuefast_time ON vuefast({col_map$timestamp})"))
  
  # Report stats
  n_rows <- dbGetQuery(con, "SELECT COUNT(*) as n FROM vuefast")$n
  cli_alert_success("Database ready: {format(n_rows, big.mark=',')} vuefast records")
  
  return(con)
}

# ------------------------------------------------------------------------------
# DATA QUERY FUNCTIONS
# ------------------------------------------------------------------------------

query_clipmarks <- function(con, config, col_map) {
  # Check if clipmarks table exists
  if (!dbExistsTable(con, "clipmarks")) {
    cli_alert_warning("No clipmarks table found")
    return(NULL)
  }
  
  # Build WHERE clause based on config
  where_clauses <- c()
  
  if (!is.null(config$mission_ids)) {
    missions_str <- paste0("'", config$mission_ids, "'", collapse = ", ")
    where_clauses <- c(where_clauses, glue("mission_id IN ({missions_str})"))
  }
  
  if (!is.null(config$date_start)) {
    where_clauses <- c(where_clauses, glue("{col_map$clip_start_time} >= '{config$date_start}'"))
  }
  
  if (!is.null(config$date_end)) {
    where_clauses <- c(where_clauses, glue("{col_map$clip_start_time} <= '{config$date_end}'"))
  }
  
  where_sql <- if (length(where_clauses) > 0) {
    paste("WHERE", paste(where_clauses, collapse = " AND "))
  } else {
    ""
  }
  
  query <- glue("
    SELECT 
      {col_map$clip_lat} as lat,
      {col_map$clip_lon} as lon,
      {col_map$clip_description} as description
    FROM clipmarks
    {where_sql}
  ")
  
  cli_alert_info("Querying clipmarks...")
  data <- dbGetQuery(con, query)
  
  if (nrow(data) == 0) {
    cli_alert_warning("No clipmarks found for query parameters")
    return(NULL)
  }
  
  cli_alert_success("Retrieved {nrow(data)} clipmarks")
  
  # Remove invalid coordinates
  data <- data[!is.na(data$lat) & !is.na(data$lon) & 
               data$lat >= -90 & data$lat <= 90 &
               data$lon >= -180 & data$lon <= 180, ]
  
  return(data)
}

query_frame_centers <- function(con, config, col_map) {
  # Build WHERE clause based on config
  where_clauses <- c()
  
  if (!is.null(config$mission_ids)) {
    missions_str <- paste0("'", config$mission_ids, "'", collapse = ", ")
    where_clauses <- c(where_clauses, glue("mission_id IN ({missions_str})"))
  }
  
  if (!is.null(config$date_start)) {
    where_clauses <- c(where_clauses, glue("{col_map$timestamp} >= '{config$date_start}'"))
  }
  
  if (!is.null(config$date_end)) {
    where_clauses <- c(where_clauses, glue("{col_map$timestamp} <= '{config$date_end}'"))
  }
  
  where_sql <- if (length(where_clauses) > 0) {
    paste("WHERE", paste(where_clauses, collapse = " AND "))
  } else {
    ""
  }
  
  query <- glue("
    SELECT 
      {col_map$frame_lat} as lat,
      {col_map$frame_lon} as lon
    FROM vuefast
    {where_sql}
  ")
  
  cli_alert_info("Querying frame centers...")
  data <- dbGetQuery(con, query)
  cli_alert_success("Retrieved {format(nrow(data), big.mark=',')} points")
  
  # Remove invalid coordinates
  data <- data[!is.na(data$lat) & !is.na(data$lon) & 
               data$lat >= -90 & data$lat <= 90 &
               data$lon >= -180 & data$lon <= 180, ]
  
  return(data)
}

query_aircraft_positions <- function(con, config, col_map) {
  #' Query aircraft position data (where the aircraft was flying)
  #' 
  #' @param con DuckDB connection

  #' @param config Configuration list
  #' @param col_map Column mapping list
  #' @return Data frame with lat/lon columns
  
  # Check if aircraft position columns exist
  table_info <- dbGetQuery(con, "DESCRIBE vuefast")
  available_cols <- table_info$column_name
  
  if (!col_map$aircraft_lat %in% available_cols || !col_map$aircraft_lon %in% available_cols) {
    cli_alert_warning("Aircraft position columns not found in data")
    return(NULL)
  }
  
  # Build WHERE clause based on config
  where_clauses <- c()
  
  if (!is.null(config$mission_ids)) {
    mission_list <- paste0("'", config$mission_ids, "'", collapse = ", ")
    where_clauses <- c(where_clauses, glue("mission_id IN ({mission_list})"))
  }
  
  if (!is.null(config$date_start)) {
    where_clauses <- c(where_clauses, glue("{col_map$timestamp} >= '{config$date_start}'"))
  }
  
  if (!is.null(config$date_end)) {
    where_clauses <- c(where_clauses, glue("{col_map$timestamp} <= '{config$date_end}'"))
  }
  
  where_sql <- if (length(where_clauses) > 0) {
    paste("WHERE", paste(where_clauses, collapse = " AND "))
  } else {
    ""
  }
  
  query <- glue("
    SELECT 
      {col_map$aircraft_lat} as lat,
      {col_map$aircraft_lon} as lon
    FROM vuefast
    {where_sql}
  ")
  
  cli_alert_info("Querying aircraft positions...")
  data <- dbGetQuery(con, query)
  cli_alert_success("Retrieved {format(nrow(data), big.mark=',')} points")
  
  # Remove invalid coordinates
  data <- data[!is.na(data$lat) & !is.na(data$lon) & 
               data$lat >= -90 & data$lat <= 90 &
               data$lon >= -180 & data$lon <= 180, ]
  
  return(data)
}

# ------------------------------------------------------------------------------
# PROJECTION FUNCTIONS
# ------------------------------------------------------------------------------
#
# These functions handle map projection selection. The auto-selection logic
# chooses an appropriate projection based on the geographic extent of the data.
#
# Projection selection rules:
#   - Small regional (<30° span): Lambert Azimuthal Equal Area
#   - Medium regional (30-60°): Albers Equal Area Conic (or Polar Stereographic)
#   - Continental (60-120°): Winkel Tripel
#   - Global (>120°): Robinson or Mercator
# ------------------------------------------------------------------------------

choose_projection <- function(points_sf, config) {
  # Wrapper for sf objects - extracts bbox and calls main function
  bbox <- st_bbox(points_sf)
  return(choose_projection_from_bbox(bbox, config))
}

choose_projection_from_bbox <- function(bbox, config) {
  
  # Extract bbox values as plain numbers
  xmin <- as.numeric(bbox["xmin"])
  xmax <- as.numeric(bbox["xmax"])
  ymin <- as.numeric(bbox["ymin"])
  ymax <- as.numeric(bbox["ymax"])
  
  center_lon <- (xmin + xmax) / 2
  center_lat <- (ymin + ymax) / 2
  
  if (!is.null(config$projection_override)) {
    proj_string <- switch(config$projection_override,
      "laea" = glue("+proj=laea +lat_0={center_lat} +lon_0={center_lon} +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"),
      "winkel3" = glue("+proj=wintri +lon_0={center_lon} +datum=WGS84 +units=m +no_defs"),
      "mercator" = "+proj=merc +datum=WGS84 +units=m +no_defs",
      "robinson" = glue("+proj=robin +lon_0={center_lon} +datum=WGS84 +units=m +no_defs"),
      "albers" = glue("+proj=aea +lat_1={center_lat - 10} +lat_2={center_lat + 10} +lat_0={center_lat} +lon_0={center_lon} +datum=WGS84 +units=m +no_defs"),
      "eqc" = glue("+proj=eqc +lat_0=0 +lon_0={center_lon} +datum=WGS84 +units=m +no_defs"),
      # Default fallback
      glue("+proj=laea +lat_0={center_lat} +lon_0={center_lon} +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
    )
    cli_alert_info("Using override projection: {config$projection_override}")
    return(proj_string)
  }
  
  # Adjust extent calculation based on extent_mode
  extent_mode <- config$extent_mode
  
  if (extent_mode == "manual" && !is.null(config$manual_center) && !is.null(config$manual_extent_deg)) {
    center_lat <- config$manual_center[1]
    center_lon <- config$manual_center[2]
    xmin <- center_lon - config$manual_extent_deg / 2
    xmax <- center_lon + config$manual_extent_deg / 2
    ymin <- center_lat - config$manual_extent_deg / 2
    ymax <- center_lat + config$manual_extent_deg / 2
    cli_alert_info("Using manual extent for projection: center ({center_lat}°N, {center_lon}°E)")
    
  } else if (extent_mode == "country" && !is.null(config$country_filter)) {
    countries <- ne_countries(scale = 50, returnclass = "sf")
    country <- countries[countries$iso_a3 == config$country_filter, ]
    
    if (nrow(country) > 0) {
      # Special handling for countries that span the dateline
      special_extents <- list(
        RUS = c(20, 41, 180, 82),
        USA = c(-125, 24, -66, 50),
        CAN = c(-141, 42, -52, 72),
        NZL = c(166, -47, 179, -34),
        FJI = c(177, -21, 180, -12)
      )
      
      if (config$country_filter %in% names(special_extents)) {
        ext <- special_extents[[config$country_filter]]
        xmin <- ext[1]
        ymin <- ext[2]
        xmax <- ext[3]
        ymax <- ext[4]
        cli_alert_info("Using continental extent for projection: {config$country_filter}")
      } else {
        bbox <- st_bbox(country)
        xmin <- as.numeric(bbox["xmin"])
        xmax <- as.numeric(bbox["xmax"])
        ymin <- as.numeric(bbox["ymin"])
        ymax <- as.numeric(bbox["ymax"])
        cli_alert_info("Using country extent for projection: {config$country_filter}")
      }
      center_lon <- (xmin + xmax) / 2
      center_lat <- (ymin + ymax) / 2
    }
  }
  # For data_countries, data_bbox, clipmark_* modes, use the passed bbox
  
  # Calculate extent characteristics
  lon_span <- xmax - xmin
  lat_span <- ymax - ymin
  
  # Adjust lon_span for wrapping if needed
  if (!is.na(lon_span) && lon_span > 180) {
    lon_span <- 360 - lon_span
  }
  
  cli_alert("Extent analysis: {round(lat_span, 1)}° lat × {round(lon_span, 1)}° lon, centered at {round(center_lat, 1)}°N, {round(center_lon, 1)}°E")
  
  # Decision logic based on extent
  # 
  # Small regional (< 30° span): Lambert Azimuthal Equal Area - best for small areas, equal area
  # Medium regional (30-60° span): Albers Equal Area Conic - good for mid-latitudes, equal area
  # Large/continental (60-120° span): Robinson - fills rectangular frame
  # Global/multi-continental (> 120° span): Robinson or Equirectangular
  #
  # Special cases:
  # - High latitudes (>65°) with SMALL extent: Use LAEA
  # - Large extents always use rectangular projections (Robinson/Equirectangular)
  
  max_span <- max(lon_span, lat_span)
  max_lat <- max(abs(ymin), abs(ymax))
  
  # Check for high-latitude data
  is_polar_adjacent <- max_lat > 65
  is_large_extent <- max_span > 60 || lat_span > 50 || lon_span > 80
  
  if (is_large_extent) {
    # Large extents - use Equirectangular (Plate Carrée) for simplicity
    # This maps lat/lon directly to x/y, avoiding complex distortion issues
    proj_string <- glue("+proj=eqc +lat_0=0 +lon_0={center_lon} +datum=WGS84 +units=m +no_defs")
    proj_name <- "Equirectangular (large extent)"
    
  } else if (is_polar_adjacent) {
    # Small extent near poles - LAEA handles this well
    proj_string <- glue("+proj=laea +lat_0={center_lat} +lon_0={center_lon} +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
    proj_name <- "Lambert Azimuthal Equal Area (high latitude)"
    
  } else if (max_span < 30) {
    # Small regional - Lambert Azimuthal Equal Area
    proj_string <- glue("+proj=laea +lat_0={center_lat} +lon_0={center_lon} +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
    proj_name <- "Lambert Azimuthal Equal Area (regional)"
    
  } else if (max_span < 60) {
    # Medium regional - Albers for mid-latitudes
    if (abs(center_lat) > 60) {
      # High latitude - use LAEA
      proj_string <- glue("+proj=laea +lat_0={center_lat} +lon_0={center_lon} +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
      proj_name <- "Lambert Azimuthal Equal Area (high latitude regional)"
    } else {
      # Mid-latitudes - Albers Equal Area
      lat1 <- center_lat - lat_span / 4
      lat2 <- center_lat + lat_span / 4
      proj_string <- glue("+proj=aea +lat_1={lat1} +lat_2={lat2} +lat_0={center_lat} +lon_0={center_lon} +datum=WGS84 +units=m +no_defs")
      proj_name <- "Albers Equal Area Conic"
    }
    
  } else {
    # Fallback - Robinson
    proj_string <- glue("+proj=robin +lon_0={center_lon} +datum=WGS84 +units=m +no_defs")
    proj_name <- "Robinson (continental)"
  }
  
  cli_alert_info("Auto-selected projection: {proj_name}")
  return(proj_string)
}

get_country_boundary <- function(country_filter) {
  if (is.null(country_filter)) return(NULL)
  
  countries <- ne_countries(scale = 50, returnclass = "sf")
  country <- countries[countries$iso_a3 == country_filter, ]
  
  if (nrow(country) == 0) {
    cli_alert_warning("Country code '{country_filter}' not found")
    return(NULL)
  }
  
  return(country)
}

# ------------------------------------------------------------------------------
# BINNING AND CLASSIFICATION (H3 hexagonal indexing via h3r)
# ------------------------------------------------------------------------------
#
# This section uses Uber's H3 hierarchical hexagonal grid system for fast
# point-to-hex assignment. H3 is MUCH faster than spatial joins because it
# computes hex membership mathematically (O(n)) rather than doing spatial
# predicate tests (O(n*m)).
#
# Performance comparison for 16M points:
#   - Spatial join (st_join): Hours
#   - H3 indexing: ~5-10 seconds
#
# H3 Resolution reference (approximate edge lengths):
#   Res 0: 1108 km  |  Res 5: 8.5 km
#   Res 1: 419 km   |  Res 6: 3.2 km
#   Res 2: 158 km   |  Res 7: 1.2 km
#   Res 3: 60 km    |  Res 8: 0.46 km
#   Res 4: 23 km    |  Res 9: 0.17 km
#
# More info: https://h3geo.org/docs/core-library/restable/
# ------------------------------------------------------------------------------

get_h3_resolution <- function(hex_diameter_km) {
  #' Select H3 resolution based on target hex diameter
  #'
  #' @param hex_diameter_km Target hexagon diameter in kilometers
 #' @return Integer H3 resolution (0-9)
  # H3 resolution edge lengths (approximate, in km)
  # Hex diameter is roughly 2x edge length
  h3_lookup <- data.frame(
    res = 0:9,
    edge_km = c(1107.71, 418.68, 158.24, 59.81, 22.61, 8.54, 3.23, 1.22, 0.46, 0.17)
  )
  
  target_edge <- hex_diameter_km / 2
  
  # Find closest resolution
  diffs <- abs(h3_lookup$edge_km - target_edge)
  best_idx <- which.min(diffs)
  best_res <- h3_lookup$res[best_idx]
  
  cli_alert_info("H3 resolution {best_res} selected (edge ~{round(h3_lookup$edge_km[best_idx], 1)} km) for target diameter {hex_diameter_km} km")
  
  return(as.integer(best_res))
}

bin_dwell_times_fast <- function(frame_data, config, sample_rate_hz = 15) {
  n_points <- nrow(frame_data)
  cli_alert("Binning {format(n_points, big.mark=',')} points using H3...")
  
  # Get H3 resolution
  h3_res <- get_h3_resolution(config$hex_diameter_km)
  
  # Compute H3 index for each point (vectorized, very fast)
  start_time <- Sys.time()
  
  h3_indices <- h3r::latLngToCell(
    lat = frame_data$lat,
    lng = frame_data$lon,
    res = h3_res
  )
  
  elapsed_index <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1)
  cli_alert_success("H3 indexing complete in {elapsed_index}s")
  
  # Count points per hex
  cli_alert("Counting points per hex...")
  hex_counts <- as.data.frame(table(h3_indices), stringsAsFactors = FALSE)
  names(hex_counts) <- c("h3_index", "point_count")
  hex_counts$point_count <- as.numeric(hex_counts$point_count)
  
  # Remove any invalid indices
  hex_counts <- hex_counts[nchar(hex_counts$h3_index) > 0 & !is.na(hex_counts$h3_index), ]
  
  cli_alert_success("Found {format(nrow(hex_counts), big.mark=',')} unique hexagons with data")
  
  # Convert counts to dwell time in minutes
  hex_counts$dwell_minutes <- hex_counts$point_count / sample_rate_hz / 60
  
  # Classify into time buckets
  time_breaks <- config$time_breaks_minutes
  break_labels <- c(
    paste0("<", time_breaks[1], " min"),
    paste0(time_breaks[-length(time_breaks)], "-", time_breaks[-1], " min"),
    paste0(">", time_breaks[length(time_breaks)], " min")
  )
  
  hex_counts$dwell_category <- cut(
    hex_counts$dwell_minutes,
    breaks = c(0, time_breaks, Inf),
    labels = break_labels,
    include.lowest = TRUE
  )
  
  # Generate hex polygons
  cli_alert("Generating {nrow(hex_counts)} hex polygons...")
  start_time <- Sys.time()
  
  # Get boundaries for all cells
  boundaries <- h3r::cellToBoundary(hex_counts$h3_index)
  
  # Convert to sf polygons
  hex_polygons <- lapply(boundaries, function(b) {
    # Close the polygon (first point = last point)
    coords <- rbind(as.matrix(b), as.matrix(b)[1, , drop = FALSE])
    st_polygon(list(coords[, c("lng", "lat")]))
  })
  
  elapsed_poly <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1)
  cli_alert_success("Polygon generation complete in {elapsed_poly}s")
  
  # Create sf object
  hex_data <- st_sf(
    hex_id = hex_counts$h3_index,
    h3_index = hex_counts$h3_index,
    point_count = hex_counts$point_count,
    dwell_minutes = hex_counts$dwell_minutes,
    dwell_category = hex_counts$dwell_category,
    geometry = st_sfc(hex_polygons, crs = 4326)
  )
  
  # Summary stats
  cat("\nDwell time distribution:\n")
  print(table(hex_data$dwell_category))
  cat("\n")
  
  total_elapsed <- round(elapsed_index + elapsed_poly, 1)
  cli_alert_success("Total binning time: {total_elapsed}s")
  
  return(hex_data)
}

flag_clipmark_hexes_fast <- function(hex_data, clipmark_data, config) {
  # Flag hexes that contain clipmarks using H3 indexing
  if (is.null(clipmark_data) || nrow(clipmark_data) == 0) {
    hex_data$has_clipmark <- FALSE
    return(hex_data)
  }
  
  cli_alert("Identifying hexes with crew clipmarks...")
  
  # Get H3 resolution (same as used for hex_data)
  h3_res <- get_h3_resolution(config$hex_diameter_km)
  
  # Get H3 indices for clipmarks
  clip_h3_indices <- h3r::latLngToCell(
    lat = clipmark_data$lat,
    lng = clipmark_data$lon,
    res = h3_res
  )
  
  # Flag hexes that match clipmark indices
  hex_data$has_clipmark <- hex_data$h3_index %in% unique(clip_h3_indices)
  
  n_flagged <- sum(hex_data$has_clipmark)
  cli_alert_success("Flagged {n_flagged} hexagons containing clipmarks")
  
  return(hex_data)
}

# ------------------------------------------------------------------------------
# PLOTTING
# ------------------------------------------------------------------------------

create_heatmap <- function(hex_data, projection, config, clipmarks_sf = NULL, title = NULL) {
  cli_alert("Creating heatmap...")
  
  # Separate hexes with clipmarks for border overlay
  hex_with_clips <- hex_data[hex_data$has_clipmark == TRUE, ]
  n_clipmark_hexes <- nrow(hex_with_clips)
  
  # Load basemap
  cli_alert("Loading basemap...")
  world <- ne_countries(scale = 50, returnclass = "sf")
  
  # Disable s2 for simpler planar operations
  sf_use_s2_setting <- sf::sf_use_s2()
  sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(sf_use_s2_setting), add = TRUE)
  
  # Determine extent based on extent_mode
  extent_mode <- config$extent_mode
  data_extent <- NULL
  
  cli_alert_info("Extent mode: {extent_mode}")
  
  if (extent_mode == "manual") {
    # Manual center + extent
    if (is.null(config$manual_center) || is.null(config$manual_extent_deg)) {
      cli_alert_warning("Manual mode requires manual_center and manual_extent_deg, falling back to data_countries")
      extent_mode <- "data_countries"
    } else {
      center_lat <- config$manual_center[1]
      center_lon <- config$manual_center[2]
      
      output_aspect <- config$output_width / config$output_height
      if (output_aspect >= 1) {
        half_lon <- config$manual_extent_deg / 2
        half_lat <- (config$manual_extent_deg / output_aspect) / 2
      } else {
        half_lat <- config$manual_extent_deg / 2
        half_lon <- (config$manual_extent_deg * output_aspect) / 2
      }
      
      data_extent <- st_bbox(c(
        xmin = center_lon - half_lon,
        ymin = center_lat - half_lat,
        xmax = center_lon + half_lon,
        ymax = center_lat + half_lat
      ), crs = st_crs(4326))
      cli_alert_info("Manual extent: center ({center_lat}°N, {center_lon}°E), span {config$manual_extent_deg}°")
    }
  }
  
  if (extent_mode == "country") {
    # Specific country filter
    if (is.null(config$country_filter)) {
      cli_alert_warning("Country mode requires country_filter, falling back to data_countries")
      extent_mode <- "data_countries"
    } else {
      filter_country <- world[world$iso_a3 == config$country_filter, ]
      if (nrow(filter_country) > 0) {
        cli_alert_info("Using country: {filter_country$name[1]} ({config$country_filter})")
        
        # Special handling for dateline-spanning countries
        special_extents <- list(
          RUS = c(20, 41, 180, 82),
          USA = c(-125, 24, -66, 50),
          CAN = c(-141, 42, -52, 72),
          NZL = c(166, -47, 179, -34),
          FJI = c(177, -21, 180, -12)
        )
        
        if (config$country_filter %in% names(special_extents)) {
          ext <- special_extents[[config$country_filter]]
          data_extent <- st_bbox(c(xmin = ext[1], ymin = ext[2], xmax = ext[3], ymax = ext[4]),
                                  crs = st_crs(4326))
          cli_alert_info("Using continental extent for {config$country_filter}")
        } else {
          data_extent <- st_bbox(filter_country)
        }
      } else {
        cli_alert_warning("Country '{config$country_filter}' not found, falling back to data_countries")
        extent_mode <- "data_countries"
      }
    }
  }
  
  if (extent_mode == "clipmark_bbox") {
    # Raw clipmark bounding box
    if (is.null(clipmarks_sf) || nrow(clipmarks_sf) == 0) {
      cli_alert_warning("No clipmarks available, falling back to data_countries")
      extent_mode <- "data_countries"
    } else {
      data_extent <- st_bbox(clipmarks_sf)
      cli_alert_info("Using clipmark bbox: {nrow(clipmarks_sf)} clipmarks")
    }
  }
  
  if (extent_mode == "clipmark_countries") {
    # Countries containing clipmarks
    if (is.null(clipmarks_sf) || nrow(clipmarks_sf) == 0) {
      cli_alert_warning("No clipmarks available, falling back to data_countries")
      extent_mode <- "data_countries"
    } else {
      clip_intersects <- suppressWarnings(st_intersects(world, clipmarks_sf, sparse = FALSE))
      clip_countries <- world[apply(clip_intersects, 1, any), ]
      if (nrow(clip_countries) > 0) {
        cli_alert_info("Clipmarks span {nrow(clip_countries)} countries: {paste(clip_countries$name[1:min(3, nrow(clip_countries))], collapse=', ')}{if(nrow(clip_countries) > 3) '...' else ''}")
        
        # Break into polygon parts to filter overseas territories
        all_polygons <- suppressWarnings(st_cast(clip_countries, "POLYGON"))
        polygon_intersects <- suppressWarnings(st_intersects(all_polygons, clipmarks_sf, sparse = FALSE))
        polygon_has_data <- apply(polygon_intersects, 1, any)
        clip_polygons <- all_polygons[polygon_has_data, ]
        
        data_extent <- st_bbox(clip_polygons)
      } else {
        cli_alert_warning("No countries found for clipmarks, falling back to data_countries")
        extent_mode <- "data_countries"
      }
    }
  }
  
  if (extent_mode == "data_bbox") {
    # Raw hex data bounding box
    data_extent <- st_bbox(hex_data)
    cli_alert_info("Using data bbox")
  }
  
  if (extent_mode == "data_countries" || is.null(data_extent)) {
    # Countries containing hex data (default)
    # Key improvement: break multipolygons into parts and only use parts that intersect data
    # This prevents France's French Guiana, UK's Falklands, etc. from stretching extent
    
    hex_centroids <- suppressWarnings(st_centroid(hex_data))
    
    # First, find which countries have any data
    countries_with_data <- suppressWarnings(st_intersects(world, hex_centroids, sparse = FALSE))
    countries_with_data <- apply(countries_with_data, 1, any)
    data_countries <- world[countries_with_data, ]
    
    if (nrow(data_countries) == 0) {
      cli_alert_warning("No countries found, using data bbox")
      data_extent <- st_bbox(hex_data)
    } else {
      cli_alert_info("Data spans {nrow(data_countries)} countries: {paste(data_countries$name[1:min(5, nrow(data_countries))], collapse=', ')}{if(nrow(data_countries) > 5) '...' else ''}")
      
      # Break multipolygons into individual polygons to handle overseas territories
      # Countries like France (French Guiana), UK (Falklands), Netherlands (Caribbean) have
      # distant territories that shouldn't affect extent when data is only in the main territory
      cli_alert_info("Filtering to polygon parts containing data (excludes overseas territories)...")
      
      all_polygons <- suppressWarnings(st_cast(data_countries, "POLYGON"))
      
      # Find which individual polygons actually contain data
      polygon_intersects <- suppressWarnings(st_intersects(all_polygons, hex_centroids, sparse = FALSE))
      polygon_has_data <- apply(polygon_intersects, 1, any)
      data_polygons <- all_polygons[polygon_has_data, ]
      
      n_total_polys <- nrow(all_polygons)
      n_data_polys <- nrow(data_polygons)
      if (n_data_polys < n_total_polys) {
        cli_alert_info("Using {n_data_polys} of {n_total_polys} polygon parts (filtered distant territories)")
      }
      
      # Now apply problematic country exclusion to remaining polygons
      problematic_iso3 <- c("RUS", "USA", "FJI", "NZL", "KIR")
      if (nrow(data_polygons) > 1) {
        extent_polygons <- data_polygons[!data_polygons$iso_a3 %in% problematic_iso3, ]
        if (nrow(extent_polygons) == 0) {
          data_extent <- st_bbox(hex_data)
        } else {
          excluded <- unique(data_polygons$name[data_polygons$iso_a3 %in% problematic_iso3])
          if (length(excluded) > 0) {
            cli_alert_info("Excluded from extent: {paste(excluded, collapse=', ')}")
          }
          data_extent <- st_bbox(extent_polygons)
        }
      } else {
        data_extent <- st_bbox(data_polygons)
      }
    }
  }
  
  # Add buffer to extent (in degrees since we're in WGS84)
  buffer_deg <- config$extent_buffer_km / 111  # ~111 km per degree
  
  xmin_wgs84 <- as.numeric(data_extent["xmin"]) - buffer_deg
  xmax_wgs84 <- as.numeric(data_extent["xmax"]) + buffer_deg
  ymin_wgs84 <- as.numeric(data_extent["ymin"]) - buffer_deg
  ymax_wgs84 <- as.numeric(data_extent["ymax"]) + buffer_deg
  
  # Clamp to valid ranges

  xmin_wgs84 <- max(xmin_wgs84, -180)
  xmax_wgs84 <- min(xmax_wgs84, 180)
  ymin_wgs84 <- max(ymin_wgs84, -85)
  ymax_wgs84 <- min(ymax_wgs84, 85)
  
  # Calculate aspect ratio expansion in WGS84
  output_aspect <- config$output_width / config$output_height
  
  lon_range <- xmax_wgs84 - xmin_wgs84
  lat_range <- ymax_wgs84 - ymin_wgs84
  center_lat <- (ymin_wgs84 + ymax_wgs84) / 2
  
  # For large extents (using Equirectangular), treat degrees as equal units
  # For smaller extents, account for latitude distortion
  is_large <- (lon_range > 80 || lat_range > 50)
  
  if (is_large) {
    # Equirectangular: degrees map directly, so simple ratio
    current_aspect <- lon_range / lat_range
    
    if (current_aspect < output_aspect) {
      target_lon_range <- lat_range * output_aspect
      expand_lon <- (target_lon_range - lon_range) / 2
      xmin_wgs84 <- max(xmin_wgs84 - expand_lon, -180)
      xmax_wgs84 <- min(xmax_wgs84 + expand_lon, 180)
    } else {
      target_lat_range <- lon_range / output_aspect
      expand_lat <- (target_lat_range - lat_range) / 2
      ymin_wgs84 <- max(ymin_wgs84 - expand_lat, -85)
      ymax_wgs84 <- min(ymax_wgs84 + expand_lat, 85)
    }
  } else {
    # Account for latitude distortion in aspect calculation
    effective_lon_range <- lon_range * cos(center_lat * pi / 180)
    current_aspect <- effective_lon_range / lat_range
    
    if (current_aspect < output_aspect) {
      target_lon_range <- lat_range * output_aspect / cos(center_lat * pi / 180)
      expand_lon <- (target_lon_range - lon_range) / 2
      xmin_wgs84 <- max(xmin_wgs84 - expand_lon, -180)
      xmax_wgs84 <- min(xmax_wgs84 + expand_lon, 180)
    } else {
      target_lat_range <- effective_lon_range / output_aspect
      expand_lat <- (target_lat_range - lat_range) / 2
      ymin_wgs84 <- max(ymin_wgs84 - expand_lat, -85)
      ymax_wgs84 <- min(ymax_wgs84 + expand_lat, 85)
    }
  }
  
  cli_alert_info("Map extent: {round(xmin_wgs84, 1)}° to {round(xmax_wgs84, 1)}°E, {round(ymin_wgs84, 1)}° to {round(ymax_wgs84, 1)}°N")
  
  # Filter basemap to countries in/near the extent for performance
  extent_box <- st_bbox(c(xmin = xmin_wgs84, ymin = ymin_wgs84, 
                          xmax = xmax_wgs84, ymax = ymax_wgs84), 
                        crs = st_crs(4326))
  extent_rect <- st_as_sfc(extent_box)
  visible_idx <- suppressWarnings(st_intersects(world, extent_rect, sparse = FALSE)[,1])
  world_visible <- world[visible_idx, ]
  
  # For large extents, clip countries that span the antimeridian to avoid artifacts
  # Instead of excluding them, we crop them to the visible extent
  lon_span <- xmax_wgs84 - xmin_wgs84
  if (lon_span > 100) {
    problematic_iso3 <- c("RUS", "FJI", "KIR", "NZL", "USA", "CAN", "NOR", "ATA")
    
    # Separate problematic and normal countries
    problem_countries <- world_visible[world_visible$iso_a3 %in% problematic_iso3, ]
    normal_countries <- world_visible[!world_visible$iso_a3 %in% problematic_iso3, ]
    
    if (nrow(problem_countries) > 0) {
      # Crop problematic countries to the visible extent
      crop_box <- st_bbox(c(xmin = xmin_wgs84, ymin = ymin_wgs84, 
                            xmax = xmax_wgs84, ymax = ymax_wgs84), 
                          crs = st_crs(4326))
      
      problem_cropped <- tryCatch({
        cropped <- suppressWarnings(st_crop(problem_countries, crop_box))
        # Remove any empty geometries after cropping
        cropped <- cropped[!st_is_empty(cropped), ]
        cropped
      }, error = function(e) {
        cli_alert_warning("Could not crop problematic countries: {e$message}")
        NULL
      })
      
      if (!is.null(problem_cropped) && nrow(problem_cropped) > 0) {
        cli_alert_info("Cropped {nrow(problem_cropped)} antimeridian-spanning countries to visible extent")
        world_visible <- suppressWarnings(rbind(normal_countries, problem_cropped))
      } else {
        cli_alert_info("Excluded {nrow(problem_countries)} problematic countries from basemap")
        world_visible <- normal_countries
      }
    }
  }
  
  # Generate title if not provided
  if (is.null(title)) {
    title <- "Camera Dwell Time Heatmap"
  }
  
  # Use configured colors
  colors <- config$dwell_colors
  
  # Validate color count matches categories
  n_categories <- length(levels(hex_data$dwell_category))
  if (length(colors) != n_categories) {
    cli_alert_warning("Color count ({length(colors)}) doesn't match categories ({n_categories}), using viridis fallback")
    colors <- viridis::plasma(n_categories, begin = 0.1, end = 0.9)
  }
  
  # Create plot - everything in WGS84, coord_sf handles projection
  p <- ggplot() +
    # Land polygons with country borders
    geom_sf(data = world_visible, 
            fill = config$land_color, 
            color = config$border_color, 
            linewidth = config$border_width) +
    
    # Hex bins
    geom_sf(data = hex_data, aes(fill = dwell_category), color = NA, alpha = config$hex_alpha) +
    
    # Hex outlines - subtle for definition
    geom_sf(data = hex_data, fill = NA, color = config$hex_border_color, 
            linewidth = config$hex_border_width, alpha = config$hex_border_alpha) +
    
    # Clipmark hex borders
    {if (n_clipmark_hexes > 0) 
      geom_sf(data = hex_with_clips, fill = NA, 
              color = config$clipmark_border_color, 
              linewidth = config$clipmark_border_width,
              alpha = config$clipmark_border_alpha)
    } +
    
    # Color scale
    scale_fill_manual(
      values = colors,
      name = "Dwell Time",
      drop = FALSE
    ) +
    
    # Projection with WGS84 limits
    # For large extents, disable graticule (datum=NA) to avoid artifacts
    coord_sf(
      xlim = c(xmin_wgs84, xmax_wgs84),
      ylim = c(ymin_wgs84, ymax_wgs84),
      expand = FALSE,
      crs = projection,
      default_crs = sf::st_crs(4326),
      datum = if (is_large) NA else sf::st_crs(4326)
    ) +
    
    # Labels
    labs(
      title = title,
      subtitle = paste0(
        format(nrow(hex_data), big.mark = ","), " hexagons with coverage",
        if (n_clipmark_hexes > 0) paste0(" | ", n_clipmark_hexes, " with crew clipmarks (", config$clipmark_border_color, " border)") else ""
      )
    ) +
    
    # Theme - minimal margins to fill frame
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = config$water_color, color = NA),
      panel.grid = element_line(color = "#2a4a6a", linewidth = 0.2),
      legend.position = c(0.98, 0.5),  # Legend inside plot on right
      legend.justification = c(1, 0.5),
      legend.background = element_rect(fill = alpha("white", 0.9), color = "gray80"),
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10, face = "bold"),
      legend.key.size = unit(0.8, "lines"),
      plot.title = element_text(size = 14, face = "bold", margin = margin(0, 0, 2, 0)),
      plot.subtitle = element_text(size = 10, color = "gray40", margin = margin(0, 0, 5, 0)),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(5, 5, 5, 5),
      axis.text = element_text(size = 8)
    )
  
  return(p)
}

# ------------------------------------------------------------------------------
# FILENAME GENERATION
# ------------------------------------------------------------------------------

generate_output_filename <- function(config) {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  if (!is.null(config$mission_ids)) {
    # Mission-based filename
    mission_str <- paste(config$mission_ids, collapse = "_")
    filename <- paste0(mission_str, "_", timestamp)
    
  } else if (!is.null(config$date_start) || !is.null(config$date_end)) {
    # Date range filename - calculate days covered
    start_date <- if (!is.null(config$date_start)) as.Date(config$date_start) else as.Date("1970-01-01")
    end_date <- if (!is.null(config$date_end)) as.Date(config$date_end) else Sys.Date()
    
    n_days <- as.integer(end_date - start_date)
    filename <- paste0(timestamp, "_", n_days, "days")
    
  } else {
    # All data - just timestamp
    filename <- paste0("all_missions_", timestamp)
  }
  
  # Add country filter if specified
  if (!is.null(config$country_filter)) {
    filename <- paste0(filename, "_", config$country_filter)
  }
  
  # Create output directory if needed
  if (!dir.exists(config$output_dir)) {
    dir.create(config$output_dir, recursive = TRUE)
  }
  
  # Full path with extension
  full_path <- file.path(config$output_dir, paste0(filename, ".", config$output_format))
  
  return(full_path)
}

# ------------------------------------------------------------------------------
# MAIN EXECUTION
# ------------------------------------------------------------------------------

process_coordinates <- function(coord_data, clipmark_data, config, projection, 
                                 title_prefix, subtitle_suffix) {
  #' Process coordinates and generate heatmap
  #' Internal helper function used by run_analysis
  
  # Validate coordinates
  cli_alert("Validating {subtitle_suffix} coordinates...")
  
  valid_rows <- is.finite(coord_data$lat) & is.finite(coord_data$lon) &
                coord_data$lat >= -90 & coord_data$lat <= 90 &
                coord_data$lon >= -180 & coord_data$lon <= 180
  
  n_invalid <- sum(!valid_rows)
  if (n_invalid > 0) {
    cli_alert_warning("Removed {format(n_invalid, big.mark=',')} rows with invalid coordinates")
    coord_data <- coord_data[valid_rows, ]
  }
  
  if (nrow(coord_data) == 0) {
    cli_alert_danger("No valid coordinates remaining")
    return(NULL)
  }
  
  # Bin dwell times
  hex_data <- bin_dwell_times_fast(coord_data, config)
  
  # Flag hexes containing clipmarks
  hex_data <- flag_clipmark_hexes_fast(hex_data, clipmark_data, config)
  
  # Convert clipmarks to sf for extent calculation (if needed)
  clipmarks_sf <- NULL
  if (nrow(clipmark_data) > 0) {
    clipmarks_sf <- st_as_sf(clipmark_data, coords = c("lon", "lat"), crs = 4326)
  }
  
  # Create plot
  plot_title <- paste0(title_prefix, " - ", subtitle_suffix)
  
  p <- create_heatmap(
    hex_data, 
    projection,
    config,
    clipmarks_sf = clipmarks_sf,
    title = plot_title
  )
  
  return(list(plot = p, hex_data = hex_data))
}

run_analysis <- function(config, col_map) {
  cli_h1("Hex Bin Dwell Time Analysis")
  
  # Initialize database
  con <- init_database(config, col_map)
  if (is.null(con)) return(NULL)
  
  # Query frame center data
  cli_h2("Loading Frame Center Data")
  frame_data <- query_frame_centers(con, config, col_map)
  
  if (nrow(frame_data) == 0) {
    cli_alert_danger("No frame center data found for query parameters")
    dbDisconnect(con)
    return(NULL)
  }
  
  # Query aircraft position data
  cli_h2("Loading Aircraft Position Data")
  aircraft_data <- query_aircraft_positions(con, config, col_map)
  
  # Query clipmarks
  clipmark_data <- query_clipmarks(con, config, col_map)
  
  # Validate clipmark data
  if (!is.null(clipmark_data) && nrow(clipmark_data) > 0) {
    valid_clips <- is.finite(clipmark_data$lat) & is.finite(clipmark_data$lon) &
                   clipmark_data$lat >= -90 & clipmark_data$lat <= 90 &
                   clipmark_data$lon >= -180 & clipmark_data$lon <= 180
    clipmark_data <- clipmark_data[valid_clips, ]
    if (nrow(clipmark_data) == 0) clipmark_data <- NULL
  }
  
  # Choose projection based on frame center data (primary dataset)
  # Use 1st/99th percentile to avoid outliers
  valid_frame <- frame_data[is.finite(frame_data$lat) & is.finite(frame_data$lon), ]
  lat_q <- quantile(valid_frame$lat, probs = c(0.01, 0.99), na.rm = TRUE)
  lon_q <- quantile(valid_frame$lon, probs = c(0.01, 0.99), na.rm = TRUE)
  
  data_bbox <- c(
    xmin = as.numeric(lon_q[1]),
    ymin = as.numeric(lat_q[1]),
    xmax = as.numeric(lon_q[2]),
    ymax = as.numeric(lat_q[2])
  )
  class(data_bbox) <- "bbox"
  attr(data_bbox, "crs") <- st_crs(4326)
  projection <- choose_projection_from_bbox(data_bbox, config)
  
  # Create title prefix
  title_prefix <- if (!is.null(config$mission_ids)) {
    paste("Missions:", paste(config$mission_ids, collapse = ", "))
  } else if (!is.null(config$date_start) || !is.null(config$date_end)) {
    paste("Period:", config$date_start, "to", config$date_end)
  } else {
    "All Missions"
  }
  
  # Generate base filename
  base_filename <- generate_output_filename(config)
  base_name <- tools::file_path_sans_ext(base_filename)
  file_ext <- tools::file_ext(base_filename)
  
  # Results storage
  results <- list(
    plots = list(),
    output_files = list()
  )
  
  # Process Frame Center data
  cli_h2("Processing Frame Center Heatmap")
  frame_result <- process_coordinates(
    frame_data, clipmark_data, config, projection,
    title_prefix, "Camera Frame Center"
  )
  
  if (!is.null(frame_result)) {
    frame_file <- paste0(base_name, "_framecenter.", file_ext)
    cli_alert("Saving to {frame_file}...")
    ggsave(
      frame_file, 
      plot = frame_result$plot, 
      width = config$output_width, 
      height = config$output_height, 
      dpi = config$output_dpi
    )
    cli_alert_success("Saved: {frame_file}")
    results$plots$frame_center <- frame_result$plot
    results$output_files$frame_center <- frame_file
  }
  
  # Process Aircraft Position data (if enabled and available)
  if (config$generate_aircraft_map) {
    if (!is.null(aircraft_data) && nrow(aircraft_data) > 0) {
      cli_h2("Processing Aircraft Position Heatmap")
      aircraft_result <- process_coordinates(
        aircraft_data, clipmark_data, config, projection,
        title_prefix, "Aircraft Position"
      )
      
      if (!is.null(aircraft_result)) {
        aircraft_file <- paste0(base_name, "_aircraft.", file_ext)
        cli_alert("Saving to {aircraft_file}...")
        ggsave(
          aircraft_file, 
          plot = aircraft_result$plot, 
          width = config$output_width, 
          height = config$output_height, 
          dpi = config$output_dpi
        )
        cli_alert_success("Saved: {aircraft_file}")
        results$plots$aircraft <- aircraft_result$plot
        results$output_files$aircraft <- aircraft_file
      }
    } else {
      cli_alert_warning("Aircraft position data not available - skipping aircraft heatmap")
    }
  } else {
    cli_alert_info("Aircraft map generation disabled (config$generate_aircraft_map = FALSE)")
  }
  
  # Cleanup
  dbDisconnect(con)
  
  # Summary
  cli_h2("Analysis Complete")
  cli_alert_success("Generated {length(results$output_files)} heatmaps:")
  for (name in names(results$output_files)) {
    cli_alert_info("  {name}: {results$output_files[[name]]}")
  }
  
  # Return results
  return(results)
}

# ==============================================================================
# RUN
# ==============================================================================

# Uncomment the line below to generate synthetic test data:
# generate_synthetic_data(output_dir = config$data_dir, n_missions = 5)

# Run the analysis:
# result <- run_analysis(config, col_map)

# To display interactively:
# print(result$plots$frame_center)
# print(result$plots$aircraft)

cat("\n")
cli_alert_info("Script loaded. To run:")
cat("
1. Generate test data:
   
   OPTION A - Random missions (original):
   generate_synthetic_data(output_dir = config$data_dir, n_missions = 5)
   
   OPTION B - From KML route files:
   # Create .kml files with LineString paths (e.g., from Google Earth)
   # Place them in data/routes/ folder, then:
   generate_from_kml(kml_dir = 'data/routes', output_dir = config$data_dir)

2. Adjust 'config' settings at top of script

3. Run analysis:
   result <- run_analysis(config, col_map)

4. Display plots:
   print(result$plots$frame_center)   # Camera frame center heatmap
   print(result$plots$aircraft)        # Aircraft position heatmap

5. Check output files:
   result$output_files

Output generates TWO maps per run:
  - *_framecenter.png  - Where the camera was looking
  - *_aircraft.png     - Where the aircraft was flying

KML Route Tips:
  - Draw paths in Google Earth Pro, save as KML (not KMZ)
  - Each .kml file becomes one mission
  - Routes can cross datelines, poles, etc. for projection testing
")
