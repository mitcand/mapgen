# ==============================================================================
# FOCUS HUNTING DETECTION - PHASE 2: VIDEO FRAME EXTRACTION
# ==============================================================================
#
# Extracts before/after frames from video files for focus events detected in Phase 1.
# Uses ffmpeg to extract specific frames based on timestamps.
#
# Prerequisites:
#   1. Run focus_analysis.R (Phase 1) first to generate focus_events.csv
#   2. Place video files in the video_dir folder
#   3. ffmpeg must be installed and accessible from command line
#
# Video filename format: <missionid>_<group>_<camera>_<HHMMSS>.mpg
#   - HHMMSS is the start time of the video
#   - Date is derived from missionid
#
# ==============================================================================

# ------------------------------------------------------------------------------
# CONFIGURATION
# ------------------------------------------------------------------------------

config <- list(
  # === INPUT ===
  events_file = "output/focus_events.csv",      # From Phase 1
  video_dir = "data/videos",                     # Folder containing .mpg files
  
  # === EXTRACTION SETTINGS ===
  before_offset_sec = 1.0,                       # Seconds before event to extract
  after_offset_sec = 1.0,                        # Seconds after event to extract
  
  # Event types to extract frames for
  extract_event_types = c("autofocus", "manual", "locked_on_target"),
  
  # Limit extractions (NULL = all, or integer for max per event type)
  max_events_per_type = NULL,
  
  # Only extract on-target events
  on_target_only = FALSE,
  
  # === OUTPUT ===
  output_dir = "output/focus_frames",
  generate_html_gallery = TRUE,
  
  # === VIDEO PARSING ===
  # Video filename pattern: <missionid>_<group>_<camera>_<HHMMSS>.mpg
  # Adjust these indices if your naming convention is different
  filename_mission_idx = 1,                      # Position of mission_id in filename (split by _)
  filename_camera_idx = 3,                       # Position of camera identifier
  filename_time_idx = 4,                         # Position of HHMMSS time
  
  # Camera name mapping (filename identifier -> display name)
  camera_mapping = list(
    "CAM1" = "Camera 1",
    "CAM2" = "Camera 2",
    "cam1" = "Camera 1",
    "cam2" = "Camera 2"
  )
)

# ------------------------------------------------------------------------------
# PACKAGES
# ------------------------------------------------------------------------------

required_packages <- c("dplyr", "lubridate", "glue", "cli", "stringr")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# ------------------------------------------------------------------------------
# UTILITY FUNCTIONS
# ------------------------------------------------------------------------------

check_ffmpeg <- function() {
  result <- system("ffmpeg -version", intern = TRUE, ignore.stderr = TRUE)
  if (length(result) == 0) {
    cli_alert_danger("ffmpeg not found. Please install ffmpeg and ensure it's in your PATH.")
    return(FALSE)
  }
  cli_alert_success("ffmpeg found: {result[1]}")
  return(TRUE)
}

parse_video_filename <- function(filename, config) {
  # Parse video filename to extract mission, camera, start time
  # Format: <missionid>_<group>_<camera>_<HHMMSS>.mpg
  
  base_name <- tools::file_path_sans_ext(basename(filename))
  parts <- str_split(base_name, "_")[[1]]
  
  if (length(parts) < config$filename_time_idx) {
    return(NULL)
  }
  
  mission_id <- parts[config$filename_mission_idx]
  camera_raw <- parts[config$filename_camera_idx]
  time_str <- parts[config$filename_time_idx]
  
  # Map camera identifier to display name
  camera <- config$camera_mapping[[camera_raw]]
  if (is.null(camera)) {
    camera <- camera_raw
  }
  
  # Parse time (HHMMSS)
  if (nchar(time_str) == 6) {
    hours <- as.integer(substr(time_str, 1, 2))
    minutes <- as.integer(substr(time_str, 3, 4))
    seconds <- as.integer(substr(time_str, 5, 6))
    start_seconds <- hours * 3600 + minutes * 60 + seconds
  } else {
    return(NULL)
  }
  
  return(list(
    filename = filename,
    mission_id = mission_id,
    camera = camera,
    time_str = time_str,
    start_seconds = start_seconds
  ))
}

find_video_for_event <- function(event, video_catalog, config) {
  # Find the video file that contains a given event timestamp
  
  # Filter to matching mission and camera
  matching <- video_catalog %>%
    filter(
      mission_id == event$mission_id,
      camera == event$camera
    )
  
  if (nrow(matching) == 0) {
    return(NULL)
  }
  
  # Convert event time to seconds since midnight
  event_time <- as.POSIXct(event$timestamp)
  event_seconds <- hour(event_time) * 3600 + minute(event_time) * 60 + second(event_time)
  
  # Find video where event falls within video duration
  # Assume videos are ~30 min each if we don't know duration
  # The video with start_seconds <= event_seconds is the candidate
  matching <- matching %>%
    filter(start_seconds <= event_seconds) %>%
    arrange(desc(start_seconds))
  
  if (nrow(matching) == 0) {
    return(NULL)
  }
  
  # Best match is the video that started most recently before the event
  best_match <- matching[1, ]
  
  # Calculate offset within video
  offset_seconds <- event_seconds - best_match$start_seconds
  
  return(list(
    video_file = best_match$filename,
    offset_seconds = offset_seconds
  ))
}

extract_frame <- function(video_file, offset_seconds, output_file) {
  # Use ffmpeg to extract a single frame
  
  # Format offset as HH:MM:SS.mmm
  offset_formatted <- sprintf("%02d:%02d:%06.3f",
                              floor(offset_seconds / 3600),
                              floor((offset_seconds %% 3600) / 60),
                              offset_seconds %% 60)
  
  cmd <- glue('ffmpeg -y -ss {offset_formatted} -i "{video_file}" -frames:v 1 -q:v 2 "{output_file}" 2>/dev/null')
  
  result <- system(cmd, intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE)
  
  return(result == 0)
}

# ------------------------------------------------------------------------------
# MAIN EXTRACTION
# ------------------------------------------------------------------------------

build_video_catalog <- function(config) {
  cli_alert("Scanning video directory...")
  
  video_files <- list.files(
    config$video_dir, 
    pattern = "\\.(mpg|mpeg|mp4|avi|mov)$", 
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  if (length(video_files) == 0) {
    cli_alert_danger("No video files found in {config$video_dir}")
    return(NULL)
  }
  
  cli_alert_info("Found {length(video_files)} video files")
  
  # Parse each video file
  catalog <- lapply(video_files, function(f) {
    parsed <- parse_video_filename(f, config)
    if (!is.null(parsed)) {
      data.frame(
        filename = parsed$filename,
        mission_id = parsed$mission_id,
        camera = parsed$camera,
        time_str = parsed$time_str,
        start_seconds = parsed$start_seconds,
        stringsAsFactors = FALSE
      )
    } else {
      NULL
    }
  })
  
  catalog <- bind_rows(catalog)
  
  if (nrow(catalog) == 0) {
    cli_alert_danger("Could not parse any video filenames. Check filename_*_idx settings.")
    return(NULL)
  }
  
  cli_alert_success("Cataloged {nrow(catalog)} videos")
  
  # Summary
  for (cam in unique(catalog$camera)) {
    n_cam <- sum(catalog$camera == cam)
    cli_alert_info("  {cam}: {n_cam} videos")
  }
  
  return(catalog)
}

run_frame_extraction <- function(config) {
  cli_h1("Focus Frame Extraction - Phase 2")
  
  # Check ffmpeg
  if (!check_ffmpeg()) {
    return(NULL)
  }
  
  # Load events from Phase 1
  if (!file.exists(config$events_file)) {
    cli_alert_danger("Events file not found: {config$events_file}")
    cli_alert_info("Run Phase 1 (focus_analysis.R) first")
    return(NULL)
  }
  
  events <- read.csv(config$events_file, stringsAsFactors = FALSE)
  events$timestamp <- as.POSIXct(events$timestamp)
  
  cli_alert_success("Loaded {nrow(events)} events from Phase 1")
  
  # Filter events
  if (config$on_target_only) {
    events <- events %>% filter(on_target == TRUE)
    cli_alert_info("Filtered to {nrow(events)} on-target events")
  }
  
  events <- events %>% filter(event_type %in% config$extract_event_types)
  cli_alert_info("Filtered to {nrow(events)} events of types: {paste(config$extract_event_types, collapse=', ')}")
  
  # Limit if configured
  if (!is.null(config$max_events_per_type)) {
    events <- events %>%
      group_by(event_type) %>%
      slice_head(n = config$max_events_per_type) %>%
      ungroup()
    cli_alert_info("Limited to {nrow(events)} events ({config$max_events_per_type} per type)")
  }
  
  if (nrow(events) == 0) {
    cli_alert_warning("No events to process after filtering")
    return(NULL)
  }
  
  # Build video catalog
  video_catalog <- build_video_catalog(config)
  if (is.null(video_catalog)) {
    return(NULL)
  }
  
  # Create output directory
  dir.create(config$output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Process each event
  cli_h2("Extracting Frames")
  
  extraction_results <- list()
  
  pb <- cli_progress_bar("Extracting frames", total = nrow(events))
  
  for (i in seq_len(nrow(events))) {
    event <- events[i, ]
    cli_progress_update()
    
    # Find video
    video_match <- find_video_for_event(event, video_catalog, config)
    
    if (is.null(video_match)) {
      cli_alert_warning("No video found for event {event$event_id} ({event$camera}, {event$timestamp})")
      next
    }
    
    # Create output filenames
    event_prefix <- glue("event_{event$event_id}_{event$event_type}")
    before_file <- file.path(config$output_dir, paste0(event_prefix, "_before.jpg"))
    after_file <- file.path(config$output_dir, paste0(event_prefix, "_after.jpg"))
    
    # Calculate offsets
    before_offset <- video_match$offset_seconds - config$before_offset_sec
    after_offset <- video_match$offset_seconds + config$after_offset_sec
    
    # Extract frames
    before_ok <- FALSE
    after_ok <- FALSE
    
    if (before_offset >= 0) {
      before_ok <- extract_frame(video_match$video_file, before_offset, before_file)
    }
    
    after_ok <- extract_frame(video_match$video_file, after_offset, after_file)
    
    extraction_results[[i]] <- data.frame(
      event_id = event$event_id,
      event_type = event$event_type,
      camera = event$camera,
      zoom_level = event$zoom_level,
      timestamp = as.character(event$timestamp),
      before_focus = event$before_focus,
      after_focus = event$after_focus,
      on_target = event$on_target,
      before_file = ifelse(before_ok, before_file, NA),
      after_file = ifelse(after_ok, after_file, NA),
      video_file = video_match$video_file,
      stringsAsFactors = FALSE
    )
  }
  
  cli_progress_done()
  
  results_df <- bind_rows(extraction_results)
  
  # Summary
  n_success <- sum(!is.na(results_df$before_file) | !is.na(results_df$after_file))
  cli_alert_success("Extracted frames for {n_success} of {nrow(events)} events")
  
  # Save results
  results_file <- file.path(config$output_dir, "extraction_results.csv")
  write.csv(results_df, results_file, row.names = FALSE)
  cli_alert_success("Saved: {results_file}")
  
  # Generate HTML gallery
  if (config$generate_html_gallery && n_success > 0) {
    generate_html_gallery(results_df, config)
  }
  
  cli_h1("Phase 2 Complete")
  
  return(results_df)
}

generate_html_gallery <- function(results, config) {
  cli_alert("Generating HTML gallery...")
  
  # Filter to successful extractions
  results <- results %>%
    filter(!is.na(before_file) | !is.na(after_file))
  
  if (nrow(results) == 0) {
    cli_alert_warning("No successful extractions for gallery")
    return(NULL)
  }
  
  # Build HTML
  html_parts <- c(
    "<!DOCTYPE html>",
    "<html>",
    "<head>",
    "<title>Focus Events Gallery</title>",
    "<style>",
    "body { font-family: Arial, sans-serif; margin: 20px; background: #f5f5f5; }",
    "h1 { color: #333; }",
    ".event { background: white; border-radius: 8px; padding: 15px; margin: 20px 0; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }",
    ".event-header { font-size: 14px; color: #666; margin-bottom: 10px; }",
    ".event-type { display: inline-block; padding: 3px 8px; border-radius: 4px; color: white; font-weight: bold; }",
    ".autofocus { background: #e41a1c; }",
    ".manual { background: #ff7f00; }",
    ".locked_on_target { background: #4daf4a; }",
    ".on-target-badge { background: #377eb8; color: white; padding: 2px 6px; border-radius: 3px; font-size: 11px; margin-left: 10px; }",
    ".frames { display: flex; gap: 20px; justify-content: center; }",
    ".frame { text-align: center; }",
    ".frame img { max-width: 600px; max-height: 400px; border: 1px solid #ddd; }",
    ".frame-label { font-weight: bold; margin-top: 5px; }",
    ".focus-value { color: #666; font-size: 12px; }",
    ".filter-bar { margin: 20px 0; padding: 10px; background: #333; border-radius: 4px; }",
    ".filter-bar button { margin: 0 5px; padding: 5px 15px; border: none; border-radius: 3px; cursor: pointer; }",
    ".filter-bar button.active { background: #4CAF50; color: white; }",
    "</style>",
    "</head>",
    "<body>",
    "<h1>Focus Events Gallery</h1>",
    glue("<p>Generated: {Sys.time()} | Total events: {nrow(results)}</p>"),
    "<div class='filter-bar'>",
    "<button onclick=\"filterEvents('all')\" class='active'>All</button>",
    "<button onclick=\"filterEvents('autofocus')\">Autofocus</button>",
    "<button onclick=\"filterEvents('manual')\">Manual</button>",
    "<button onclick=\"filterEvents('locked_on_target')\">Locked</button>",
    "</div>"
  )
  
  # Add each event
  for (i in seq_len(nrow(results))) {
    r <- results[i, ]
    
    on_target_badge <- if (r$on_target) "<span class='on-target-badge'>ON TARGET</span>" else ""
    
    before_img <- if (!is.na(r$before_file) && file.exists(r$before_file)) {
      glue("<div class='frame'><img src='{basename(r$before_file)}' /><div class='frame-label'>Before</div><div class='focus-value'>Focus: {round(r$before_focus, 1)}%</div></div>")
    } else {
      "<div class='frame'><div style='width:300px;height:200px;background:#eee;display:flex;align-items:center;justify-content:center;'>No frame</div></div>"
    }
    
    after_img <- if (!is.na(r$after_file) && file.exists(r$after_file)) {
      glue("<div class='frame'><img src='{basename(r$after_file)}' /><div class='frame-label'>After</div><div class='focus-value'>Focus: {round(r$after_focus, 1)}%</div></div>")
    } else {
      "<div class='frame'><div style='width:300px;height:200px;background:#eee;display:flex;align-items:center;justify-content:center;'>No frame</div></div>"
    }
    
    event_html <- glue("
      <div class='event' data-type='{r$event_type}'>
        <div class='event-header'>
          <span class='event-type {r$event_type}'>{toupper(r$event_type)}</span>
          {on_target_badge}
          Event #{r$event_id} | {r$camera} | {r$zoom_level} | {r$timestamp}
        </div>
        <div class='frames'>
          {before_img}
          {after_img}
        </div>
      </div>
    ")
    
    html_parts <- c(html_parts, event_html)
  }
  
  # Add JavaScript for filtering
  html_parts <- c(html_parts,
    "<script>",
    "function filterEvents(type) {",
    "  document.querySelectorAll('.filter-bar button').forEach(b => b.classList.remove('active'));",
    "  event.target.classList.add('active');",
    "  document.querySelectorAll('.event').forEach(e => {",
    "    if (type === 'all' || e.dataset.type === type) {",
    "      e.style.display = 'block';",
    "    } else {",
    "      e.style.display = 'none';",
    "    }",
    "  });",
    "}",
    "</script>",
    "</body>",
    "</html>"
  )
  
  # Write HTML
  html_file <- file.path(config$output_dir, "gallery.html")
  writeLines(html_parts, html_file)
  
  cli_alert_success("Saved: {html_file}")
  cli_alert_info("Open in browser to view the focus event gallery")
  
  return(html_file)
}

# ==============================================================================
# RUN
# ==============================================================================

# Uncomment to run:
# results <- run_frame_extraction(config)
