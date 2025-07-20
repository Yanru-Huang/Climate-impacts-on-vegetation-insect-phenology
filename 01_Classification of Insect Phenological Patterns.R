# Required packages
library(lme4)      
library(mclust)    
library(pastecs)  
library(lubridate) 
library(ggplot2)  
library(gridExtra) 
library(maps)      
library(dplyr)   
library(signal)   
library(pracma)    
library(cluster)  

# Input / output directories
in_dir  <- ".../data"     # Input data directory
out_dir <- ".../output" # Output directory
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Color palette (Tableau10)
tableau10 <- c(
  "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
  "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC"
)

# Load world map data
world_map <- map_data("world")

# Prepare summary collector
summary_list <- list()

# Process each CSV
csv_files <- list.files(in_dir, pattern = "\\.csv$", full.names = TRUE)
for (file_path in csv_files) {
  # Read data
  tab <- read.csv(file_path, stringsAsFactors = FALSE)
  speciesname <- tools::file_path_sans_ext(basename(file_path))
  
  # Convert columns\ n  tab$julian_ins           <- as.numeric(tab$julian_ins)   # Julian day for insects
  tab$Latitude  <- as.numeric(tab$Latitude)   # Latitude as numeric
  tab$Longitude <- as.numeric(tab$Longitude)  # Longitude as numeric
  tab$year      <- factor(tab$year)           # Year as factor
  
  # Filter to Europe bounding box
  tab <- subset(
    tab,
    Latitude  >= 32 & Latitude  <= 72 &
      Longitude >= -15 & Longitude <= 36
  )
  
  # Remove outliers in julian_ins (±3 standard deviations)
  mu    <- mean(tab$julian_ins, na.rm = TRUE)
  sigma <- sd(tab$julian_ins, na.rm = TRUE)
  tab   <- subset(tab, julian_ins >= (mu - 3*sigma) & julian_ins <= (mu + 3*sigma))
  
  # 1) Remove spatial + year trend using mixed-effects model
  model1    <- lmer(julian_ins ~ Latitude + Longitude + (1|year),
                    data = tab, REML = FALSE)
  tab$resi <- residuals(model1)  # Extract residuals
  
  # Kernel density estimation
  dens_obj <- density(tab$resi, kernel = "gaussian")
  dens_x   <- dens_obj$x
  dens_y   <- dens_obj$y
  
  # Compute smoothing window parameters
  dx            <- dens_x[2] - dens_x[1]
  min_dist_pts  <- round(30 / dx)
  window_len    <- round(min_dist_pts / 2)
  if (window_len %% 2 == 0) window_len <- window_len + 1
  
  # 2) Savitzky-Golay smoothing
  sg_y <- sgolayfilt(dens_y, p = 3, n = window_len)
  
  # 3) Peak detection
  peaks <- findpeaks(
    sg_y,
    minpeakheight   = 0.2 * max(sg_y),
    minpeakdistance = min_dist_pts
  )
  nbmax <- if (!is.null(peaks)) nrow(peaks) else 1  # Number of modes detected
  
  # Initialize clustering result
  tab$cate <- NA
  kopt     <- NA
  sc_opt   <- NA
  
  # Only perform K-means clustering if more than one mode
  if (nbmax > 1) {
    # Preserve row indices for mapping back
    tab$id <- seq_len(nrow(tab))
    
    resi_mu <- mean(tab$resi, na.rm = TRUE)
    resi_sd <- sd(tab$resi, na.rm = TRUE)
    tab2 <- tab  # Use full data without further filtering
    
    # Prepare silhouette scores
    sil_scores <- rep(NA, nbmax)
    
    # Subsample if dataset is large, else use all data
    N <- nrow(tab2)
    if (N > 5000) {
      set.seed(123)
      idx_samp  <- sample.int(N, 5000)
      resi_samp <- tab2$resi[idx_samp]
    } else {
      idx_samp  <- seq_len(N)
      resi_samp <- tab2$resi
    }
    
    # Compute distance matrix on sampled data
    dmat <- dist(matrix(resi_samp, ncol = 1))
    
    # Evaluate k = 2..nbmax using silhouette
    for (k in 2:nbmax) {
      km_res      <- kmeans(resi_samp, centers = k, nstart = 25)
      sil_obj     <- silhouette(km_res$cluster, dmat)
      sil_scores[k] <- mean(sil_obj[, "sil_width"] )
    }
    
    # Select optimal k
    kopt   <- which.max(sil_scores)
    sc_opt <- sil_scores[kopt]
    
    # Final clustering on full data with optimal k
    final_km <- kmeans(tab2$resi, centers = kopt, nstart = 25)
    tab$cate[tab2$id] <- final_km$cluster
  } else {
    # Single mode: assign all to cluster 1
    tab$cate <- factor(1)
    kopt      <- 1
  }
  
  # 4) Reorder cluster labels by increasing mean julian_ins
  mean_j <- tapply(tab$julian_ins, tab$cate, mean, na.rm = TRUE)
  new_lv <- names(sort(mean_j))
  mapping <- setNames(seq_along(new_lv), new_lv)
  tab$cate <- mapping[as.character(tab$cate)]
  tab <- tab[order(tab$cate), ]
  
  # 5) Plotting
  tab <- droplevels(tab[!is.na(tab$cate), ])
  ncat <- length(unique(tab$cate))
  cols <- tableau10[1:ncat]
  
  pl1 <- ggplot(tab, aes(x = julian_ins, fill = factor(cate))) +
    geom_histogram(binwidth = 5, alpha = 0.7, position = 'identity') +
    scale_fill_manual("Mode", values = cols) +
    theme_bw() + labs(x = "Julian day", y = "Count") +
    theme(legend.position = "right")
  
  pl2 <- ggplot() +
    geom_polygon(data = world_map,
                 aes(x = long, y = lat, group = group),
                 fill = "gray95", color = "gray70", size = 0.3) +
    geom_point(data = tab,
               aes(x = Longitude, y = Latitude,
                   color = factor(cate)),
               size = 1, alpha = 0.7) +
    scale_color_manual("Mode", values = cols) +
    coord_quickmap(xlim = c(-15, 36), ylim = c(32, 72)) +
    theme_bw() + labs(x = "Longitude", y = "Latitude") +
    theme(legend.position = "none", panel.grid = element_blank())
  
  pl3 <- ggplot(tab, aes(x = Latitude, y = julian_ins,
                         color = factor(cate))) +
    geom_point(size = 1, alpha = 0.7) +
    scale_color_manual("Mode", values = cols) +
    theme_bw() + labs(x = "Latitude", y = "Julian day") +
    theme(legend.position = "none")
  
  g <- arrangeGrob(
    pl1, pl2, pl3,
    layout_matrix = rbind(c(1,1), c(2,3)),
    heights = c(1, 1.2)
  )
  
  # 6) Save outputs
  ggsave(
    filename = file.path(out_dir, paste0(speciesname, ".png")),
    plot = g, width = 4, height = 4, dpi = 300
  )
  
  # Save classified table
  write.csv(
    tab,
    file = file.path(out_dir, paste0(speciesname, ".csv")),
    row.names = FALSE
  )
  
  # Collect summary
  summary_list[[speciesname]] <- data.frame(
    species = speciesname,
    nbmax   = nbmax,
    kopt    = kopt,
    sc_opt  = sc_opt,
    stringsAsFactors = FALSE
  )
}

