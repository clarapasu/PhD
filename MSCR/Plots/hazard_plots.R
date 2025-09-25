# ==============================================================================
# MSCR HAZARD PLOTS 
# By Clara Panchaud
# ==============================================================================
#This code produces an example MSCR hazard plot, as found in Figure 3.1 of the thesis. 
# ==============================================================================

# Load required libraries
library(secr)
library(ggplot2)
library(viridis)

# Source the hazard function code 
source("Functions/Sim_Func.R")

# ==============================================================================
# SETUP AND PARAMETERS
# ==============================================================================

h0<-0 # Fixed baseline hazard parameter (log-scale)
sigma_values <- c(-0.7,0,0.7) # Three different spatial scales (log-scale)
beta_values <- c(-1,0,1) # Three different memory parameters

t<- 1 # Time after the capture occurred (days)
      # To obtain Figure 3.1, run this code with t=1/24, 0.5, 1 and 7

# Generate a regular 3x3 trap grid
traps <- make.grid(3, 3, spacex = 1, detector = "multi")
trap <- as.data.frame(traps)

# Define a fine mesh for plotting the hazard
mesh <- make.mask(traps, buffer = 2, type = "traprect", nx = 70)

z <- data.frame(x =trap[1,1], y =trap[1,2]) # Location of the previous capture
s <- data.frame(x = 1.5, y = 0.5) # Activity centre location
all_results <- data.frame() # Data frame initialisation

# ==============================================================================
# HAZARD CALCULATION
# ==============================================================================

# Loop over the combinations of parameters
for (beta in beta_values) {
  for (sigma in sigma_values) {
    theta <- c(h0,sigma,beta) # Define the current parameters
    
    # Calculate the hazard function over the whole mesh
    f <- apply(mesh, 1, hazard, theta = theta, t = t, z = z, s = s, m = 1)
    
    # Combine the results into a data frame
    df <- data.frame(as.data.frame(mesh), unlist(f))
    df$beta_label <- paste("Beta = ", beta, sep = "")
    df$sigma_label <- paste("Sigma = ", sigma, sep = "")
    
    # Add to all results
    all_results <- rbind(all_results, df)
  }
}

# Rename the column containing the hazard value 
names(all_results)[names(all_results) == "unlist.f."] <- "f"


# Combine and name the points that we want on the plot
trap$Type <- "Camera traps"
s$Type <- "Activity centre"
z$Type <- "Last observed location"
points <- rbind(trap, s, z)
all_results$sigma_label <- factor(all_results$sigma_label,
                                  levels = paste("Sigma = ", sigma_values, sep = ""))



g4 <- ggplot(data = all_results, aes(x = x, y = y)) +
  # Raster for smooth color gradients
  geom_raster(aes(fill = f)) +
  scale_fill_viridis(
    option = "mako", 
    direction = -1,
  ) +
  # Camera traps, activity centre, and last observation points
  geom_point(data = points,
             aes(x = x, y = y, shape = Type, color = Type),
             size = 5, stroke = 1) +
  scale_shape_manual(values = c(
    "Camera traps" = 4,
    "Activity centre" = 20,
    "Last observed location" = 4
  )) +
  scale_color_manual(values = c(
    "Camera traps" = "black",
    "Activity centre" = "orange2",
    "Last observed location" = "#DC143C"
  )) +
  guides(shape = guide_legend(title = NULL),
         color = guide_legend(title = NULL)) +
  
  coord_fixed() +
  scale_x_continuous(limits = c(-1, 3), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-1, 3), expand = c(0, 0)) +
  
  # Faceting for all beta × sigma combinations
  facet_grid(beta_label ~ sigma_label) +
  
  theme_minimal(base_size = 15) +
  theme(
    panel.grid = element_blank(),
    panel.spacing = unit(1.2, "lines"),
    legend.position = "bottom",
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Hazard function after 1 day", # Change if "t" was modified in the setup
    x = NULL,
    y = NULL
  )  + theme(legend.position = "none")

g4



# Save the figure 
#ggsave(
#  filename = "Figures/hazard_1day.pdf",  
#  plot = g4,                               
#  device = cairo_pdf,                            
#  width = 9,                                    
#  height = 7,                                   
#  units = "in",
#  dpi = 600
#)
