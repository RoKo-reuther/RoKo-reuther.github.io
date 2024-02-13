save_and_upload_results <- function(result, result_name, result_type) {

  # result_type: "std" | "trans"

  # save as .rds
  saveRDS(result, paste0("last_results/", result_name, ".rds"))

  # save results as .csv
  ## std
  if (result_type == "std") {
    # species concentrations
    write.csv(result$y, paste0("../temp_drop/", result_name, "_species.csv"))
    # transport
    write.csv(result$transport, paste0("../temp_drop/", result_name, "_transport.csv"))
    # sumR
    write.csv(result$sumR, paste0("../temp_drop/", result_name, "_sumR.csv"))
    # reaction rates
    write.csv(result$rates, paste0("../temp_drop/", result_name, "_rates.csv"))
    # omega
    write.csv(result$omega, paste0("../temp_drop/", result_name, "_omega.csv"))
    # speciation
    try(write.csv(result$solute_equilibrium$species_conc, paste0("../temp_drop/", result_name, "_speciation.csv")), silent = TRUE)
    
  } else if (result_type == "trans") {

  }

  # upload to github repo: git@github.com:RoKo-reuther/temp_drop.git
  system("cd ../temp_drop && git add . && git commit -m 'update' && git push")
}

std0 <- readRDS("last_results/std0.rds")
std1 <- readRDS("last_results/std1.rds")
std2 <- readRDS("last_results/std2.rds")
#save_and_upload_results(std0, "std0", "std")
#save_and_upload_results(std1, "std1", "std")
#save_and_upload_results(std2, "std2", "std")


download_model_result <- function(url) {
  download.file(url, "temp.csv")
  read.csv("temp.csv", sep = ",")
}

grid <- download_model_result("https://raw.githubusercontent.com/RoKo-reuther/temp_drop/master/grid.csv")

std0 <- list(
  species    = download_model_result("https://raw.githubusercontent.com/RoKo-reuther/temp_drop/master/std0_species.csv"),
  transport  = download_model_result("https://raw.githubusercontent.com/RoKo-reuther/temp_drop/master/std0_transport.csv"),
  sumR       = download_model_result("https://raw.githubusercontent.com/RoKo-reuther/temp_drop/master/std0_sumR.csv"),
  rates      = download_model_result("https://raw.githubusercontent.com/RoKo-reuther/temp_drop/master/std0_rates.csv"),
  omega      = download_model_result("https://raw.githubusercontent.com/RoKo-reuther/temp_drop/master/std0_omega.csv")
)

std1 <- list(
  species    = download_model_result("https://raw.githubusercontent.com/RoKo-reuther/temp_drop/master/std1_species.csv"),
  transport  = download_model_result("https://raw.githubusercontent.com/RoKo-reuther/temp_drop/master/std1_transport.csv"),
  sumR       = download_model_result("https://raw.githubusercontent.com/RoKo-reuther/temp_drop/master/std1_sumR.csv"),
  rates      = download_model_result("https://raw.githubusercontent.com/RoKo-reuther/temp_drop/master/std1_rates.csv"),
  omega      = download_model_result("https://raw.githubusercontent.com/RoKo-reuther/temp_drop/master/std1_omega.csv")
)

std2 <- list(
  species    = download_model_result("https://raw.githubusercontent.com/RoKo-reuther/temp_drop/master/std2_species.csv"),
  transport  = download_model_result("https://raw.githubusercontent.com/RoKo-reuther/temp_drop/master/std2_transport.csv"),
  sumR       = download_model_result("https://raw.githubusercontent.com/RoKo-reuther/temp_drop/master/std2_sumR.csv"),
  rates      = download_model_result("https://raw.githubusercontent.com/RoKo-reuther/temp_drop/master/std2_rates.csv"),
  omega      = download_model_result("https://raw.githubusercontent.com/RoKo-reuther/temp_drop/master/std2_omega.csv"),
  speciation = download_model_result("https://raw.githubusercontent.com/RoKo-reuther/temp_drop/master/std2_speciation.csv")
)


# HELPER FUNCTIONS TO FIND VALUES

index_by_depth <- function(start, end) {
  which(grid[,2] >= start & grid[,2] <= end)
}

minmaxdepth <- function() {c(min(grid[,2]), max(grid[,2]))}

std.conc <- function(std, species, depth = minmaxdepth()) {
  # std.conc(std = std1, species = "Fe2", depth = c(0.02,0.03))
  
  data <- as.matrix(std$species)
  if (!is.null(std$speciation)) {
    data <- cbind(as.matrix(std$species), as.matrix(std$speciation))
  }

  depth <- index_by_depth(depth[1], depth[2])

  data[depth, species]
}

std.transport <- function(std, species, depth = minmaxdepth()) {
  # std.transport(std = std1, species = "Fe2", depth = c(0.02,0.03))
  std$transport[[paste0("tran.", species, ".dC")]][index_by_depth(depth[1], depth[2])]
}

std.flux_up <- function(std, species) {
  # std.flux_up(std = std1, species = "Fe2")
  std$transport[[paste0("tran.", species, ".flux.up")]][1]
}

std.flux_down <- function(std, species) {
  # std.flux_down(std = std1, species = "Fe2")
  std$transport[[paste0("tran.", species, ".flux.down")]][1]
}

std.sumR <- function(std, species, depth = minmaxdepth()) {
  # std.sumR(std = std1, species = "Fe2", depth = c(0.02,0.03))
  std$sumR[index_by_depth(depth[1], depth[2]), paste0("R.", species)]
}

std.rate <- function(std, rate, depth = minmaxdepth()) {
  # std.rate(std = std1, rate = "R12", depth = c(0.02,0.03))
  std$rates[index_by_depth(depth[1], depth[2]), rate]
}

std.omega <- function(std, species, depth = minmaxdepth()) {
  # std.omega(std = std1, species = "Siderite", depth = c(0.02, 0.03))
  std$omega[index_by_depth(depth[1], depth[2]), species]
}

std.plot <- function(..., depth = minmaxdepth(), legendpos = "bottomright") {
  # plots an vector of x values and more over the depth specified for this elements
  # do not specify a depth for single elements!!!
  # std.plot(std1 = std.conc(std1, "Fe2"), std2 = std.conc(std2, "Fe2"))

  index <- index_by_depth(depth[1], depth[2])

  depth <- grid[index, 2]

  elements <- list(...)

  elements <- as.data.frame(elements)

  cols <- palette.colors(ncol(elements))

  plot(
    x    = elements[index, 1],
    y    = depth,
    ylim = c(max(depth), min(depth)),
    xlim = c(min(elements), max(elements)),
    type = "l",
    lwd  = 2,
    col  = cols[1],
    ylab = "depth",
    xlab = ""
  )

  for (i in seq_len(ncol(elements))[-1]) {

    points(
      x    = elements[index, i],
      y    = depth,
      type = "l",
      lwd  = 2,
      col  = cols[i]
    )

  }

  legend(
    legendpos,
    legend = colnames(elements),
    col    = cols,
    lwd = 2
  )
  
}








