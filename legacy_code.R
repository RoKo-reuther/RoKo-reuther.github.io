# Runge-Kutta Method Definition
rk <- function(func, x0, y0, maxsteps, stol = 1e-8, nrow_concs, names, N_grid, tableau_species, ...) {
  
  times <- c()
  h_log <- c()
  concs <- matrix(ncol = maxsteps+1, nrow = nrow_concs)
  omega_siderite <- matrix(ncol = maxsteps+1, nrow = N_grid)
  tran.Fe2 <- matrix(ncol = maxsteps+1, nrow = N_grid)
  timescale.Fe2 <- matrix(ncol = maxsteps+1, nrow = N_grid)
  dCdt.Fe2 <- matrix(ncol = maxsteps+1, nrow = N_grid)

  for (i in seq(maxsteps+1)) {
    
    # calculate model at current time
    model_return <- func(x0, y0, ...)
    
    # set integration timestep
    h <- min(model_return$time_scales) * 0.05
    # maximum timestep: 1e-6
    h <- min(h, 1e-6)
    
    # logging
    if (i == 1) {
      rates <- as.data.frame(model_return$rates)
      speciation <- model_return$solute_equilibrium$species_conc
    } else {
      rates <- rbind(rates, as.data.frame(model_return$rates))
      speciation <- rbind(speciation, model_return$solute_equilibrium$species_conc)
    }
    
    times[i]  <- x0
    h_log[i]  <- h
    concs[,i] <- y0
    omega_siderite[,i] <- model_return$omega_siderite
    tran.Fe2[,i] <- model_return$transport$tran.Fe2$dC
    dCdt.Fe2[,i] <- model_return$dCdt_Fe2
    timescale.Fe2[,i] <- model_return$timescale_Fe2
    
    dCdt  <- model_return[[1]]
    
    # check if to break loop
    #if (mean(abs(dCdt)) < stol) break
    
    # 1st order R-K method
    x0 <- x0 + h
    y0 <- y0 + h * dCdt
    
    y0[y0 < 0] <- 0
  
  }
  
  concs_split <- list()
  for (j in seq_len(length(names))) {
    concs_split[[j]] <- concs[((j-1) * N_grid + 1) : (j * N_grid),]
  }
  names(concs_split) <- names
  
  rates_split <- list()
  for (j in seq_len(ncol(rates))) {
    rates_split[[j]] <- matrix(rates[,j], ncol = i)
  }
  names(rates_split) <- colnames(rates)
  
  speciation_split <- list()
  for (j in seq_len(ncol(speciation))) {
    speciation_split[[j]] <- matrix(speciation[,j], ncol = i)
  }
  names(speciation_split) <- tableau_species
  
  return(list(
    # method log
    times = times,
    h = h_log,
    # output of final step according to steady.1D
    y = matrix(y0, nrow = N_grid, byrow = FALSE),
    transport = model_return$transport,
    sumR = model_return$sumR,
    rates = model_return$rates,
    solute_equilibrium = model_return$solute_equilibrium,
    time_scales = model_return$time_scales,
    # time series output
    concs_split = concs_split,
    rates_split = rates_split,
    speciation_split = speciation_split,
    # to be cleaned up ...
    omega_siderite = omega_siderite,
    timescale.Fe2 = timescale.Fe2,
    tran.Fe2 = tran.Fe2,
    dCdt.Fe2 = dCdt.Fe2,
    dCdt = dCdt
  ))
}


# Runge-Kutta Method Application
timesteps <- 200

std2 <- rk(
  func = model,
  x0 = 0,
  y0 = as.vector(std3$y), # !!!
  maxsteps = timesteps,
  stol = 1e-8,
  nrow_concs = length(model_species) * N_grid,
  names = model_species,
  N_grid = N_grid,
  tableau_species = tableau_species,
  parms = parameter,
  diff_coeffs = diffusion_coefficients,
  adv_vel = advective_velocities,
  solve_equilibrium = TRUE,
  precipitation = TRUE
)

# Runge-Kutta plot Timesteps
plot(
  x = seq_along(std2$h),
  y = std2$h, type = "b",
  xlab = "timestep n",
  ylab = "integration timestep h"
)

# Runge-Kutta plot Concentration Profiles
 for (species in model_species) {
  plot_std_conc_evolution(std2$concs_split, species, 25, grid)
}

# Runge-Kutta plot Reaction Rate Development
for (rate in names(std2$rates)) {
  plot_std_rate_evolution(std2, rate, 50)
}

# Runge-Kutta Problem Analysis (Try I)
```{r, fig.width=10, fig.height=10, echo=FALSE}
#| eval: FALSE

par(mfrow = c(2, 2))

#spots <- seq(1, max_steps, 5)
spots <- seq(max_steps-5, max_steps, 1)

omega_sid <- std2$omega_siderite[, spots]
omega_sid2 <- ifelse(omega_sid>1, (omega_sid-1)^2, -(omega_sid-1)^2)
xlim1 <- c(0,50)
matplot(
  omega_sid,
  type = "l",
  #log = "y",
  ylim = c(-1,30),
  xlim = xlim1,
  main = "omega siderite"
)

matplot(
  omega_sid2,
  type = "l",
  #log = "y",
  #ylim = c(-1,2),
  xlim = xlim1,
  main = "(omega siderite-1)^2*sign"
)

matplot(
  std2$concs_split$Fe2[, spots],
  type = "l",
  #log = "y",
  xlim = xlim1,
  main = "Fe2"
)

matplot(
  std2$CO3[, spots],
  type = "l",
  log = "y",
  main = "CO3"
)

#matplot(
#  std2$concs_split$FeCO3[, spots],
#  type = "l",
#  log = "y",
#  main = "FeCO3"
#)

plot_std_omega_evolution <- function(data, N_lines, grid) {
  
  # number of available timesteps
  available_timesteps <- ncol(data)
  
  # selection of N_lines timesteps
  selected_timesteps <- floor(seq.int(1, available_timesteps, length.out = N_lines))
  
  # plot
  matplot(y = grid$x.mid,
          x = data,
          type = "l",
          ylim = c(max(grid$x.mid), min(grid$x.mid)),
          xlab = "concentration (mol/m^3_phase)",
          ylab = "depth (m)",
          log = "x",
          main = "omega siderite")
}

#plot_std_omega_evolution(std2$omega_siderite, 100, grid)
```
---

# Runge-Kutta Problem Analysis (Shiny)
```{r}
#| eval: FALSE

sliderInput("steps", "Time steps:", width = "500px", min = 1, max = timesteps, value = c(1, timesteps))
```

```{r}
#| eval: FALSE

sliderInput("depth", "Depth:", width = "500px", min = round(min(grid$x.mid), 3), max = round(max(grid$x.mid), 3), value = c(round(min(grid$x.mid), 3), round(max(grid$x.mid), 3)))
```
:::

```{r}
#| eval: FALSE
#| column: screen-inset
#| layout-ncol: 4

#plotOutput("plotOmega", height = "550px")
plotOutput("plotOmega2", height = "550px")
plotOutput("plotFe2", height = "550px")
plotOutput("plotFeCO3", height = "550px")
plotOutput("plotPrecip", height = "550px")
plotOutput("plotDiss", height = "550px")
plotOutput("plotR4", height = "550px")
plotOutput("plotTran.Fe2", height = "550px")
plotOutput("plotdCdt.Fe2", height = "550px")
#plotOutput("plottimescale.Fe2", height = "550px")
```

```{r}
#| eval: FALSE
#| context: data

timesteps <- timesteps
std <- std2
grid <- grid
omega_sid <- std$omega_siderite
omega_sid2 <- ifelse(omega_sid>1, (omega_sid-1)^2, -(omega_sid-1)^2)
s2p <- conversion_factors$s2p
```

```{r}
#| eval: FALSE
#| context: server

index_depth <- reactive({grid$x.mid > input$depth[1] & grid$x.mid < input$depth[2]})

output$plotOmega <- renderPlot({
  timespots <- seq(input$steps[1], input$steps[2], 1)
  matplot(
    x = omega_sid[index_depth(), timespots],
    y = grid$x.mid[index_depth()],
    type = "l",
    lwd = 2,
    ylim = c(input$depth[2], input$depth[1]),
    main = "omega siderite"
  )
  },
  res = 120
)

output$plotOmega2 <- renderPlot({
  timespots <- seq(input$steps[1], input$steps[2], 1)
  matplot(
    x = omega_sid2[index_depth(), timespots] + 1,
    y = grid$x.mid[index_depth()],
    type = "l",
    lwd = 2,
    ylim = c(input$depth[2], input$depth[1]),
    log = "x",
    main = "omega siderite^2"
  )
  },
  res = 120
)

output$plotFe2 <- renderPlot({
  timespots <- seq(input$steps[1], input$steps[2], 1)
  matplot(
    x = std$concs_split$Fe2[index_depth(), timespots],
    y = grid$x.mid[index_depth()],
    type = "l",
    lwd = 2,
    ylim = c(input$depth[2], input$depth[1]),
    main = "Fe2"
  )
  },
  res = 120
)

output$plotFeCO3 <- renderPlot({
  timespots <- seq(input$steps[1], input$steps[2], 1)
  matplot(
    x = std$concs_split$FeCO3[index_depth(), timespots],
    y = grid$x.mid[index_depth()],
    type = "l",
    lwd = 2,
    ylim = c(input$depth[2], input$depth[1]),
    main = "FeCO3"
  )
  },
  res = 120
)

output$plotPrecip <- renderPlot({
  timespots <- seq(input$steps[1], input$steps[2], 1)
  matplot(
    x = cbind(std$rates_split$R4[index_depth(), 1] * 4, std$rates_split$R7[index_depth(), timespots]),
    y = grid$x.mid[index_depth()],
    type = "l",
    lwd = 2,
    ylim = c(input$depth[2], input$depth[1]),
    main = "Precipitation"
  )
  },
  res = 120
)

output$plotR4 <- renderPlot({
  timespots <- seq(input$steps[1], input$steps[2], 1)
  matplot(
    x = std$rates_split$R4[index_depth(), timespots] * 4,
    y = grid$x.mid[index_depth()],
    type = "l",
    lwd = 2,
    ylim = c(input$depth[2], input$depth[1]),
    main = "R4"
  )
  },
  res = 120
)

output$plotDiss <- renderPlot({
  timespots <- seq(input$steps[1], input$steps[2], 1)
  matplot(
    x = std$rates_split$R8[index_depth(), timespots],
    y = grid$x.mid[index_depth()],
    type = "l",
    lwd = 2,
    ylim = c(input$depth[2], input$depth[1]),
    main = "Dissolution"
  )
  },
  res = 120
)

output$plotTran.Fe2 <- renderPlot({
  timespots <- seq(input$steps[1], input$steps[2], 1)
  matplot(
    x = std$tran.Fe2[index_depth(), timespots],
    y = grid$x.mid[index_depth()],
    type = "l",
    lwd = 2,
    ylim = c(input$depth[2], input$depth[1]),
    main = "Transport Fe2"
  )
  },
  res = 120
)

output$plotdCdt.Fe2 <- renderPlot({
  timespots <- seq(input$steps[1], input$steps[2], 1)
  matplot(
    x = std$dCdt.Fe2[index_depth(), timespots],
    y = grid$x.mid[index_depth()],
    type = "l",
    lwd = 2,
    ylim = c(input$depth[2], input$depth[1]),
    main = "dCdt Fe2"
  )
  },
  res = 120
)

output$plottimescale.Fe2 <- renderPlot({
  timespots <- seq(input$steps[1], input$steps[2], 1)
  matplot(
    x = std$timescale.Fe2[index_depth(), timespots],
    y = grid$x.mid[index_depth()],
    type = "l",
    lwd = 2,
    ylim = c(input$depth[2], input$depth[1]),
    main = "timescale Fe2"
  )
  },
  res = 120
)
```
---

# Search for steady state (plots std1, plots std2 etc...)

### Step 1) Steady State without Speciation

Solve steady state without equilibrium solver.
Precipitation reactions are controlled by component-totals.

```{r std1}
#| echo: true

initial <- rep(1e-4, N_grid * N_species)

std1 <- steady.1D(
  y = initial,
  func = model,
  parms = parameter, 
  dimens = N_grid,
  nspec = N_species,
  names = model_species,
  positive = TRUE,
  method = "stodes",
  diff_coeffs = diffusion_coefficients,
  adv_vel = advective_velocities,
  solve_equilibrium = FALSE,
  precipitation = TRUE
)
```

The speciation is calculated for the resulting steady-state.

```{r std1 speciation}
#| echo: true

# the first 7 model-species are tableau components
# In the model we have alkalinity (column 1) - but here we need TOT_H = -ALK
component_total <- c(-1*std1$y[,1], std1$y[,2:7])

microbenchmark(
  equilibrium  <- solve_tableau(component_total, tableau, logK, N_grid),
  times = 1
)

# Check "success" state for every layer
## 0 = at least one total_component differs more than 1e-12 (atol)
## 1 = atol reached at first try with pcfm
## 2 = atol reached at first try with newton
## 3 = atol reached at second try with pcfm
## 4 = atol reached at second try with newton
## 5 = atol reached at third try with pcfm
table(equilibrium$success)

# Check for NA results, indicating something went wrong
sum(is.na(equilibrium$species_conc))
```

#### Profiles: Component-Totals

```{r plot std1}
#| column: screen-inset
#| layout-ncol: 5
#| fig.height: 3.5
#| fig.width: 3.5

plot_std_profiles(std1)
```

#### Profiles: Speciation

```{r plot std1 speciation}
#| column: screen-inset
#| layout-ncol: 3
#| fig.width: 6

speciation <- equilibrium$species_conc
colnames(speciation) <- tableau_species

rcols <- viridis::viridis(6)

matplot(y=grid$x.mid, x = speciation[,3:5],
        type="l", lwd=2, lty=1,
        ylim=c(length,0), 
        col=c(rcols[3], rcols[6], rcols[1]),
        ylab="depth (m)", 
        xlab="mol/m3_phase"
        )
legend(
  x = "topright",
  legend = c("H2CO3", "HCO3", "CO3"),
  lty = 1,
  col = c(rcols[3], rcols[6], rcols[1]),
  lwd = 2
)

matplot(y=grid$x.mid, x = speciation[,6:9],
        type="l", lwd=2, lty=1,
        ylim=c(length,0), 
        col=c(rcols[3], rcols[6], rcols[1], rcols[4]),
        ylab="depth (m)", 
        xlab="mol/m3_phase"
        )
legend(
  x = "topright",
  legend = c("H3PO4", "H2PO4-", "HPO42-", "PO43-"),
  lty = 1,
  col = c(rcols[3], rcols[6], rcols[1], rcols[4]),
  lwd = 2
)

matplot(y=grid$x.mid, x = speciation[,10:11],
        type="l", lwd=2, lty=1,
        ylim=c(length,0), 
        col=c(rcols[3], rcols[6]),
        ylab="depth (m)", 
        xlab="mol/m3_phase"
        )
legend(
  x = "topright",
  legend = c("HNO3", "NO3-"),
  lty = 1,
  col = c(rcols[3], rcols[6]),
  lwd = 2
)

matplot(y=grid$x.mid, x = speciation[,12:13],
        type="l", lwd=2, lty=1,
        ylim=c(length,0), 
        col=c(rcols[3], rcols[6]),
        ylab="depth (m)", 
        xlab="mol/m3_phase"
        )
legend(
  x = "topright",
  legend = c("NH4+", "NH3"),
  lty = 1,
  col = c(rcols[3], rcols[6]),
  lwd = 2
)

matplot(y=grid$x.mid, x = speciation[,14:16],
        type="l", lwd=2, lty=1,
        ylim=c(length,0), 
        col=c(rcols[3], rcols[6], rcols[1]),
        ylab="depth (m)", 
        xlab="mol/m3_phase"
        )
legend(
  x = "topright",
  legend = c("H2SO4", "HSO4-", "SO42-"),
  lty = 1,
  col = c(rcols[3], rcols[6], rcols[1]),
  lwd = 2
)

matplot(y=grid$x.mid, x = speciation[,17:19],
        type="l", lwd=2, lty=1,
        ylim=c(length,0), 
        col=c(rcols[3], rcols[6], rcols[1]),
        ylab="depth (m)", 
        xlab="mol/m3_phase"
        )
legend(
  x = "topright",
  legend = c("H2S", "HS-", "S2-"),
  lty = 1,
  col = c(rcols[3], rcols[6], rcols[1]),
  lwd = 2
)
```

#### Profile: pH

```{r plot pH}
#| out.width: 100%

pH <- -log10(speciation[,1] * 10^-3)
plot(grid$x.mid ~ pH, ylim = c(0.1, 0), type = "l", ylab = "depth (m)", lwd = 2)
```

#### Profiles: Saturation

```{r plot saturation}
#| out.width: 100%

rcols <- viridis::viridis(2)

omega_siderite <- speciation[,"CO3"] * std1$y[, "Fe2"] / parameter$Ksp_siderite

plot(
    x = omega_siderite,
    y = grid$x.mid,
    ylim = c(max(grid$x.int), min(grid$x.int)),
    xlim = c(min(std1$omega_siderite, omega_siderite), max(std1$omega_siderite, omega_siderite)),
    type = "l",
    lwd = 2,
    col = rcols[1],
    xlab = "omega siderite",
    ylab = "depth (m)",
    log = "x"
  )

points(
    x = std1$omega_siderite,
    y = grid$x.mid,
    type = "l",
    lwd = 2,
    col = rcols[2]
  )

legend("bottomright", c("with CO32-", "with DIC * 0.1"), col = rcols, lwd = 2)
```

#### Profiles: Reaction Rates

```{r plot rates}
#| out.width: 100%

plot_std_rates(std1)
```


### Step 2) Steady State with Speciation

Solve steady state with equilibrium solver.
Precipitation reactions are controlled by the actual species concentrations.

```{r std3}
#| echo: true
std3 <- steady.1D( # !!!
  y = std1$y,
  func = model,
  parms = parameter,
  dimens = N_grid,
  nspec = N_species,
  names = model_species,
  positive = TRUE,
  method = "stode",
  maxiter = 100,
  diff_coeffs = diffusion_coefficients,
  adv_vel = advective_velocities,
  solve_equilibrium = TRUE,
  precipitation = TRUE
)
```

#### Profiles: Speciation

```{r}
#| column: screen-inset
#| layout-ncol: 3
#| fig.width: 6

speciation <- std2$solute_equilibrium$species_conc
colnames(speciation) <- tableau_species

rcols <- viridis::viridis(6)

matplot(y=grid$x.mid, x = speciation[,3:5],
        type="l", lwd=2, lty=1,
        ylim=c(length,0), 
        col=c(rcols[3], rcols[6], rcols[1]),
        ylab="depth (m)", 
        xlab="mol/m3_phase"
        )
legend(
  x = "topright",
  legend = c("H2CO3", "HCO3", "CO3"),
  lty = 1,
  col = c(rcols[3], rcols[6], rcols[1]),
  lwd = 2
)

matplot(y=grid$x.mid, x = speciation[,6:9],
        type="l", lwd=2, lty=1,
        ylim=c(length,0), 
        col=c(rcols[3], rcols[6], rcols[1], rcols[4]),
        ylab="depth (m)", 
        xlab="mol/m3_phase"
        )
legend(
  x = "topright",
  legend = c("H3PO4", "H2PO4-", "HPO42-", "PO43-"),
  lty = 1,
  col = c(rcols[3], rcols[6], rcols[1], rcols[4]),
  lwd = 2
)

matplot(y=grid$x.mid, x = speciation[,10:11],
        type="l", lwd=2, lty=1,
        ylim=c(length,0), 
        col=c(rcols[3], rcols[6]),
        ylab="depth (m)", 
        xlab="mol/m3_phase"
        )
legend(
  x = "topright",
  legend = c("HNO3", "NO3-"),
  lty = 1,
  col = c(rcols[3], rcols[6]),
  lwd = 2
)

matplot(y=grid$x.mid, x = speciation[,12:13],
        type="l", lwd=2, lty=1,
        ylim=c(length,0), 
        col=c(rcols[3], rcols[6]),
        ylab="depth (m)", 
        xlab="mol/m3_phase"
        )
legend(
  x = "topright",
  legend = c("NH4+", "NH3"),
  lty = 1,
  col = c(rcols[3], rcols[6]),
  lwd = 2
)

matplot(y=grid$x.mid, x = speciation[,14:16],
        type="l", lwd=2, lty=1,
        ylim=c(length,0), 
        col=c(rcols[3], rcols[6], rcols[1]),
        ylab="depth (m)", 
        xlab="mol/m3_phase"
        )
legend(
  x = "topright",
  legend = c("H2SO4", "HSO4-", "SO42-"),
  lty = 1,
  col = c(rcols[3], rcols[6], rcols[1]),
  lwd = 2
)

matplot(y=grid$x.mid, x = speciation[,17:19],
        type="l", lwd=2, lty=1,
        ylim=c(length,0), 
        col=c(rcols[3], rcols[6], rcols[1]),
        ylab="depth (m)", 
        xlab="mol/m3_phase"
        )
legend(
  x = "topright",
  legend = c("H2S", "HS-", "S2-"),
  lty = 1,
  col = c(rcols[3], rcols[6], rcols[1]),
  lwd = 2
)
```

#### Profile: pH

```{r}
#| out.width: 100%

pH <- -log10(speciation[,1] * 10^-3)
plot(grid$x.mid ~ pH, ylim = c(0.1, 0), type = "l", ylab = "depth (m)", lwd = 2)
```

### Steady State Visualization

```{mermaid}
flowchart LR
  a([SO4]) -- 10 --> b(SO4);
  b -- 10 --> c(H2S);
  c -- 10 --> b;
  c -- 10 --> d(S0);
  b -- 10 --> e([SO4]);
```