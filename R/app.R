# =============================================================================
# Potato GWAS Interactive Dashboard
# Optimizing GWAS to Map the Genetic Basis of Complex Traits in Crops
# Author: Atharv Yeole | MS Bioinformatics
# =============================================================================

library(shiny)
library(shinydashboard)
library(plotly)
library(DT)
library(ggplot2)
library(dplyr)
library(tidyr)

# =============================================================================
# GENERATE DEMO DATA (mirrors real analysis outputs)
# Replace these sections with read.csv() calls to your actual output files
# =============================================================================

set.seed(42)

# --- Variety names (sample from real dataset naming conventions) ---
varieties <- paste0("Variety_", sprintf("%03d", 1:282))

# --- BLUEs data for 5 traits across 282 varieties ---
blues_data <- data.frame(
  ID = varieties,
  yield_kg = rnorm(282, mean = 4.5, sd = 1.8),
  Max_Max_height_m = rnorm(282, mean = 0.65, sd = 0.12),
  Max_Ave_height_m = rnorm(282, mean = 0.50, sd = 0.10),
  Max_area_plant_m2 = rnorm(282, mean = 0.15, sd = 0.05),
  Max_volume_plant_m3 = rnorm(282, mean = 0.08, sd = 0.03),
  stringsAsFactors = FALSE
)
blues_data$yield_kg <- pmax(blues_data$yield_kg, 0.5)

# --- Raw phenotype data (simulated trial structure) ---
raw_pheno <- expand.grid(
  Variety = varieties,
  System = c("Conv", "Org"),
  Replicate = 1:2,
  Year = 2016:2018
) %>%
  mutate(
    Yield_Kg = rnorm(n(), mean = 4.5, sd = 2.0) %>% pmax(0.3),
    Max_Max_height_m = ifelse(Year == 2016, NA, rnorm(n(), 0.65, 0.15)),
    Max_Ave_height_m = ifelse(Year == 2016, NA, rnorm(n(), 0.50, 0.12)),
    Max_area_plant_m2 = ifelse(Year == 2016, NA, rnorm(n(), 0.15, 0.06)),
    Max_volume_plant_m3 = ifelse(Year == 2016, NA, rnorm(n(), 0.08, 0.04))
  )

# --- PCA data ---
pca_scores <- data.frame(
  ID = varieties,
  PC1 = rnorm(282, 0, 2.5),
  PC2 = rnorm(282, 0, 1.8),
  PC3 = rnorm(282, 0, 1.2),
  PC4 = rnorm(282, 0, 1.0),
  PC5 = rnorm(282, 0, 0.8)
)
# Create some clustering structure
cluster_assignment <- sample(1:4, 282, replace = TRUE, prob = c(0.3, 0.25, 0.25, 0.2))
for (i in 1:282) {
  offset <- switch(cluster_assignment[i],
                   `1` = c(3, 2), `2` = c(-3, 2), `3` = c(-2, -3), `4` = c(3, -2))
  pca_scores$PC1[i] <- pca_scores$PC1[i] + offset[1]
  pca_scores$PC2[i] <- pca_scores$PC2[i] + offset[2]
}

# --- PCA eigenvalues ---
pca_eigenvalues <- data.frame(
  PC = 1:20,
  Eigenvalue = c(8.2, 5.1, 3.4, 2.8, 2.1, 1.5, 1.2, 1.0, 0.9, 0.85,
                 0.8, 0.75, 0.72, 0.70, 0.68, 0.65, 0.63, 0.61, 0.59, 0.57),
  Cumulative_Var = cumsum(c(8.2, 5.1, 3.4, 2.8, 2.1, 1.5, 1.2, 1.0, 0.9, 0.85,
                            0.8, 0.75, 0.72, 0.70, 0.68, 0.65, 0.63, 0.61, 0.59, 0.57))
)
pca_eigenvalues$Pct_Var <- pca_eigenvalues$Eigenvalue / sum(pca_eigenvalues$Eigenvalue) * 100
pca_eigenvalues$Cum_Pct <- cumsum(pca_eigenvalues$Pct_Var)

# --- DAPC cluster membership probabilities ---
dapc_membership <- data.frame(
  ID = varieties,
  Cluster = factor(cluster_assignment)
)
# Generate membership probabilities
membership_probs <- matrix(0, nrow = 282, ncol = 4)
for (i in 1:282) {
  dominant <- cluster_assignment[i]
  probs <- runif(4, 0.02, 0.15)
  probs[dominant] <- runif(1, 0.5, 0.95)
  probs <- probs / sum(probs)
  membership_probs[i, ] <- probs
}
dapc_membership <- cbind(dapc_membership, as.data.frame(membership_probs))
colnames(dapc_membership)[3:6] <- paste0("Cluster_", 1:4)

# --- GWAS results (SNP effect tables) ---
n_snps <- 3000
chromosomes <- sample(1:12, n_snps, replace = TRUE)
positions <- sapply(chromosomes, function(ch) sample(1e6:8e7, 1))

generate_gwas_results <- function(method_name) {
  models <- c("additive", "general", "1-dom-alt", "1-dom-ref", "2-dom-alt", "2-dom-ref")
  results_list <- lapply(models, function(model) {
    log10_pvals <- rexp(n_snps, rate = 1.5)
    # Inject significant signals on chr 6 and 7
    sig_idx_chr6 <- which(chromosomes == 6)[1:3]
    sig_idx_chr7 <- which(chromosomes == 7)[1:3]
    sig_idx <- c(sig_idx_chr6, sig_idx_chr7)
    sig_idx <- sig_idx[!is.na(sig_idx)]
    if (model %in% c("general", "additive")) {
      log10_pvals[sig_idx] <- runif(length(sig_idx), 5.5, 8.0)
    }
    bonf_threshold <- -log10(0.05 / n_snps)
    data.frame(
      SNP = paste0("solcap_snp_", sprintf("%05d", 1:n_snps)),
      Chromosome = chromosomes,
      Position = positions,
      Model = model,
      Effect = rnorm(n_snps, 0, 0.5),
      Score = log10_pvals,
      Significant = log10_pvals >= bonf_threshold,
      Method = method_name,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, results_list)
}

gwas_pca <- generate_gwas_results("PCA")
gwas_dapc <- generate_gwas_results("DAPC")
gwas_structure <- generate_gwas_results("STRUCTURE")
gwas_all <- rbind(gwas_pca, gwas_dapc, gwas_structure)

# --- Significant SNPs summary ---
sig_snps <- gwas_all %>% filter(Significant == TRUE)

# --- QTL boxplot data ---
sig_snp_names <- unique(sig_snps$SNP)[1:min(6, length(unique(sig_snps$SNP)))]
qtl_data <- expand.grid(
  SNP = sig_snp_names,
  Genotype = c("AAAA", "AAAB", "AABB", "ABBB", "BBBB")
) %>%
  mutate(
    Yield = rnorm(n(), mean = 4.5, sd = 1.2) +
      as.numeric(factor(Genotype)) * 0.3
  )


# =============================================================================
# UI
# =============================================================================

ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(
    title = span(
      icon("dna"), " Potato GWAS Dashboard",
      style = "font-size: 16px; font-weight: 600;"
    ),
    titleWidth = 320
  ),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      id = "tabs",
      menuItem("Project Overview", tabName = "overview", icon = icon("house")),
      menuItem("Phenotype Explorer", tabName = "phenotype", icon = icon("chart-bar")),
      menuItem("Population Structure", tabName = "population", icon = icon("diagram-project")),
      menuItem("GWAS Results", tabName = "gwas", icon = icon("magnifying-glass-chart")),
      menuItem("Significant SNPs", tabName = "significant", icon = icon("star")),
      menuItem("Methods Summary", tabName = "methods", icon = icon("book"))
    ),
    hr(),
    div(
      style = "padding: 10px 15px; color: #718096; font-size: 12px;",
      p(icon("user"), " Atharv Yeole"),
      p(icon("graduation-cap"), " MS Bioinformatics"),
      p(icon("calendar"), " 2025")
    )
  ),

  dashboardBody(
    tags$head(
      tags$style(HTML("
        /* Global styling */
        .content-wrapper { background-color: #edf2f7; }
        .main-header .logo { background-color: #1e3a5f !important; color: #fff !important; }
        .main-header .navbar { background-color: #2b6cb0 !important; }
        .main-sidebar { background-color: #1e3a5f !important; }
        .sidebar-menu > li > a { color: #bee3f8 !important; }
        .sidebar-menu > li.active > a { background-color: #3182ce !important; color: #fff !important; }
        .sidebar-menu > li:hover > a { background-color: #2c5282 !important; }

        /* Cards & boxes */
        .box { background-color: #ffffff; border: 1px solid #e2e8f0; border-top: 3px solid #3182ce; box-shadow: 0 2px 8px rgba(0,0,0,0.06); }
        .box-header { color: #1e3a5f; }
        .box-body { color: #2d3748; }
        .box-title { font-weight: 600; color: #1e3a5f; }
        .small-box { border-radius: 8px; }
        .small-box .icon { font-size: 70px; top: 5px; }

        /* Info boxes */
        .info-box { background: #ffffff; border: 1px solid #e2e8f0; border-radius: 8px; }
        .info-box-number { color: #2b6cb0; font-size: 24px; }
        .info-box-text { color: #2d3748; font-size: 13px; }

        /* Value boxes */
        .small-box { background: linear-gradient(135deg, #2b6cb0 0%, #3182ce 100%) !important; }
        .small-box h3 { color: #ffffff; font-size: 32px; }
        .small-box p { color: #bee3f8; font-size: 14px; }
        .small-box .icon { color: rgba(255,255,255,0.15); }

        /* Tab content */
        .nav-tabs-custom > .tab-content { background: #ffffff; color: #2d3748; }
        .nav-tabs-custom > .nav-tabs > li.active > a { background-color: #fff; color: #2b6cb0; border-top: 3px solid #3182ce; }
        .nav-tabs-custom > .nav-tabs > li > a { color: #718096; }

        /* DataTables */
        .dataTables_wrapper { color: #2d3748; }
        table.dataTable { color: #2d3748 !important; }
        table.dataTable thead th { color: #fff !important; background: #2b6cb0 !important; }
        table.dataTable tbody tr { background: #ffffff !important; }
        table.dataTable tbody tr:hover { background: #ebf4ff !important; }
        .dataTables_filter input { border: 1px solid #bee3f8 !important; border-radius: 4px; }

        /* Select inputs */
        .selectize-input { background: #ffffff !important; color: #2d3748 !important; border: 1px solid #bee3f8 !important; }
        .selectize-dropdown { background: #ffffff !important; color: #2d3748 !important; }
        .selectize-dropdown-content .option { color: #2d3748 !important; }
        .selectize-dropdown-content .active { background: #ebf4ff !important; color: #1e3a5f !important; }

        /* Slider */
        .irs--shiny .irs-bar { background: #3182ce; }
        .irs--shiny .irs-handle { border-color: #3182ce; }

        /* Custom classes */
        .stat-highlight { color: #2b6cb0; font-weight: 700; font-size: 28px; }
        .pipeline-step { background: linear-gradient(135deg, #ebf4ff, #ffffff); border-radius: 10px;
                         padding: 20px; margin: 8px; border-left: 4px solid #3182ce; box-shadow: 0 2px 6px rgba(0,0,0,0.05); }
        .pipeline-step h4 { color: #2b6cb0; margin-top: 0; }
        .pipeline-step p { color: #4a5568; margin-bottom: 0; }
        .overview-text { color: #2d3748; font-size: 15px; line-height: 1.7; }
        .overview-text strong { color: #2b6cb0; }

        /* Plotly theme fix */
        .js-plotly-plot .plotly .modebar { background: transparent !important; }

        /* Scrollbar */
        ::-webkit-scrollbar { width: 8px; }
        ::-webkit-scrollbar-track { background: #edf2f7; }
        ::-webkit-scrollbar-thumb { background: #90cdf4; border-radius: 4px; }

        /* HR line in sidebar */
        .main-sidebar hr { border-color: #2c5282; }
      "))
    ),

    tabItems(

      # =====================================================================
      # TAB 1: PROJECT OVERVIEW
      # =====================================================================
      tabItem(
        tabName = "overview",
        fluidRow(
          valueBox(282, "Potato Varieties", icon = icon("seedling"), color = "red", width = 3),
          valueBox("~8K", "SNP Markers", icon = icon("dna"), color = "blue", width = 3),
          valueBox(5, "Traits Analyzed", icon = icon("chart-line"), color = "green", width = 3),
          valueBox(3, "Pop. Structure Methods", icon = icon("diagram-project"), color = "yellow", width = 3)
        ),
        fluidRow(
          box(
            title = "Project Summary", width = 12, solidHeader = TRUE, status = "primary",
            div(
              class = "overview-text",
              h3("Optimizing GWAS to Map the Genetic Basis of Complex Traits in Crops",
                 style = "color: #2b6cb0; margin-bottom: 20px;"),
              p("This project performs ", strong("Genome-Wide Association Studies (GWAS)"),
                " on a diversity panel of ~282 tetraploid potato (", em("Solanum tuberosum"),
                ") varieties using ", strong("~8,000 SNP markers"), " (Illumina Infinium array)."),
              p("Phenotype data comes from ", strong("field trials conducted 2016-2018"),
                " at Nafferton Farm, Newcastle University, including UAV-derived canopy architecture measurements."),
              p("The study evaluates whether incorporating high-throughput phenotyping data with
                multiple population structure correction methods improves genetic mapping resolution
                in this autotetraploid crop.")
            )
          )
        ),
        fluidRow(
          box(
            title = "Analysis Pipeline", width = 12, solidHeader = TRUE, status = "primary",
            fluidRow(
              column(3, div(class = "pipeline-step",
                h4(icon("filter"), " Step 1: Data QC"),
                p("Filter SNPs (≥20% missing), calculate PIC values, handle phenotype NAs via replicate imputation.")
              )),
              column(3, div(class = "pipeline-step",
                h4(icon("calculator"), " Step 2: BLUEs"),
                p("Mixed-effects models: Trait ~ Genotype + (1|System) + (1|Year). Model selection via LRT, AIC/BIC.")
              )),
              column(3, div(class = "pipeline-step",
                h4(icon("diagram-project"), " Step 3: Pop. Structure"),
                p("Three methods compared: PCA, DAPC (K=4 clusters), and STRUCTURE (K=3). Population covariates exported.")
              )),
              column(3, div(class = "pipeline-step",
                h4(icon("magnifying-glass-chart"), " Step 4: GWAS"),
                p("GWASpoly with 6 genetic models × 3 pop. structure methods × 5 traits. Bonferroni correction (α=0.05).")
              ))
            )
          )
        ),
        fluidRow(
          box(
            title = "Key Findings", width = 12, solidHeader = TRUE, status = "primary",
            div(class = "overview-text",
              p(icon("check-circle", style = "color: #3182ce;"),
                " 6 significant SNPs identified for Yield on ", strong("chromosomes 6 and 7")),
              p(icon("check-circle", style = "color: #3182ce;"),
                " Results consistent across all three population structure correction methods (PCA, DAPC, STRUCTURE)"),
              p(icon("check-circle", style = "color: #3182ce;"),
                " General model showed the most significant associations in the DAPC-corrected analysis"),
              p(icon("check-circle", style = "color: #3182ce;"),
                " DAPC identified 4 genetic clusters in the potato diversity panel with optimized PC retention (n=5)")
            )
          )
        )
      ),

      # =====================================================================
      # TAB 2: PHENOTYPE EXPLORER
      # =====================================================================
      tabItem(
        tabName = "phenotype",
        fluidRow(
          box(
            title = "Controls", width = 3, solidHeader = TRUE, status = "primary",
            selectInput("pheno_trait", "Select Trait:",
              choices = c("Yield (kg)" = "yield_kg",
                          "Max Height (m)" = "Max_Max_height_m",
                          "Average Height (m)" = "Max_Ave_height_m",
                          "Canopy Area (m²)" = "Max_area_plant_m2",
                          "Canopy Volume (m³)" = "Max_volume_plant_m3")),
            hr(),
            selectInput("pheno_trait_x", "Trait Correlation X:",
              choices = c("Yield (kg)" = "yield_kg",
                          "Max Height (m)" = "Max_Max_height_m",
                          "Average Height (m)" = "Max_Ave_height_m",
                          "Canopy Area (m²)" = "Max_area_plant_m2",
                          "Canopy Volume (m³)" = "Max_volume_plant_m3")),
            selectInput("pheno_trait_y", "Trait Correlation Y:",
              choices = c("Max Height (m)" = "Max_Max_height_m",
                          "Yield (kg)" = "yield_kg",
                          "Average Height (m)" = "Max_Ave_height_m",
                          "Canopy Area (m²)" = "Max_area_plant_m2",
                          "Canopy Volume (m³)" = "Max_volume_plant_m3")),
            hr(),
            sliderInput("n_top_bottom", "Top/Bottom N varieties:", min = 3, max = 15, value = 5)
          ),
          column(
            width = 9,
            fluidRow(
              box(title = "BLUEs Distribution", width = 6, solidHeader = TRUE, status = "primary",
                  plotlyOutput("blues_histogram", height = "350px")),
              box(title = "Trait Correlation", width = 6, solidHeader = TRUE, status = "primary",
                  plotlyOutput("trait_correlation", height = "350px"))
            ),
            fluidRow(
              box(title = "Top Performing Varieties", width = 6, solidHeader = TRUE, status = "primary",
                  plotlyOutput("top_varieties", height = "350px")),
              box(title = "Bottom Performing Varieties", width = 6, solidHeader = TRUE, status = "primary",
                  plotlyOutput("bottom_varieties", height = "350px"))
            )
          )
        ),
        fluidRow(
          box(title = "Trait Correlation Heatmap", width = 12, solidHeader = TRUE, status = "primary",
              plotlyOutput("corr_heatmap", height = "400px"))
        )
      ),

      # =====================================================================
      # TAB 3: POPULATION STRUCTURE
      # =====================================================================
      tabItem(
        tabName = "population",
        fluidRow(
          box(
            title = "PCA Controls", width = 3, solidHeader = TRUE, status = "primary",
            selectInput("pc_x", "X Axis:", choices = paste0("PC", 1:5), selected = "PC1"),
            selectInput("pc_y", "Y Axis:", choices = paste0("PC", 1:5), selected = "PC2"),
            checkboxInput("color_by_cluster", "Color by DAPC Cluster", value = TRUE)
          ),
          box(title = "PCA Scatter Plot", width = 9, solidHeader = TRUE, status = "primary",
              plotlyOutput("pca_scatter", height = "450px"))
        ),
        fluidRow(
          box(title = "PCA Scree Plot (Eigenvalues)", width = 6, solidHeader = TRUE, status = "primary",
              plotlyOutput("scree_plot", height = "350px")),
          box(title = "DAPC Cluster Sizes", width = 6, solidHeader = TRUE, status = "primary",
              plotlyOutput("cluster_pie", height = "350px"))
        ),
        fluidRow(
          box(title = "DAPC Membership Probabilities (STRUCTURE-like Plot)", width = 12,
              solidHeader = TRUE, status = "primary",
              plotlyOutput("membership_bar", height = "350px"))
        )
      ),

      # =====================================================================
      # TAB 4: GWAS RESULTS
      # =====================================================================
      tabItem(
        tabName = "gwas",
        fluidRow(
          box(
            title = "GWAS Controls", width = 3, solidHeader = TRUE, status = "primary",
            selectInput("gwas_method", "Population Structure Method:",
              choices = c("PCA", "DAPC", "STRUCTURE"), selected = "DAPC"),
            selectInput("gwas_model", "Genetic Model:",
              choices = c("additive", "general", "1-dom-alt", "1-dom-ref", "2-dom-alt", "2-dom-ref"),
              selected = "general"),
            hr(),
            checkboxInput("show_threshold", "Show Bonferroni Threshold", value = TRUE),
            checkboxInput("highlight_sig", "Highlight Significant SNPs", value = TRUE),
            hr(),
            p(style = "color: #b0b0b0; font-size: 12px;",
              icon("info-circle"), " Bonferroni threshold = -log10(0.05 / n_SNPs)")
          ),
          column(
            width = 9,
            box(title = "Manhattan Plot", width = 12, solidHeader = TRUE, status = "primary",
                plotlyOutput("manhattan_plot", height = "400px")),
            box(title = "QQ Plot", width = 12, solidHeader = TRUE, status = "primary",
                plotlyOutput("qq_plot", height = "350px"))
          )
        ),
        fluidRow(
          box(title = "Compare Methods: Significant SNP Counts", width = 12,
              solidHeader = TRUE, status = "primary",
              plotlyOutput("method_comparison", height = "300px"))
        )
      ),

      # =====================================================================
      # TAB 5: SIGNIFICANT SNPs
      # =====================================================================
      tabItem(
        tabName = "significant",
        fluidRow(
          valueBox(
            nrow(sig_snps %>% distinct(SNP)), "Unique Significant SNPs",
            icon = icon("star"), color = "red", width = 4),
          valueBox(
            length(unique(sig_snps$Chromosome)), "Chromosomes with Hits",
            icon = icon("dna"), color = "blue", width = 4),
          valueBox(
            "Bonferroni", "Correction Method (α=0.05)",
            icon = icon("shield-halved"), color = "green", width = 4)
        ),
        fluidRow(
          box(title = "Significant SNPs Table", width = 12, solidHeader = TRUE, status = "primary",
              DTOutput("sig_table"))
        ),
        fluidRow(
          box(
            title = "QTL Boxplots (Genotype vs Yield)", width = 12,
            solidHeader = TRUE, status = "primary",
            selectInput("qtl_snp", "Select SNP:",
              choices = sig_snp_names, width = "300px"),
            plotlyOutput("qtl_boxplot", height = "400px")
          )
        )
      ),

      # =====================================================================
      # TAB 6: METHODS SUMMARY
      # =====================================================================
      tabItem(
        tabName = "methods",
        fluidRow(
          box(
            title = "Detailed Methods", width = 12, solidHeader = TRUE, status = "primary",
            div(
              class = "overview-text",
              h3("1. Data & Study Design", style = "color: #2b6cb0;"),
              p("A diversity panel of ", strong("282 tetraploid potato varieties"),
                " was phenotyped across 3 years (2016-2018) at Nafferton Farm, Newcastle University,
                 under two cultivation systems (Organic and Conventional), with 2 replicates each.
                 Genotyping was performed using the ", strong("Illumina Infinium 8K SNP array"),
                " (Sharma et al., 2018)."),

              h3("2. Phenotype Preprocessing", style = "color: #2b6cb0;"),
              p("Five traits were selected: Yield (kg), Max Max Height (m), Max Average Height (m),
                 Canopy Area (m²), and Canopy Volume (m³). Yield was measured in all 3 years;
                 canopy traits only in 2017-2018 (UAV-derived). Missing values were imputed using
                 same-variety replicate data. ", strong("BLUEs (Best Linear Unbiased Estimates)"),
                " were calculated using lme4 mixed-effects models with Genotype as fixed effect
                 and System + Year as random effects. Model selection was validated with
                 log-likelihood ratio tests, AIC, and BIC."),

              h3("3. Genotype QC & Filtering", style = "color: #2b6cb0;"),
              p("SNPs with ≥20% missing data were removed (consistent with Sharma et al., 2018).
                 For DAPC and STRUCTURE comparison, a ", strong("PIC (Polymorphism Information Content) > 0.4"),
                " filter was applied to reduce computational demands while retaining informative markers."),

              h3("4. Population Structure Analysis", style = "color: #2b6cb0;"),
              p(strong("PCA:"), " Performed using glPca() from the adegenet package on the genlight object.
                 Scree plot analysis identified the elbow at 4 PCs; 5 PCs retained for GWAS."),
              p(strong("DAPC:"), " K-means clustering (K=4 based on BIC) followed by DAPC with
                 optimized PC retention. A-score optimization suggested 7 PCs, visual PCA inspection
                 suggested 4-5, and Thia (2023) k-1 criterion suggested 3. Final decision: 5 PCs
                 (compromise). 3 discriminant functions retained."),
              p(strong("STRUCTURE:"), " External analysis with K=3 using the STRUCTURE software.
                 Q-matrix (membership probabilities) imported for GWAS correction."),

              h3("5. GWAS with GWASpoly", style = "color: #2b6cb0;"),
              p("GWASpoly was used for tetraploid-aware GWAS with 6 genetic models:
                 additive, general, 1-dom-alt, 1-dom-ref, 2-dom-alt, and 2-dom-ref.
                 Each trait was analyzed with all three population structure corrections.
                 Kinship (K) matrix was computed internally. Significance determined by
                 ", strong("Bonferroni correction (α = 0.05)"), "."),

              h3("6. Post-GWAS Analysis", style = "color: #2b6cb0;"),
              p("Significant SNPs were compiled into effect tables. QTL boxplots were generated
                 for significant markers showing genotype-phenotype relationships across the
                 5 tetraploid dosage classes (AAAA, AAAB, AABB, ABBB, BBBB).
                 Chromosome-specific QQ plots were generated for chromosomes with significant hits."),

              h3("Software & Packages", style = "color: #2b6cb0;"),
              p("R v4.x | GWASpoly | adegenet | lme4 | tidyverse | ggplot2 | ggpubr | STRUCTURE")
            )
          )
        )
      )
    )
  )
)


# =============================================================================
# SERVER
# =============================================================================

server <- function(input, output, session) {

  # Custom light blue plotly theme
  dark_layout <- function(p, ...) {
    p %>% layout(
      paper_bgcolor = "rgba(255,255,255,0)",
      plot_bgcolor = "rgba(237,242,247,0.5)",
      font = list(color = "#2d3748"),
      xaxis = list(gridcolor = "rgba(0,0,0,0.07)", zerolinecolor = "rgba(0,0,0,0.12)"),
      yaxis = list(gridcolor = "rgba(0,0,0,0.07)", zerolinecolor = "rgba(0,0,0,0.12)"),
      ...
    )
  }

  # ----- PHENOTYPE TAB -----

  output$blues_histogram <- renderPlotly({
    trait <- input$pheno_trait
    trait_label <- names(which(c("yield_kg" = "Yield (kg)", "Max_Max_height_m" = "Max Height (m)",
      "Max_Ave_height_m" = "Average Height (m)", "Max_area_plant_m2" = "Canopy Area (m²)",
      "Max_volume_plant_m3" = "Canopy Volume (m³)") == trait))
    if (length(trait_label) == 0) trait_label <- trait

    plot_ly(blues_data, x = ~get(trait), type = "histogram",
            marker = list(color = "#3182ce", line = list(color = "#2b6cb0", width = 0.5)),
            nbinsx = 30) %>%
      dark_layout(xaxis = list(title = trait_label), yaxis = list(title = "Count"),
                  title = list(text = paste("BLUEs Distribution:", trait_label), font = list(size = 14)))
  })

  output$trait_correlation <- renderPlotly({
    plot_ly(blues_data, x = ~get(input$pheno_trait_x), y = ~get(input$pheno_trait_y),
            type = "scatter", mode = "markers",
            marker = list(color = "#4299e1", size = 5, opacity = 0.6),
            text = ~ID, hoverinfo = "text+x+y") %>%
      dark_layout(xaxis = list(title = input$pheno_trait_x),
                  yaxis = list(title = input$pheno_trait_y),
                  title = list(text = "Trait Correlation", font = list(size = 14)))
  })

  output$top_varieties <- renderPlotly({
    n <- input$n_top_bottom
    trait <- input$pheno_trait
    top <- blues_data %>% arrange(desc(!!sym(trait))) %>% head(n)
    plot_ly(top, x = ~reorder(ID, get(trait)), y = ~get(trait), type = "bar",
            marker = list(color = "#2b6cb0")) %>%
      dark_layout(xaxis = list(title = "", tickangle = -45),
                  yaxis = list(title = trait),
                  title = list(text = paste("Top", n, "Varieties"), font = list(size = 14)))
  })

  output$bottom_varieties <- renderPlotly({
    n <- input$n_top_bottom
    trait <- input$pheno_trait
    bottom <- blues_data %>% arrange(!!sym(trait)) %>% head(n)
    plot_ly(bottom, x = ~reorder(ID, -get(trait)), y = ~get(trait), type = "bar",
            marker = list(color = "#90cdf4")) %>%
      dark_layout(xaxis = list(title = "", tickangle = -45),
                  yaxis = list(title = trait),
                  title = list(text = paste("Bottom", n, "Varieties"), font = list(size = 14)))
  })

  output$corr_heatmap <- renderPlotly({
    trait_cols <- c("yield_kg", "Max_Max_height_m", "Max_Ave_height_m",
                    "Max_area_plant_m2", "Max_volume_plant_m3")
    trait_labels <- c("Yield", "Max Height", "Avg Height", "Area", "Volume")
    cor_mat <- cor(blues_data[, trait_cols], use = "complete.obs")
    plot_ly(x = trait_labels, y = trait_labels, z = cor_mat, type = "heatmap",
            colorscale = list(c(0, "#bee3f8"), c(0.5, "#4299e1"), c(1, "#1e3a5f")),
            text = round(cor_mat, 3), hoverinfo = "text") %>%
      dark_layout(title = list(text = "Trait Correlation Matrix (BLUEs)", font = list(size = 14)))
  })

  # ----- POPULATION STRUCTURE TAB -----

  output$pca_scatter <- renderPlotly({
    pc_x <- input$pc_x
    pc_y <- input$pc_y
    df <- merge(pca_scores, dapc_membership[, c("ID", "Cluster")], by = "ID")
    if (input$color_by_cluster) {
      colors <- c("#D81B60", "#1E88E5", "#FFC107", "#004D40")
      plot_ly(df, x = ~get(pc_x), y = ~get(pc_y), color = ~Cluster,
              colors = colors, type = "scatter", mode = "markers",
              marker = list(size = 7, opacity = 0.7), text = ~ID) %>%
        dark_layout(xaxis = list(title = pc_x), yaxis = list(title = pc_y),
                    title = list(text = paste(pc_x, "vs", pc_y), font = list(size = 14)))
    } else {
      plot_ly(df, x = ~get(pc_x), y = ~get(pc_y), type = "scatter", mode = "markers",
              marker = list(size = 7, color = "#3182ce", opacity = 0.6), text = ~ID) %>%
        dark_layout(xaxis = list(title = pc_x), yaxis = list(title = pc_y),
                    title = list(text = paste(pc_x, "vs", pc_y), font = list(size = 14)))
    }
  })

  output$scree_plot <- renderPlotly({
    plot_ly(pca_eigenvalues, x = ~PC, y = ~Eigenvalue, type = "scatter", mode = "lines+markers",
            line = list(color = "#2b6cb0"), marker = list(color = "#2b6cb0", size = 8),
            name = "Eigenvalue") %>%
      add_trace(y = ~Cum_Pct / max(Cum_Pct) * max(Eigenvalue), name = "Cumulative %",
                line = list(color = "#90cdf4", dash = "dash"),
                marker = list(color = "#90cdf4", size = 5)) %>%
      dark_layout(xaxis = list(title = "Principal Component"),
                  yaxis = list(title = "Eigenvalue"),
                  title = list(text = "PCA Scree Plot", font = list(size = 14)))
  })

  output$cluster_pie <- renderPlotly({
    counts <- dapc_membership %>% count(Cluster)
    colors <- c("#D81B60", "#1E88E5", "#FFC107", "#004D40")
    plot_ly(counts, labels = ~paste("Cluster", Cluster), values = ~n, type = "pie",
            marker = list(colors = colors),
            textinfo = "label+percent", textfont = list(color = "#333")) %>%
      dark_layout(title = list(text = "DAPC Cluster Distribution", font = list(size = 14)))
  })

  output$membership_bar <- renderPlotly({
    long_df <- dapc_membership %>%
      arrange(Cluster, ID) %>%
      mutate(order = row_number()) %>%
      pivot_longer(cols = starts_with("Cluster_"), names_to = "Cluster_name", values_to = "Probability")
    colors <- c("#D81B60", "#1E88E5", "#FFC107", "#004D40")
    plot_ly(long_df, x = ~order, y = ~Probability, color = ~Cluster_name,
            colors = colors, type = "bar") %>%
      layout(barmode = "stack") %>%
      dark_layout(xaxis = list(title = "Cultivar (sorted by cluster)", showticklabels = FALSE),
                  yaxis = list(title = "Membership Probability"),
                  title = list(text = "DAPC Membership Probabilities", font = list(size = 14)))
  })

  # ----- GWAS TAB -----

  gwas_filtered <- reactive({
    gwas_all %>%
      filter(Method == input$gwas_method, Model == input$gwas_model)
  })

  output$manhattan_plot <- renderPlotly({
    df <- gwas_filtered()
    df <- df %>%
      arrange(Chromosome, Position) %>%
      mutate(index = row_number())

    chr_centers <- df %>% group_by(Chromosome) %>%
      summarise(center = mean(index), .groups = "drop")

    bonf <- -log10(0.05 / n_snps)

    chr_colors <- rep(c("#2b6cb0", "#90cdf4"), length.out = 12)

    p <- plot_ly()

    for (chr in sort(unique(df$Chromosome))) {
      chr_df <- df %>% filter(Chromosome == chr)
      color <- chr_colors[chr]
      if (input$highlight_sig) {
        chr_df$color <- ifelse(chr_df$Significant, "#d62828", color)
        chr_df$size <- ifelse(chr_df$Significant, 8, 4)
      } else {
        chr_df$color <- color
        chr_df$size <- 4
      }
      p <- p %>% add_trace(data = chr_df, x = ~index, y = ~Score,
                            type = "scatter", mode = "markers",
                            marker = list(color = ~color, size = ~size, opacity = 0.7),
                            text = ~paste("SNP:", SNP, "<br>Chr:", Chromosome,
                                          "<br>Pos:", Position, "<br>Score:", round(Score, 3)),
                            hoverinfo = "text", showlegend = FALSE)
    }

    if (input$show_threshold) {
      p <- p %>% add_trace(x = c(min(df$index), max(df$index)), y = c(bonf, bonf),
                            type = "scatter", mode = "lines",
                            line = list(color = "#d62828", dash = "dash", width = 2),
                            showlegend = FALSE)
    }

    p %>% dark_layout(
      xaxis = list(title = "Chromosome",
                   tickvals = chr_centers$center, ticktext = chr_centers$Chromosome),
      yaxis = list(title = "-log10(p-value)"),
      title = list(text = paste("Manhattan Plot |", input$gwas_method, "|", input$gwas_model),
                   font = list(size = 14))
    )
  })

  output$qq_plot <- renderPlotly({
    df <- gwas_filtered()
    observed <- sort(df$Score, decreasing = TRUE)
    n <- length(observed)
    expected <- -log10(ppoints(n))
    expected <- sort(expected, decreasing = TRUE)

    qq_df <- data.frame(Expected = expected, Observed = observed)

    plot_ly(qq_df, x = ~Expected, y = ~Observed, type = "scatter", mode = "markers",
            marker = list(color = "#3182ce", size = 3, opacity = 0.5)) %>%
      add_trace(x = c(0, max(expected)), y = c(0, max(expected)),
                type = "scatter", mode = "lines",
                line = list(color = "#888", dash = "dash"), showlegend = FALSE) %>%
      dark_layout(xaxis = list(title = "Expected -log10(p)"),
                  yaxis = list(title = "Observed -log10(p)"),
                  title = list(text = paste("QQ Plot |", input$gwas_method, "|", input$gwas_model),
                               font = list(size = 14)))
  })

  output$method_comparison <- renderPlotly({
    comparison <- gwas_all %>%
      filter(Significant == TRUE) %>%
      group_by(Method, Model) %>%
      summarise(Count = n_distinct(SNP), .groups = "drop")

    colors <- c("PCA" = "#1e3a5f", "DAPC" = "#3182ce", "STRUCTURE" = "#90cdf4")

    plot_ly(comparison, x = ~Model, y = ~Count, color = ~Method,
            colors = colors, type = "bar") %>%
      layout(barmode = "dodge") %>%
      dark_layout(xaxis = list(title = "Genetic Model"),
                  yaxis = list(title = "# Significant SNPs"),
                  title = list(text = "Significant SNPs by Method & Model", font = list(size = 14)))
  })

  # ----- SIGNIFICANT SNPs TAB -----

  output$sig_table <- renderDT({
    sig_display <- sig_snps %>%
      select(SNP, Chromosome, Position, Model, Method, Score, Effect) %>%
      mutate(Score = round(Score, 3), Effect = round(Effect, 4)) %>%
      arrange(Chromosome, Position)

    datatable(sig_display, options = list(
      pageLength = 15,
      dom = 'fltip',
      initComplete = JS("function(settings, json) {",
        "$(this.api().table().header()).css({'background-color': '#2b6cb0', 'color': '#ffffff'});",
        "$(this.api().table().body()).css({'background-color': '#ffffff', 'color': '#333333'});",
        "}")
    ), filter = "top", rownames = FALSE)
  })

  output$qtl_boxplot <- renderPlotly({
    snp <- input$qtl_snp
    df <- qtl_data %>% filter(SNP == snp)

    plot_ly(df, x = ~Genotype, y = ~Yield, type = "box",
            marker = list(color = "#2b6cb0"),
            line = list(color = "#2b6cb0"),
            fillcolor = "rgba(64,145,108,0.3)") %>%
      dark_layout(
        xaxis = list(title = paste("Genotype for", snp)),
        yaxis = list(title = "Yield (kg)"),
        title = list(text = paste("QTL Boxplot:", snp), font = list(size = 14))
      )
  })
}


# =============================================================================
# RUN APP
# =============================================================================
shinyApp(ui = ui, server = server)
