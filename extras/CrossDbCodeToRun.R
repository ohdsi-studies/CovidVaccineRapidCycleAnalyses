resultsFolder <- "../Results"
maxCores <- 8

# Merge cohort diagnostics ---------------------------------
shinyFolder <- "c:/temp/RcaCovidVaccineCohortDiagnostics"
dir.create(shinyFolder)
shinyDataFolder <- file.path(shinyFolder, "data")
dir.create(shinyDataFolder)

tempFolder <- tempfile("cdMergeTempFolder")
dir.create(tempFolder)

files <- list.files(resultsFolder, pattern = "Results_.*.zip", full.names = TRUE, recursive = TRUE)
file.copy(files, tempFolder)
CohortDiagnostics::preMergeDiagnosticsFiles(tempFolder)

file.copy(file.path(tempFolder, "PreMerged.RData"), shinyDataFolder)
unlink(tempFolder, recursive = TRUE)

shinyAppFolder <- system.file("shiny", "DiagnosticsExplorer", package = "CohortDiagnostics")
files <- list.files(shinyAppFolder, include.dirs = TRUE, full.names = TRUE)
file.copy(files, shinyFolder, recursive = TRUE)


# Perform meta-analysis ----------------------------------------------------
library(RcaCovidVaccine)
performMetaAnalysis(resultsFolder, maxCores)


# Combine discovery across databases and methods -----------------------------
library(dplyr)
loadFile <- function(file) {
  data <- readr::read_csv(file.path(resultsFolder, file), col_types = readr::cols())
  data$database <- dirname(dirname(file))
  return(data)
}

files <- list.files(resultsFolder, pattern = "cmDiscovery.csv", include.dirs = TRUE, recursive = TRUE)
cmDiscovery <- lapply(files, loadFile)
cmDiscovery <- bind_rows(cmDiscovery)
cmDiscovery <- cmDiscovery %>%
  rename(exposureId = targetId,
         exposureName = targetName) %>%
  mutate(method = "Cohort method")

files <- list.files(resultsFolder, pattern = "sccsDiscovery.csv", include.dirs = TRUE, recursive = TRUE)
sccsDiscovery <- lapply(files, loadFile)
sccsDiscovery <- bind_rows(sccsDiscovery)
sccsDiscovery <- sccsDiscovery %>%
  mutate(method = "SCCS")

discovery <- bind_rows(cmDiscovery, sccsDiscovery) %>%
  arrange(outcomeName, method, analysisId, database)
readr::write_excel_csv(discovery, file.path(resultsFolder, "Discovery.csv"))

discoverySummary <- discovery %>%
  group_by(outcomeId, outcomeName) %>%
  summarise(exceeds = max(exceeds), .groups = "drop") %>%
  mutate(exceeds = exceeds == 1)
readr::write_excel_csv(discoverySummary, file.path(resultsFolder, "DiscoverySummary.csv"))

# Combine estimation across databases and methods -----------------------------
library(dplyr)
loadFile <- function(file) {
  data <- readr::read_csv(file.path(resultsFolder, file), col_types = readr::cols())
  data$database <- dirname(dirname(file))
  return(data)
}

files <- list.files(resultsFolder, pattern = "cmEstimates.csv", include.dirs = TRUE, recursive = TRUE)
cmEstimation <- lapply(files, loadFile)
cmEstimation <- bind_rows(cmEstimation)

effectiveAlpha <- sum(cmEstimation$cvAlpha, na.rm = TRUE)

cmEstimation <- cmEstimation %>%
  rename(exposureId = targetId,
         exposureName = targetName) %>%
  mutate(database = ifelse(database == ".", "Meta-analysis", database)) %>%
  mutate(method = "Cohort method",
         rr = ifelse(database == "Meta-analysis", calibratedRrBayesian, calibratedRr),
         ci95Lb = ifelse(database == "Meta-analysis", calibratedCi95LbBayesian, calibratedci95Lb),
         ci95Ub = ifelse(database == "Meta-analysis", calibratedCi95UbBayesian, calibratedci95Ub),
         discovery = ifelse(is.na(discovery), FALSE, discovery),
         exceeds = calibratedLlr > cv) %>%
  select(database,
         method,
         analysisId,
         analysisDescription,
         exposureId,
         exposureName,
         comparatorId,
         comparatorName,
         outcomeId,
         outcomeName,
         eventsTarget,
         eventsComparator,
         rr,
         ci95Lb,
         ci95Ub,
         tau = tauBayesian,
         discovery,
         exceeds)

files <- list.files(resultsFolder, pattern = "sccsEstimates.csv", include.dirs = TRUE, recursive = TRUE)
sccsEstimation <- lapply(files, loadFile)
sccsEstimation <- bind_rows(sccsEstimation)

# effectiveAlpha <- effectiveAlpha + sum(cmEstimation$cvAlpha, na.rm = TRUE)

sccsEstimation <- sccsEstimation %>%
  mutate(database = ifelse(database == ".", "Meta-analysis", database)) %>%
  mutate(method = "SCCS",
         rr = ifelse(database == "Meta-analysis", calibratedRrBayesian, calibratedRr),
         ci95Lb = ifelse(database == "Meta-analysis", calibratedCi95LbBayesian, calibratedci95Lb),
         ci95Ub = ifelse(database == "Meta-analysis", calibratedCi95UbBayesian, calibratedci95Ub),
         discovery = ifelse(is.na(discovery), FALSE, discovery),
         exceeds = calibratedLlr > cv) %>%
  select(database,
         method,
         analysisId,
         analysisDescription,
         exposureId,
         exposureName,
         outcomeId,
         outcomeName,
         eventsTarget = outcomesExposedDuringTar,
         eventsComparator = outcomesExposedOutsideTar,
         rr,
         ci95Lb,
         ci95Ub,
         tau = tauBayesian,
         discovery,
         exceeds)

estimation <- bind_rows(cmEstimation, sccsEstimation) %>%
  arrange(outcomeName, method, analysisId, database)
readr::write_excel_csv(estimation, file.path(resultsFolder, "Estimation.csv"))


# writeLines(sprintf("Effective alpha: %0.3f", effectiveAlpha))

# Forest plot per outcome -----------------------------------------
library(ggplot2)
estimation <- readr::read_csv(file.path(resultsFolder, "Estimation.csv"), show_col_types = FALSE)
estimation <- estimation[estimation$exposureId == 1082, ]
estimation$tar <- paste("TAR:", gsub(",.*$", "", gsub("^.*TAR = ", "", estimation$analysisDescription)))
estimation$analysis <- estimation$method
estimation$analysis[grepl("post-exposure", estimation$analysisDescription)] <- "SCCS, post-exposure only"
estimation$analysis[grepl("longer lookback", estimation$analysisDescription)] <- "SCCS, longer lookback"
estimation$analysis[!is.na(estimation$comparatorId) & estimation$comparatorId == 5373] <- "Cohort method, vs Pfizer"
estimation$analysis[!is.na(estimation$comparatorId) & estimation$comparatorId == 5374] <- "Cohort method, vs Moderna"
estimation$database <- factor(estimation$database, levels = c("CCAE", "Optum", "Pharmetrics", "Meta-analysis"))


subset <- estimation[estimation$outcomeId == 609, ]
createPlotForOutcome <- function(subset) {
  fileName <- file.path(resultsFolder, "metaAnalysis", sprintf("Forest_o%d.png", subset$outcomeId[1]))    
  if (!file.exists(fileName)) {
    breaks <- c(0.1, 0.25, 0.5, 1, 2, 4, 6, 8)
    plot <- ggplot(subset, aes(x = rr, y = analysis)) +
      geom_vline(xintercept = 1, size = 0.5) +
      geom_point(color = "#000088", alpha = 0.8) +
      geom_errorbarh(aes(xmin = ci95Lb, xmax = ci95Ub), height = 0.5, color = "#000088", alpha = 0.8) +
      scale_x_log10("Effect size (Hazard Ratio or Incidence Rate Ratio)", breaks = breaks, labels = breaks) +
      coord_cartesian(xlim = c(0.1, 10)) + 
      facet_grid(database ~ tar)  +
      ggtitle(subset$outcomeName[1]) +
      theme(axis.title.y = element_blank(),
            panel.grid.minor = element_blank())
    
    ggsave(filename = fileName, plot = plot, width = 9, height = 5, dpi = 300)
  }
  return(NULL)
}
lapply(split(estimation, estimation$outcomeId), createPlotForOutcome)

# Exposure counts across databases ---------------------------------------------
library(dplyr)
files <- list.files(resultsFolder, pattern = "CohortCounts.csv", recursive = TRUE)
tcosOfInterest <- readr::read_csv("inst/settings/TcosOfInterest.csv", show_col_types = FALSE)
exposureIds <- unique(c(tcosOfInterest$targetId, tcosOfInterest$comparatorId))

# file = files[1]
getExposureCounts <- function(file) {
  database <- gsub("/CohortCounts.csv", "", file)
  counts <- readr::read_csv(file.path(resultsFolder, file), show_col_types = FALSE)
  counts %>% 
    filter(cohortDefinitionId %in% exposureIds) %>%
    select(-cohortCount) %>%
    mutate(database = !!database,
           cohortName = gsub("^.*exposure to ", "", cohortName),
           prettyCount = format(personCount, big.mark = ",", scientific = FALSE)) %>%
    return()
}
counts <- lapply(files, getExposureCounts)
counts <- bind_rows(counts)
readr::write_csv(counts, file.path(resultsFolder, "ExposureCounts.csv"))

# Count of estimates passing diagnostics ---------------------------------------
library(dplyr)

loadFile <- function(file) {
  print(file)
  data <- readr::read_csv(file.path(resultsFolder, file), col_types = readr::cols()) #%>%
    #SqlRender::snakeCaseToCamelCaseNames()
  data$database <- dirname(dirname(file))
  return(data)
}
files <- list.files(resultsFolder, pattern = "CmDiagnosticsOverview.csv", recursive = TRUE)
diagSumCm <- lapply(files, loadFile)
diagSumCm <- bind_rows(diagSumCm)

files <- list.files(resultsFolder, pattern = "SccsDiagnosticsOverview.csv", include.dirs = TRUE, recursive = TRUE)
diagSumSccs <- lapply(files, loadFile)
diagSumSccs <- bind_rows(diagSumSccs)

diagSumCm %>%
  summarise(estimates = n(), passingAll = sum(unblind), fraction = mean(unblind))

diagSumSccs %>%
  summarise(estimates = n(), passingAll = sum(unblind, na.rm = TRUE), fraction = mean(unblind, na.rm = TRUE))


bind_rows(
  diagSumCm %>% 
    select(unblind),
  diagSumSccs %>% 
    select(unblind)
) %>%
  summarise(estimates = n(), passingAll = sum(unblind, na.rm = TRUE), fraction = mean(unblind, na.rm = TRUE))



