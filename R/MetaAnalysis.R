# Copyright 2022 Observational Health Data Sciences and Informatics
#
# This file is part of RcaCovidVaccine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# resultsFolder <- "../Results"
# maxCores = 3
# method = "SCCS"
# library(dplyr)
# source("~/git/epi_1071/Code/R/MetaAnalysis.R")

#' Title
#'
#' @param resultsFolder   Folder containing the per-database results, and where
#'                        The meta-analytic estimates will be written. Each DB
#'                        is assumed to have its own sub-folder.
#' @param maxCores        Maximum number of cores used for parallel processing.
#'
#' @return
#' Does not return anything. Produces CSV files in the results folder.
#' 
#' @export
performMetaAnalysis <- function(resultsFolder, maxCores = 8) {
  performMetaAnalysisForMethod(resultsFolder = resultsFolder,
                               method = "CohortMethod",
                               maxCores = maxCores)
  performMetaAnalysisForMethod(resultsFolder = resultsFolder,
                               method = "SCCS",
                               maxCores = maxCores)
}

performMetaAnalysisForMethod <- function(resultsFolder, method, maxCores) {
  message(sprintf("Processing results for method '%s'", method))
  maFolder <- file.path(resultsFolder, "metaAnalysis")
  if (!file.exists(maFolder)) {
    dir.create(maFolder)
  }
  if (method == "CohortMethod") {
    pattern <- "cmResults"
  } else {
    pattern <- "sccsResults"
  }
  folders <- list.files(resultsFolder, pattern = pattern, include.dirs = TRUE, recursive = TRUE)
  loadProfiles <- function(folder) {
    profiles <- readRDS(file.path(resultsFolder, folder, "likelihoodProfiles.rds"))
    return(profiles)
  }
  allProfiles <- lapply(folders, loadProfiles)
  allProfileKeys <- unique(c(unlist(sapply(allProfiles, names))))
  
  loadEstimates <- function(folder) {
    estimates <- readr::read_csv(file.path(resultsFolder, folder, "rawEstimates.csv"), show_col_types = FALSE) 
    if (method == "CohortMethod") {
      estimates$key <- paste(estimates$analysisId, estimates$targetId, estimates$comparatorId, estimates$outcomeId)
    } else {
      estimates$key <- paste(estimates$analysisId, estimates$exposureId, estimates$outcomeId)
    }
    return(estimates)
  }
  allEstimates <- lapply(folders, loadEstimates)
  allEstimateKeys <- unique(c(unlist(sapply(allEstimates, pull, var = "key"))))
  keys <- unique(c(allProfileKeys, allEstimateKeys)) 
  
  pathToCsv <- system.file("settings", "NegativeControls.csv", package = "RcaCovidVaccine")
  negativeControls <- readr::read_csv(pathToCsv, col_types = readr::cols())
  
  # Perform meta-analysis -------------------------
  message("Performing meta-analysis")
  ParallelLogger::addDefaultFileLogger(file.path(maFolder, "log.txt"))
  ParallelLogger::addDefaultErrorReportLogger(file.path(maFolder, "errorReportR.txt"))
  on.exit(ParallelLogger::unregisterLogger("DEFAULT_FILE_LOGGER", silent = TRUE))
  on.exit(ParallelLogger::unregisterLogger("DEFAULT_ERRORREPORT_LOGGER", silent = TRUE), add = TRUE)
  
  cluster <- ParallelLogger::makeCluster(maxCores)
  on.exit(ParallelLogger::stopCluster(cluster), add = TRUE)
  ParallelLogger::clusterRequire(cluster, "dplyr")
  ParallelLogger::clusterApply(cluster, seq_len(maxCores), setThreadData, allProfiles = allProfiles, allEstimates = allEstimates)
  # rm(allProfiles)
  # rm(allEstimates)
  maEstimates <- ParallelLogger::clusterApply(cluster, keys, computeSingleMetaAnalysis, method = method)
  maEstimates <- bind_rows(maEstimates)
  
  if (method == "CohortMethod") {
    # Using negative controls for one outcome lookback window for all outcome lookback windows:
    maEstimates <- bind_rows(
      maEstimates,
      maEstimates %>%
        filter(.data$analysisId <= 19 & .data$outcomeId %in% negativeControls$outcomeId) %>%
        mutate(analysisId = .data$analysisId + 4),
      maEstimates %>%
        filter(.data$analysisId <= 19 & .data$outcomeId %in% negativeControls$outcomeId) %>%
        mutate(analysisId = .data$analysisId + 8),
      maEstimates %>%
        filter(.data$analysisId <= 19 & .data$outcomeId %in% negativeControls$outcomeId) %>%
        mutate(analysisId = .data$analysisId + 12)
    )
    maEstimates$key <- paste(maEstimates$analysisId, maEstimates$targetId, maEstimates$comparatorId, maEstimates$outcomeId)
  }
  
  
  # Perform calibration ---------------------------
  message("Calibrating meta-analysis")
  maEstimates$key <- gsub(" [0-9]+$", "", maEstimates$key)
  maEstimates$nc <- maEstimates$outcomeId %in% negativeControls$outcomeId
  subsets <- split(maEstimates, maEstimates$key)
  
  maEstimates <- ParallelLogger::clusterApply(cluster, subsets, doMaCalibration, resultsFolder)
  maEstimates <- bind_rows(maEstimates)
  
  ncEstimates <- maEstimates %>%
    filter(.data$nc)

  hoiEstimates <- maEstimates %>%
    filter(!.data$nc) %>%
    select(-.data$nc)
  
  if (method == "CohortMethod") {
    readr::write_excel_csv(ncEstimates, file.path(maFolder, "cmNegativeControlEstimates.csv"))
    readr::write_excel_csv(hoiEstimates, file.path(maFolder, "cmEstimates.csv"))
  } else {
    readr::write_excel_csv(ncEstimates, file.path(maFolder, "sccsNegativeControlEstimates.csv"))
    readr::write_excel_csv(hoiEstimates, file.path(maFolder, "sccsEstimates.csv"))
  }
}

setThreadData <- function(index, allProfiles, allEstimates) {
  ParallelLogger::logTrace(sprintf("Settings data for thread %d", index))
  
  allProfiles <<- allProfiles
  allEstimates <<- allEstimates
  return(NULL)
}

# key = keys[1]
# key = "101 5374 137977"
computeSingleMetaAnalysis <- function(key, method) {
  ParallelLogger::logTrace(sprintf("Key: %s", key))
  getProfile <- function(profiles, key) {
    profile <- profiles[[key]]
    return(profile)
  }
  allNa <- function(profile) {
    return(all(is.na(profile$value)))
  }
  profilesSubset <- lapply(allProfiles, getProfile, key)  
  profilesSubset <- profilesSubset[!sapply(profilesSubset, is.null)]
  profilesSubset <- profilesSubset[!sapply(profilesSubset, allNa)]
  
  getEstimate <- function(estimates, key) {
    estimate <- estimates %>%
      filter(.data$key == !!key)
    return(estimate)
  }
  estimatesSubset <- lapply(allEstimates, getEstimate, key)  
  estimatesSubset <- estimatesSubset[!sapply(estimatesSubset, is.null)]
  estimatesSubset <- bind_rows(estimatesSubset)
  
  if (method == "CohortMethod") {
    maEstimate <- estimatesSubset %>%
      group_by(.data$analysisId,
               .data$analysisDescription,
               .data$targetId,
               .data$targetName,
               .data$comparatorId,
               .data$comparatorName,
               .data$outcomeId,
               .data$outcomeName,
               .data$key) %>%
      summarise(target = sum(.data$target),
                comparator = sum(.data$comparator),
                eventsTarget = sum(.data$eventsTarget),
                eventsComparator = sum(.data$eventsComparator),
                .groups = "drop")
  } else {
    maEstimate <- estimatesSubset %>%
      group_by(.data$analysisId,
               .data$analysisDescription,
               .data$exposureId,
               .data$exposureName,
               .data$outcomeId,
               .data$outcomeName,
               .data$key) %>%
      summarise(outcomeSubjects = sum(.data$outcomeSubjects),
                outcomeEvents = sum(.data$outcomeEvents),
                .groups = "drop")
  }
  
  if (length(profilesSubset) <= 1 && nrow(estimatesSubset) > 0) {
    maEstimate$logRrBayesian <- estimatesSubset$logRr[1]
    maEstimate$seLogRrBayesian <- estimatesSubset$seLogRr[1]
    maEstimate$ci95LbBayesian <- estimatesSubset$ci95lb[1]
    maEstimate$ci95UbBayesian <- estimatesSubset$ci95ub[1]
    maEstimate$nSitesBayesian <- 1
    maEstimate$tauBayesian <- NA
  } else {
    meta <- EvidenceSynthesis::computeBayesianMetaAnalysis(profilesSubset)
    maEstimate$logRrBayesian <- meta$mu
    maEstimate$seLogRrBayesian <- meta$muSe
    maEstimate$ci95LbBayesian <- exp(meta$mu95Lb)
    maEstimate$ci95UbBayesian <- exp(meta$mu95Ub)
    maEstimate$nSitesBayesian <- length(profilesSubset) 
    maEstimate$tauBayesian <- meta$tau
    ParallelLogger::logTrace("Finished Bayesian meta-analysis")
  }
  if (nrow(estimatesSubset) == 1) {
    maEstimate$logRrDl <- estimatesSubset$logRr[1]
    maEstimate$seLogRrDl <- estimatesSubset$seLogRr[1]
    maEstimate$ci95LbDl <- estimatesSubset$ci95lb[1]
    maEstimate$ci95UbDl <- estimatesSubset$ci95ub[1]
    maEstimate$nSitesDl <- 1
    maEstimate$i2Dl <- NA
  } else {
    meta <- meta::metagen(TE = estimatesSubset$logRr,
                          seTE = estimatesSubset$seLogRr,
                          sm = "RR",
                          hakn = FALSE,
                          control = list(maxiter=1000))
    s <- summary(meta)
    rnd <- s$random
    maEstimate$logRrDl <- rnd$TE
    maEstimate$seLogRrDl <- rnd$seTE
    maEstimate$ci95LbDl <- exp(rnd$lower)
    maEstimate$ci95UbDl <- exp(rnd$upper)
    maEstimate$nSitesDl <- nrow(estimatesSubset)
    maEstimate$i2Dl <- s$I2
    ParallelLogger::logTrace("Finished non-Bayesian meta-analysis")
  }
  ParallelLogger::logTrace("Finished key")
  return(maEstimate)
}

# subset = subsets[[1]]
# subset <- maEstimates[maEstimates$key == "20 1082 5374", ]
doMaCalibration <- function(subset, resultsFolder) {
  ParallelLogger::logTrace(sprintf("Analysis ID:%s, exposureId: %s", subset$analysisId[1], subset$exposureId))
  ncs <- subset %>%
    filter(.data$nc)
  
  ncsBayesian <- ncs %>%
    filter(!is.na(.data$seLogRrBayesian))
  
  if (nrow(ncsBayesian) >= 5) {
    null <- EmpiricalCalibration::fitMcmcNull(ncsBayesian$logRrBayesian, ncsBayesian$seLogRrBayesian)
    if ("targetId" %in% names(subset)) {
      fileName <- file.path(resultsFolder, "metaAnalysis", sprintf("NcDistribution_a%d_t%d_c%d.png",
                                                                   subset$analysisId[1],
                                                                   subset$targetId[1],
                                                                   subset$comparatorId[1]))
    } else {
      fileName <- file.path(resultsFolder, "metaAnalysis", sprintf("NcDistribution_a%d_e%d.png",
                                                                   subset$analysisId[1],
                                                                   subset$exposureId[1]))
    }
    EmpiricalCalibration::plotCalibrationEffect(logRrNegatives = ncsBayesian$logRrBayesian,
                                                seLogRrNegatives = ncsBayesian$seLogRrBayesian,
                                                null = null,
                                                showCis = TRUE,
                                                title = subset$analysisDescription[1],
                                                fileName = fileName)
    model <- EmpiricalCalibration::convertNullToErrorModel(null)
    calibratedCi <- EmpiricalCalibration::calibrateConfidenceInterval(logRr = subset$logRrBayesian, 
                                                                      seLogRr = subset$seLogRrBayesian, 
                                                                      model = model)
    subset$calibratedRrBayesian <- exp(calibratedCi$logRr)
    subset$calibratedCi95LbBayesian <- exp(calibratedCi$logLb95Rr)
    subset$calibratedCi95UbBayesian <- exp(calibratedCi$logUb95Rr)
  }
  
  ncsDl <- ncs %>%
    filter(!is.na(.data$seLogRrDl))
  
  if (nrow(ncsDl) >= 5) {
    null <- EmpiricalCalibration::fitMcmcNull(ncsDl$logRrDl, ncsDl$seLogRrDl)
    if ("targetId" %in% names(subset)) {
      fileName <- file.path(resultsFolder, "metaAnalysis", sprintf("NcDistribution_a%d_t%d_c%d_DL.png",
                                                                   subset$analysisId[1],
                                                                   subset$targetId[1],
                                                                   subset$comparatorId[1]))
    } else {
      fileName <- file.path(resultsFolder, "metaAnalysis", sprintf("NcDistribution_a%d_e%d_DL.png",
                                                                   subset$analysisId[1],
                                                                   subset$exposureId[1]))
    }
    EmpiricalCalibration::plotCalibrationEffect(logRrNegatives = ncsDl$logRrDl,
                                                seLogRrNegatives = ncsDl$seLogRrDl,
                                                null = null,
                                                showCis = TRUE,
                                                title = subset$analysisDescription[1],
                                                fileName = fileName)
    model <- EmpiricalCalibration::convertNullToErrorModel(null)
    calibratedCi <- EmpiricalCalibration::calibrateConfidenceInterval(logRr = subset$logRrDl, 
                                                                      seLogRr = subset$seLogRrDl, 
                                                                      model = model)
    subset$calibratedRrDl <- exp(calibratedCi$logRr)
    subset$calibratedCi95LbDl <- exp(calibratedCi$logLb95Rr)
    subset$calibratedCi95UbDl <- exp(calibratedCi$logUb95Rr)
  }
  return(subset)
}
