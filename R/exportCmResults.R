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

# source("R/HelperFunctions.R")
# library(dplyr)
exportCmResults <- function(outputFolder, databaseId) {
  cmResultsFolder <- file.path(outputFolder, "cmResults")
  if (!file.exists(cmResultsFolder)) {
    dir.create(cmResultsFolder)
  }
  cmAnalysisListFile <- system.file("settings",
                                    "cmAnalysisList.json",
                                    package = "RcaCovidVaccine"
  )
  cmAnalysisList <- CohortMethod::loadCmAnalysisList(cmAnalysisListFile)
  
  unblind <- readr::read_csv(file.path(outputFolder, "CmDiagnosticsOverview.csv"), col_types = readr::cols())
  unblind <- unblind %>%
    filter(.data$unblind == TRUE) %>%
    select(
      .data$analysisId,
      .data$targetId,
      .data$comparatorId,
      .data$outcomeId
    )
  
  discovery <- getDiscoveryAnalyses() %>%
    filter(.data$method == "CohortMethod") %>%
    select(-.data$method) %>%
    mutate(discovery = TRUE)
  
  criticalValuesFile <- system.file("settings",
                                    "criticalValuesCohortMethod.csv",
                                    package = "RcaCovidVaccine")
  criticalValues <- readr::read_csv(criticalValuesFile, col_types = readr::cols()) %>%
    filter(.data$database == databaseId & .data$method == "Cohort method") %>%
    select(-.data$database, -.data$method)
  
  # Export raw estimates for meta-analysis ---------------------------------------
  cmSummaryFile <- file.path(outputFolder, "analysisSummary.csv")
  cmSummary <- readr::read_csv(cmSummaryFile, col_types = readr::cols())
  
  pathToCsv <- system.file("settings", "NegativeControls.csv", package = "RcaCovidVaccine")
  negativeControls <- readr::read_csv(pathToCsv, col_types = readr::cols())
  
  rawEstimates <- cmSummary %>%
    left_join(unblind %>% mutate(unblind = TRUE), by = c("analysisId", "outcomeId", "targetId", "comparatorId")) %>%
    filter(.data$outcomeId %in% negativeControls$outcomeId | !is.na(.data$unblind)) %>%
    select(-.data$unblind)
  readr::write_excel_csv(rawEstimates, file.path(cmResultsFolder, "rawEstimates.csv"))
  
  # Calibration -------------
  cmSummaryFile <- file.path(outputFolder, "analysisSummary.csv")
  cmSummary <- readr::read_csv(cmSummaryFile, col_types = readr::cols())
  
  profiles <- readRDS(file.path(outputFolder, "cmOutput", "profiles.rds"))
  
  pathToCsv <- system.file("settings", "NegativeControls.csv", package = "RcaCovidVaccine")
  negativeControls <- readr::read_csv(pathToCsv, col_types = readr::cols())
  
  # First compute empirical null distributions
  ncs <- cmSummary %>%
    filter(.data$outcomeId %in% negativeControls$outcomeId)
  
  ncsDatas <- new.env()
  # subset <- split(ncs, paste(ncs$analysisId, ncs$targetId, ncs$comparatorId))[[1]]
  computeNull <- function(subset) {
    profilesSubset <- profiles[paste(subset$analysisId, subset$targetId, subset$comparatorId, subset$outcomeId)]
    profilesSubset <- profilesSubset[!sapply(profilesSubset, is.null)]
    null <- EmpiricalCalibration::fitNullNonNormalLl(profilesSubset)
    
    # null <- EmpiricalCalibration::fitMcmcNull(subset$logRr, subset$seLogRr)
    model <- EmpiricalCalibration::convertNullToErrorModel(null)
    assign(sprintf("t%d_c%d_a%d", subset$targetId[1], subset$comparatorId[1], subset$analysisId[1]),
           list(
             null = null,
             model = model,
             ncs = subset
           ),
           envir = ncsDatas
    )
  }
  plyr::l_ply(split(ncs, paste(ncs$analysisId, ncs$targetId, ncs$comparatorId)), computeNull)
  
  # Calibrate using precomputed nulls
  hois <- cmSummary %>%
    filter(!.data$outcomeId %in% !!negativeControls$outcomeId)
  
  # subset <- split(hois, paste(hois$analysisId, hois$targetId, hois$comparatorId))[[8]]
  performCalibration <- function(subset) {
    # fileName <- file.path(cmResultsFolder, sprintf("NcHoiDistribution_a%d_t%d_c%d.png", subset$analysisId[1], subset$targetId[1], subset$comparatorId[1]))
    # title <- subset$analysisDescription[1]
    
    analysisId <- subset$analysisId[1]
    # Warning: very specific code for current analysis design matrix. Needs updating when analyses are changed:
    while (analysisId > 19 & analysisId < 32) {
      analysisId <- analysisId - 4
    }
    ncsData <- get(sprintf("t%d_c%d_a%d", subset$targetId[1], subset$comparatorId[1], analysisId), envir = ncsDatas)
    
    # Only include unblinded estimates:
    subset <- subset %>%
      inner_join(unblind, by = c("analysisId", "outcomeId", "comparatorId", "targetId"))
    if (nrow(subset) == 0) {
      return(NULL)
    }
    
    null <- ncsData$null
    # EmpiricalCalibration::plotCalibrationEffect(
    #   logRrNegatives = ncsData$ncs$logRr,
    #   seLogRrNegatives = ncsData$ncs$seLogRr,
    #   logRrPositives = subset$logRr,
    #   seLogRrPositives = subset$seLogRr,
    #   null = null,
    #   xLabel = "Hazard Ratio",
    #   title = title,
    #   showCis = FALSE,
    #   fileName = fileName
    # )
    
    results <- subset %>%
      select(-.data$llr) %>%
      left_join(discovery, by = c("analysisId", "outcomeId")) %>%
      mutate(discovery = ifelse(is.na(.data$discovery), FALSE, TRUE)) %>%
      left_join(criticalValues, by = c("analysisId", "outcomeId", "comparatorId", "targetId")) %>%
      select(-.data$totalExpectedEvents)
    
    calibratedP <- EmpiricalCalibration::calibrateP(
      logRr = results$logRr,
      seLogRr = results$seLogRr,
      null = null,
      twoSided = FALSE,
      upper = TRUE
    )
    calibratedEstimates <- EmpiricalCalibration::calibrateConfidenceInterval(
      logRr = results$logRr,
      seLogRr = results$seLogRr,
      model = ncsData$model
    )
    results$calibratedLlr <- NA
    profilesSubset <- profiles[paste(results$analysisId, results$targetId, results$comparatorId, results$outcomeId)]
    for (i in 1:length(profilesSubset)) {
      if (!is.null(profilesSubset[[i]]) && results$discovery[i]) {
        results$calibratedLlr[i] <- EmpiricalCalibration::calibrateLlr(null = null,
                                                                       likelihoodApproximation = profilesSubset[[i]],
                                                                       twoSided = FALSE,
                                                                       upper = TRUE)
      }
    }
    results$calibratedRr <- exp(calibratedEstimates$logRr)
    results$calibratedci95Lb <- exp(calibratedEstimates$logLb95Rr)
    results$calibratedci95Ub <- exp(calibratedEstimates$logUb95Rr)
    results$calibratedP <- calibratedP
    results$calibratedPFromLlr <- EmpiricalCalibration:::computePFromLlr(llr = results$calibratedLlr, results$logRr)
    return(results)
  }
  results <- plyr::llply(split(hois, paste(hois$analysisId, hois$targetId, hois$comparatorId)), performCalibration)
  results <- bind_rows(results)
  
  readr::write_excel_csv(results, file.path(cmResultsFolder, "cmEstimates.csv"))
  
  discoveryResults <- results %>%
    filter(.data$discovery) %>%
    select(.data$analysisId,
           .data$analysisDescription,
           .data$targetId,
           .data$targetName,
           .data$comparatorId,
           .data$comparatorName,
           .data$outcomeId,
           .data$outcomeName,
           .data$calibratedLlr,
           .data$cv) %>%
    mutate(exceeds = .data$calibratedLlr > .data$cv)
  readr::write_excel_csv(discoveryResults, file.path(cmResultsFolder, "cmDiscovery.csv"))
  
  # KM plots ---------------------------
  outcomeIds <- getOutcomesOfInterest()
  subset <- readRDS(file.path(outputFolder, "cmOutput", "outcomeModelReference.rds")) %>%
    filter(.data$outcomeId %in% outcomeIds) %>%
    distinct(.data$analysisId, .data$targetId, .data$comparatorId, .data$outcomeId, .data$strataFile)
  
  # Only include unblinded estimates:
  subset <- subset %>%
    inner_join(unblind, by = c("analysisId", "outcomeId", "comparatorId", "targetId"))
  
  # row <- split(subset, 1:nrow(subset))[[1]]
  plotKm <- function(row) {
    strataPop <- readRDS(file.path(outputFolder, "cmOutput", row$strataFile))
    row <- row %>%
      addCohortNames(IdColumnName = "targetId", nameColumnName = "targetName") %>%
      addCohortNames(IdColumnName = "comparatorId", nameColumnName = "comparatorName") %>%
      addCohortNames(IdColumnName = "outcomeId", nameColumnName = "outcomeName") %>%
      addAnalysisDescription(cmAnalysisList)
    
    fileName <- file.path(cmResultsFolder, sprintf("KaplanMeier_a%s_t%s_c%s_o%s.png", row$analysisId, row$targetId, row$comparatorId, row$outcomeId))
    if (!file.exists(fileName)) {
      title <- sprintf("%s\n%s", row$analysisDescription, gsub(".*\\] ", "", row$outcomeName))
      CohortMethod::plotKaplanMeier(
        population = strataPop,
        targetLabel = gsub(".*exposure to ", "", row$targetName),
        comparatorLabel = gsub(".*exposure to ", "", row$comparatorName),
        title = title,
        fileName = fileName
      )
    }
    # fileName <- file.path(cmResultsFolder, sprintf("Attrition_a%s_o%s.png", row$analysisId, row$outcomeId))
    # if (!file.exists(fileName)) {
    #   CohortMethod::drawAttritionDiagram(object = strataPop,
    #                                      targetLabel = gsub("COVID-19 ", "", gsub(".*exposure to ", "", row$targetName)),
    #                                      comparatorLabel = gsub("\\(Pfizer or Moderna\\) COVID-19 ", "", gsub(".*exposure to ", "", row$comparatorName)),
    #                                      fileName = fileName)
    # }
  }
  plyr::l_ply(split(subset, 1:nrow(subset)), plotKm)
  
  # Profiles ------------------------------------------
  profiles <- readRDS(file.path(outputFolder, "cmOutput", "profiles.rds"))
  cmSummaryFile <- file.path(outputFolder, "analysisSummary.csv")
  cmSummary <- readr::read_csv(cmSummaryFile, col_types = readr::cols())
  pathToCsv <- system.file("settings", "NegativeControls.csv", package = "RcaCovidVaccine")
  negativeControls <- readr::read_csv(pathToCsv, col_types = readr::cols())
  
  toExport <- cmSummary %>%
    left_join(unblind %>% mutate(unblind = TRUE), by = c("analysisId", "outcomeId", "comparatorId", "targetId")) %>%
    mutate(unblind = ifelse(is.na(.data$unblind), FALSE, TRUE)) %>%
    filter(.data$outcomeId %in% negativeControls$outcomeId | .data$unblind)
  
  toExport <- paste(toExport$analysisId, toExport$targetId, toExport$comparatorId, toExport$outcomeId)
  profiles <- profiles[toExport]
  profiles <- profiles[!sapply(profiles, is.null)]
  saveRDS(profiles, file.path(cmResultsFolder, "likelihoodProfiles.rds"))
}
