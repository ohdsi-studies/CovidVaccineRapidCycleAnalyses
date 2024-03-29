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
exportSccsResults <- function(outputFolder, databaseId) {
  sccsResultsFolder <- file.path(outputFolder, "sccsResults")
  if (!file.exists(sccsResultsFolder)) {
    dir.create(sccsResultsFolder)
  }
  sccsSummaryFile <- file.path(outputFolder, "sccsSummary.csv")
  sccsSummary <- readr::read_csv(sccsSummaryFile, col_types = readr::cols())
  outcomeRef <- sccsSummary %>%
    distinct(
      .data$outcomeId,
      .data$outcomeName
    )
  # #these give the not used warning on checkUsagePackage
  # exposureRef <- sccsSummary %>%
  #   distinct(
  #     .data$exposureId,
  #     .data$exposureName
  #   )
  # analysisRef <- sccsSummary %>%
  #   distinct(
  #     .data$analysisId,
  #     .data$analysisDescription
  #   )
  # outcomeIds <- getOutcomesOfInterest()
  # 
  # sccsAnalysisListFile <- system.file("settings", "sccsAnalysisList.json", package = "RcaCovidVaccine")
  # sccsAnalysisList <- SelfControlledCaseSeries::loadSccsAnalysisList(sccsAnalysisListFile)

  unblind <- readr::read_csv(file.path(outputFolder, "sccsDiagnosticsOverview.csv"), col_types = readr::cols())
  unblind <- unblind %>%
    filter(.data$unblind == TRUE) %>%
    select(
      .data$analysisId,
      .data$exposureId,
      .data$outcomeId
    )
  
  discovery <- getDiscoveryAnalyses() %>%
    filter(.data$method == "SCCS") %>%
    select(-.data$method) %>%
    mutate(exposureId = 1082, discovery = TRUE) 
  
  criticalValuesFile <- system.file("settings",
                                    "criticalValuesSccs.csv",
                                    package = "RcaCovidVaccine")
  criticalValues <- readr::read_csv(criticalValuesFile, col_types = readr::cols()) %>%
    filter(.data$database == databaseId & .data$method == "SCCS") %>%
    select(-.data$database, -.data$method)

  profiles <- readRDS(file.path(outputFolder, "sccsOutput", "profiles.rds"))
  
  # Export raw estimates for meta-analysis and signal detection ------------------------
  pathToCsv <- system.file("settings", "NegativeControls.csv", package = "RcaCovidVaccine")
  negativeControls <- readr::read_csv(pathToCsv, col_types = readr::cols())
  
  rawEstimates <- sccsSummary %>%
    left_join(unblind %>% mutate(unblind = TRUE), by = c("analysisId", "outcomeId", "exposureId")) %>%
    filter(.data$outcomeId %in% negativeControls$outcomeId | !is.na(.data$unblind)) %>%
    select(.data$analysisId,
           .data$analysisDescription,
           .data$exposureId,
           .data$exposureName,
           .data$outcomeId,
           .data$outcomeName,
           .data$outcomeSubjects,
           .data$outcomeEvents,
           logRr = .data$`logRr(Main)`,
           seLogRr = .data$`seLogRr(Main)`,
           ci95lb = .data$`ci95lb(Main)`,
           ci95ub = .data$`ci95ub(Main)`)
  readr::write_excel_csv(rawEstimates, file.path(sccsResultsFolder, "rawEstimates.csv"))

  # Empirical calibration --------------------------------
  pathToCsv <- system.file("settings", "NegativeControls.csv", package = "RcaCovidVaccine")
  negativeControls <- readr::read_csv(pathToCsv, col_types = readr::cols())

  sccsSummary$nc <- sccsSummary$outcomeId %in% negativeControls$outcomeId
  
  omr <- readRDS(file.path(outputFolder, "sccsOutput", "outcomeModelReference.rds"))

  # subset <- split(sccsSummary, paste(sccsSummary$analysisId, sccsSummary$exposureId))[[4]]
  # subset <- sccsSummary[sccsSummary$analysisId == 2 & sccsSummary$exposureId == 5374, ]
  performCalibration <- function(subset) {
    # print(paste(subset$analysisId[1], subset$exposureId[1]))
    # fileName <- file.path(
    #   sccsResultsFolder,
    #   sprintf("NcHoiDistribution_a%d_e%d.png", subset$analysisId[1], subset$exposureId[1])
    # )
    # title <- sprintf("%s\n%s", gsub(".*exposure to", "", subset$exposureName), subset$analysisDescription)

    ncs <- subset %>%
      filter(.data$nc)

    # Only include unblinded estimates:
    hois <- subset %>%
      filter(!.data$nc) %>%
      inner_join(unblind, by = c("analysisId", "outcomeId", "exposureId"))
    if (nrow(hois) == 0) {
      return(NULL)
    }
    profilesSubset <- profiles[paste(ncs$analysisId, ncs$exposureId, ncs$outcomeId)]
    profilesSubset <- profilesSubset[!sapply(profilesSubset, is.null)]
    isInvalid <- function(profile) {
      return(any(is.na(profile$point) | is.na(profile$value)))
    }
    invalid <- sapply(profilesSubset, isInvalid)
    profilesSubset <- profilesSubset[!invalid]
    null <- EmpiricalCalibration::fitNullNonNormalLl(profilesSubset)
    # null <- EmpiricalCalibration::fitMcmcNull(logRr = ncs$`logRr(Main)`,
    #                                           seLogRr = ncs$`seLogRr(Main)`)

    # EmpiricalCalibration::plotCalibrationEffect(
    #   logRrNegatives = ncs$`logRr(Main)`,
    #   seLogRrNegatives = ncs$`seLogRr(Main)`,
    #   logRrPositives = hois$`logRr(Main)`,
    #   seLogRrPositives = hois$`seLogRr(Main)`,
    #   null = null,
    #   xLabel = "Incidence Rate Ratio"
    # )

    model <- EmpiricalCalibration::convertNullToErrorModel(null)
    omrSubset <- omr[omr$exposureId == subset$exposureId[1] & omr$outcomeId %in% hois$outcomeId & omr$analysisId == subset$analysisId[1], ]
    # row <- split(omrSubset, 1:nrow(omrSubset))[[1]]
    getCounts <- function(row) {
      outcomeModel <- readRDS(file.path(outputFolder, "sccsOutput", row$sccsModelFile))
      stats <- outcomeModel$metaData$covariateStatistics[outcomeModel$metaData$covariateStatistics$covariateId == 1000, ]
      return(tibble(
        outcomeId = row$outcomeId,
        outcomesExposedDuringTar = stats$outcomeCount,
        outcomesExposedOutsideTar = stats$personCount - stats$outcomeCount
      ))
    }
    counts <- lapply(split(omrSubset, 1:nrow(omrSubset)), getCounts)
    counts <- do.call(rbind, counts)

    hois$p <- mapply(function(logRr, seLogRr) {
      EmpiricalCalibration::computeTraditionalP(logRr, seLogRr)
    }, logRr = hois$`logRr(Main)`, seLogRr = hois$`seLogRr(Main)`)
    results <- hois %>%
      inner_join(counts, by = "outcomeId") %>%
      select(.data$analysisId,
        .data$analysisDescription,
        .data$outcomeId,
        .data$outcomeName,
        .data$exposureId,
        .data$exposureName,
        .data$outcomesExposedDuringTar,
        .data$outcomesExposedOutsideTar,
        rr = .data$`rr(Main)`,
        ci95Lb = .data$`ci95lb(Main)`,
        ci95Ub = .data$`ci95ub(Main)`,
        p = .data$p,
        logRr = .data$`logRr(Main)`,
        seLogRr = .data$`seLogRr(Main)`
      ) %>%
      mutate(
        rr = if_else(is.na(.data$seLogRr), as.numeric(NA), .data$rr),
        logRr = if_else(is.na(.data$seLogRr), as.numeric(NA), .data$logRr),
        ci95Lb = if_else(is.na(.data$seLogRr), as.numeric(NA), .data$ci95Lb),
        ci95Ub = if_else(is.na(.data$seLogRr), as.numeric(NA), .data$ci95Ub)
      ) %>%
      left_join(discovery, by = c("analysisId", "outcomeId", "exposureId")) %>%
      mutate(discovery = ifelse(is.na(.data$discovery), FALSE, TRUE)) %>%
      left_join(criticalValues, by = c("analysisId", "outcomeId", "exposureId")) %>%
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
      model = model
    )
    results$calibratedLlr <- NA
    profilesSubset <- profiles[paste(results$analysisId, results$exposureId, results$outcomeId)]
    for (i in 1:length(profilesSubset)) {
      if (!is.null(profilesSubset[[i]]) && results$discovery[i])
        results$calibratedLlr[i] <- EmpiricalCalibration::calibrateLlr(null = null,
                                                                       likelihoodApproximation = profilesSubset[[i]],
                                                                       twoSided = FALSE,
                                                                       upper = TRUE)
    }
    results$calibratedRr <- exp(calibratedEstimates$logRr)
    results$calibratedci95Lb <- exp(calibratedEstimates$logLb95Rr)
    results$calibratedci95Ub <- exp(calibratedEstimates$logUb95Rr)
    results$calibratedP <- calibratedP
    results$calibratedPFromLlr <- EmpiricalCalibration:::computePFromLlr(llr = results$calibratedLlr, results$logRr)
    return(results)
  }
  results <- plyr::llply(split(sccsSummary, paste(sccsSummary$analysisId, sccsSummary$exposureId)), performCalibration)
  results <- do.call(rbind, results)

  readr::write_excel_csv(results, file.path(sccsResultsFolder, "sccsEstimates.csv"))
  
  discoveryResults <- results %>%
    filter(.data$discovery) %>%
    select(.data$analysisId,
           .data$analysisDescription,
           .data$exposureId,
           .data$exposureName,
           .data$outcomeId,
           .data$outcomeName,
           .data$calibratedLlr,
           .data$cv) %>%
    mutate(exceeds = .data$calibratedLlr > .data$cv)
  readr::write_excel_csv(discoveryResults, file.path(sccsResultsFolder, "sccsDiscovery.csv"))
  
  # Profiles ------------------------------------------
  profiles <- readRDS(file.path(outputFolder, "sccsOutput", "profiles.rds"))
  sccsSummaryFile <- file.path(outputFolder, "sccsSummary.csv")
  sccsSummary <- readr::read_csv(sccsSummaryFile, col_types = readr::cols())
  pathToCsv <- system.file("settings", "NegativeControls.csv", package = "RcaCovidVaccine")
  negativeControls <- readr::read_csv(pathToCsv, col_types = readr::cols())
  
  toExport <- sccsSummary %>%
    left_join(unblind %>% mutate(unblind = TRUE), by = c("analysisId", "exposureId", "outcomeId")) %>%
    mutate(unblind = ifelse(is.na(.data$unblind), FALSE, TRUE)) %>%
    filter(.data$outcomeId %in% negativeControls$outcomeId | .data$unblind)
  
  toExport <- paste(toExport$analysisId, toExport$exposureId, toExport$outcomeId)
  
  profiles <- profiles[toExport]
  profiles <- profiles[!sapply(profiles, is.null)]
  saveRDS(profiles, file.path(sccsResultsFolder, "likelihoodProfiles.rds"))
}
