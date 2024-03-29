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
# library(survival)
generateSccsDiagnostics <- function(outputFolder, maxCores) {
  mainExposureId <- 1082
  sccsDiagnosticsFolder <- file.path(outputFolder, "sccsDiagnostics")
  if (!file.exists(sccsDiagnosticsFolder)) {
    dir.create(sccsDiagnosticsFolder)
  }
  sccsSummaryFile <- file.path(outputFolder, "sccsSummary.csv")
  sccsSummary <- readr::read_csv(sccsSummaryFile, col_types = readr::cols())
  outcomeRef <- sccsSummary %>%
    distinct(
      .data$outcomeId,
      .data$outcomeName
    )
  exposureRef <- sccsSummary %>%
    distinct(
      .data$exposureId,
      .data$exposureName
    )
  analysisRef <- sccsSummary %>%
    distinct(
      .data$analysisId,
      .data$analysisDescription
    )
  outcomeIds <- getOutcomesOfInterest()
  
  sccsAnalysisListFile <- system.file("settings", "sccsAnalysisList.json", package = "RcaCovidVaccine")
  sccsAnalysisList <- SelfControlledCaseSeries::loadSccsAnalysisList(sccsAnalysisListFile)
  
  # Negative control plots --------------------------------
  pathToCsv <- system.file("settings", "NegativeControls.csv", package = "RcaCovidVaccine")
  negativeControls <- readr::read_csv(pathToCsv, col_types = readr::cols())
  
  sccsSummary <- sccsSummary %>%
    inner_join(negativeControls, by = "outcomeId")
  
  # subset <- split(sccsSummary, paste(sccsSummary$analysisId, sccsSummary$exposureId))[[5]]
  plotCalibration <- function(subset) {
    fileName <- file.path(
      sccsDiagnosticsFolder,
      sprintf("NcDistribution_a%d_e%d.png", subset$analysisId[1], subset$exposureId[1])
    )
    if (!file.exists(fileName)) {
      title <- sprintf("%s\n%s", gsub(".*exposure to", "", subset$exposureName), subset$analysisDescription)
      
      # Assume likelihood normally distributed
      EmpiricalCalibration::plotCalibrationEffect(
        logRrNegatives = subset$`logRr(Main)`,
        seLogRrNegatives = subset$`seLogRr(Main)`,
        xLabel = "Incidence Rate Ratio",
        title = title,
        showCis = TRUE,
        fileName = fileName
      )
    }
    null <- EmpiricalCalibration::fitMcmcNull(subset$`logRr(Main)`, subset$`seLogRr(Main)`)
    ease <- EmpiricalCalibration::computeExpectedAbsoluteSystematicError(null)
    return(data.frame(
      exposureId = subset$exposureId[1],
      analysisId = subset$analysisId[1],
      ease = ease$ease
    ))
  }
  eases <- plyr::llply(split(sccsSummary, paste(sccsSummary$analysisId, sccsSummary$exposureId)), plotCalibration)
  eases <- bind_rows(eases)
  eases <- eases %>%
    mutate(easeDiagnostic = case_when(
      abs(.data$ease) < 0.1 ~ "PASS",
      abs(.data$ease) < 0.25 ~ "WARNING",
      TRUE ~ "FAIL"
    ))
  
  # Event dependent observation end and calendar time ------------------
  subset <- readRDS(file.path(outputFolder, "sccsOutput", "outcomeModelReference.rds")) %>%
    filter(.data$outcomeId %in% outcomeIds & .data$exposureId == mainExposureId) %>%
    distinct(.data$analysisId, .data$outcomeId, .data$studyPopFile, .data$sccsModelFile) %>%
    inner_join(outcomeRef, by = "outcomeId") %>%
    inner_join(analysisRef, by = "analysisId")
  
  # which(subset$analysisId == 2 & subset$outcomeId == 5831)
  # row <- split(subset, 1:nrow(subset))[[12]]
  createPlots <- function(row, outputFolder, sccsDiagnosticsFolder) {
    # print(row)
    studyPop <- readRDS(file.path(file.path(outputFolder, "sccsOutput", row$studyPopFile)))
    if (nrow(studyPop$outcomes) ==0) {
      return(data.frame(
        outcomeId = row$outcomeId,
        analysisId = row$analysisId,
        stable = FALSE
      ))
    }
    model <- readRDS(file.path(file.path(outputFolder, "sccsOutput", row$sccsModelFile)))
    
    fileName <- file.path(
      sccsDiagnosticsFolder,
      sprintf("EventObservationDependence_a%d_o%s.png", row$analysisId, row$outcomeId)
    )
    if (!file.exists(fileName)) {
      title <- sprintf("%s\n%s", gsub(".*\\] ", "", row$outcomeName), row$analysisDescription)
      SelfControlledCaseSeries::plotEventObservationDependence(
        studyPopulation = studyPop,
        title = title,
        fileName = fileName
      )
    }
    fileName <- file.path(
      sccsDiagnosticsFolder,
      sprintf("EventToCalendarTime_a%d_o%s.png", row$analysisId, row$outcomeId)
    )
    if (!file.exists(fileName)) {
      title <- sprintf("%s\n%s", gsub(".*\\] ", "", row$outcomeName), row$analysisDescription)
      SelfControlledCaseSeries::plotEventToCalendarTime(
        studyPopulation = studyPop,
        sccsModel = model,
        title = title,
        fileName = fileName
      )
    }
    
    # outcomes <- studyPop$outcomes %>%
    #   group_by(.data$caseId) %>%
    #   summarise(outcomeDay = min(.data$outcomeDay), .groups = "drop_last") %>%
    #   inner_join(studyPop$cases, by = "caseId") %>%
    #   transmute(
    #     daysFromEvent = .data$endDay - .data$outcomeDay,
    #     censoring = case_when(
    #       .data$noninformativeEndCensor == 1 ~ "Uncensored",
    #       TRUE ~ "Censored"
    #     )
    #   )
    # 
    # medianDaysToEnd <- outcomes %>%
    #   group_by(.data$censoring) %>%
    #   summarise(medianDaysToEnd = median(.data$daysFromEvent))
    # Wilcoxon Rank Sum Test always significant due to large sample, so not informative:
    # wilcox.test(daysFromEvent ~ censoring, data = outcomes)
    
    timeStability <- SelfControlledCaseSeries::computeTimeStability(
      studyPopulation = studyPop,
      sccsModel = model,
      maxRatio = 1.25,
      alpha = 0.05
    )
    return(data.frame(
      outcomeId = row$outcomeId,
      analysisId = row$analysisId,
      # medianDaysToEndCensored = medianDaysToEnd$medianDaysToEnd[medianDaysToEnd$censoring == "Censored"],
      # medianDaysToEndUncensored = medianDaysToEnd$medianDaysToEnd[medianDaysToEnd$censoring == "Uncensored"],
      stable = all(timeStability$stable)
    ))
  }
  cluster <- ParallelLogger::makeCluster(min(10, maxCores))
  medianDaysToEnd <- ParallelLogger::clusterApply(cluster, split(subset, 1:nrow(subset)), createPlots, outputFolder = outputFolder, sccsDiagnosticsFolder = sccsDiagnosticsFolder)
  ParallelLogger::stopCluster(cluster)
  
  # medianDaysToEnd <- plyr::llply(split(subset, 1:nrow(subset)), createPlots, outputFolder = outputFolder, sccsDiagnosticsFolder = sccsDiagnosticsFolder)
  medianDaysToEnd <- bind_rows(medianDaysToEnd)
  medianDaysToEnd <- medianDaysToEnd %>%
    mutate(
      # medianDaysToEndDiagnostic = case_when(
      #   .data$medianDaysToEndCensored / .data$medianDaysToEndUncensored > 0.9 ~ "PASS",
      #   .data$medianDaysToEndCensored / .data$medianDaysToEndUncensored > 0.5 ~ "WARNING",
      #   TRUE ~ "FAIL"
      # ),
      timeTrendDiagnostic = case_when(
        .data$stable ~ "PASS",
        TRUE ~ "FAIL"
      )
    )
  
  # Exposure-centered plots ---------------------------
  # Pick one analysis ID per outcome:
  analyses <- subset %>%
    filter(.data$outcomeId %in% outcomeIds) %>%
    group_by(.data$outcomeId) %>%
    summarise(analysisId = min(analysisId))
  
  subset <- readRDS(file.path(outputFolder, "sccsOutput", "outcomeModelReference.rds")) %>%
    inner_join(analyses, by = c("outcomeId", "analysisId")) %>%
    distinct(.data$exposureId, .data$outcomeId, .data$studyPopFile, .data$sccsDataFile) %>%
    inner_join(outcomeRef, by = "outcomeId") %>%
    inner_join(exposureRef, by = "exposureId")
  
  # row <- split(subset, 1:nrow(subset))[[1]]
  createCenteredPlots <- function(row, outputFolder, sccsDiagnosticsFolder) {
    # print(row)
    
    computePreExposureGainP <- function(sccsData, studyPopulation, exposureEraId) {
      # Copied from plotExposureCentered:
      cases <- studyPopulation$cases %>%
        select(.data$caseId, caseEndDay = .data$endDay, .data$offset)
      
      exposures <- sccsData$eras %>%
        filter(.data$eraId == exposureEraId & .data$eraType == "rx") %>%
        group_by(.data$caseId) %>%
        inner_join(cases, by = "caseId", copy = TRUE) %>%
        mutate(
          startDay = .data$startDay - .data$offset,
          endDay = .data$endDay - .data$offset
        ) %>%
        filter(.data$startDay >= 0, .data$startDay < .data$caseEndDay) %>%
        collect()
      
      if (nrow(exposures) == 0) {
        warning("No exposures found with era ID ", exposureEraId)
        return(NA)
      }
      firstExposures <- exposures %>%
        group_by(.data$caseId, .data$caseEndDay) %>%
        summarise(
          startDay = min(.data$startDay, na.rm = TRUE),
          endDay = min(.data$endDay, na.rm = TRUE),
          .groups = "drop_last"
        )
      
      outcomes <- studyPopulation$outcomes %>%
        inner_join(firstExposures, by = "caseId") %>%
        mutate(delta = .data$outcomeDay - .data$startDay) %>%
        select(.data$caseId, .data$outcomeDay, .data$delta)
      
      # end copy
      
      # Restrict to 30 days before and after exposure start:
      outcomes <- outcomes %>%
        filter(.data$delta >= -30 & .data$delta <= 30) %>%
        mutate(
          beforeExposure = .data$delta < 0,
          y = 1
        ) %>%
        select(
          .data$caseId,
          .data$beforeExposure,
          .data$y
        )
      
      observed <- bind_rows(
        firstExposures %>%
          mutate(
            daysObserved = if_else(.data$startDay > 30, 30, .data$startDay),
            beforeExposure = TRUE
          ) %>%
          select(
            .data$caseId,
            .data$daysObserved,
            .data$beforeExposure
          ),
        firstExposures %>%
          mutate(daysAfterExposure = .data$caseEndDay - .data$startDay) %>%
          mutate(
            daysObserved = if_else(.data$daysAfterExposure > 30, 30, .data$daysAfterExposure),
            beforeExposure = FALSE
          ) %>%
          select(
            .data$caseId,
            .data$daysObserved,
            .data$beforeExposure
          )
      ) %>%
        filter(.data$daysObserved > 0)
      
      poissonData <- observed %>%
        left_join(outcomes, by = c("caseId", "beforeExposure")) %>%
        mutate(
          y = if_else(is.na(.data$y), 0, .data$y),
          logDays = log(.data$daysObserved)
        )
      # library(survival)
      cyclopsData <- Cyclops::createCyclopsData(y ~ beforeExposure + strata(caseId) + offset(logDays),
                                                modelType = "cpr",
                                                data = poissonData
      )
      fit <- Cyclops::fitCyclopsModel(cyclopsData)
      if (fit$return_flag == "ILLCONDITIONED") {
        return(NA)
      }
      # compute one-sided p-value:
      llNull <- Cyclops::getCyclopsProfileLogLikelihood(
        object = fit,
        parm = "beforeExposureTRUE",
        x = 0
      )$value
      llr <- fit$log_likelihood - llNull
      p <- EmpiricalCalibration:::computePFromLlr(llr, coef(fit))
      return(p)
    }
    
    sccsData <- SelfControlledCaseSeries::loadSccsData(file.path(outputFolder, "sccsOutput", row$sccsDataFile))
    studyPopulation <- readRDS(file.path(file.path(outputFolder, "sccsOutput", row$studyPopFile)))
    
    fileName <- file.path(sccsDiagnosticsFolder, sprintf("ExposureCentered_e%s_o%s.png", row$exposureId, row$outcomeId))
    if (!file.exists(fileName)) {
      title <- sprintf(
        "%s\n%s",
        gsub(".*exposure to ", "", row$exposureName),
        gsub(".*\\] ", "", row$outcomeName)
      )
      
      SelfControlledCaseSeries::plotExposureCentered(
        studyPopulation = studyPopulation,
        sccsData = sccsData,
        exposureEraId = row$exposureId,
        highlightExposedEvents = FALSE,
        title = title,
        fileName = fileName
      )
    }
    p <- computePreExposureGainP(
      sccsData = sccsData,
      studyPopulation = studyPopulation,
      exposureEraId = row$exposureId
    )
    return(data.frame(
      exposureId = row$exposureId,
      outcomeId = row$outcomeId,
      preExposureP = p
    ))
  }
  
  cluster <- ParallelLogger::makeCluster(min(10, maxCores))
  ParallelLogger::clusterRequire(cluster, "dplyr")
  ParallelLogger::clusterRequire(cluster, "survival")
  preExposure <- ParallelLogger::clusterApply(cluster, split(subset, 1:nrow(subset)), createCenteredPlots, outputFolder = outputFolder, sccsDiagnosticsFolder = sccsDiagnosticsFolder)
  ParallelLogger::stopCluster(cluster)
  
  # preExposure <- plyr::llply(split(subset, 1:nrow(subset)), createCenteredPlots)
  preExposure <- bind_rows(preExposure)
  preExposure <- preExposure %>%
    mutate(preExposureDiagnostic = case_when(
      .data$preExposureP >= 0.05 ~ "PASS",
      .data$preExposureP >= 0.01 ~ "WARNING",
      TRUE ~ "FAIL"
    ))
  
  # MDRR ---------------------------
  subset <- readRDS(file.path(outputFolder, "sccsOutput", "outcomeModelReference.rds")) %>%
    filter(.data$outcomeId %in% outcomeIds) %>%
    distinct(.data$analysisId, .data$exposureId, .data$outcomeId, .data$sccsIntervalDataFile) %>%
    inner_join(outcomeRef, by = "outcomeId") %>%
    inner_join(exposureRef, by = "exposureId") %>%
    inner_join(analysisRef, by = "analysisId")
  
  # row <- split(subset, 1:nrow(subset))[[1]]
  computeMdrr <- function(row, outputFolder) {
    # print(row)
    sccsIntervalData <- SelfControlledCaseSeries::loadSccsIntervalData(file.path(outputFolder, "sccsOutput", row$sccsIntervalDataFile))
    
    result <- SelfControlledCaseSeries::computeMdrr(
      sccsIntervalData = sccsIntervalData,
      exposureCovariateId = 1000,
      method = "SRL1"
    )
    result %>%
      mutate(
        exposureId = row$exposureId,
        exposureName = row$exposureName,
        outcomeId = row$outcomeId,
        outcomeName = row$outcomeName,
        analysisId = row$analysisId,
        analysisDescription = row$analysisDescription
      ) %>%
      return()
  }
  cluster <- ParallelLogger::makeCluster(min(10, maxCores))
  ParallelLogger::clusterRequire(cluster, "dplyr")
  results <- ParallelLogger::clusterApply(cluster, split(subset, 1:nrow(subset)), computeMdrr, outputFolder = outputFolder)
  ParallelLogger::stopCluster(cluster)
  
  # results <- plyr::llply(split(subset, 1:nrow(subset)), computeMdrr)
  results <- bind_rows(results)
  readr::write_csv(results, file.path(sccsDiagnosticsFolder, "MDRR.csv"))
  
  # Write overview ------------------------
  overview <- results %>%
    mutate(mdrrDiagnostic = case_when(
      .data$mdrr < 2 ~ "PASS",
      .data$mdrr < 10 ~ "WARNING",
      TRUE ~ "FAIL"
    )) %>%
    select(
      .data$analysisId,
      .data$analysisDescription,
      .data$outcomeId,
      .data$outcomeName,
      .data$exposureId,
      .data$exposureName,
      .data$mdrr,
      .data$mdrrDiagnostic
    ) %>%
    inner_join(eases, by = c("analysisId", "exposureId")) %>%
    inner_join(medianDaysToEnd, by = c("outcomeId", "analysisId")) %>%
    inner_join(preExposure, by = c("exposureId", "outcomeId")) %>%
    mutate(unblind = .data$mdrrDiagnostic != "FAIL" &
             .data$easeDiagnostic != "FAIL" &
             # .data$medianDaysToEndDiagnostic != "FAIL" &
             .data$timeTrendDiagnostic != "FAIL" &
             .data$preExposureP != "FAIL")
  readr::write_csv(overview, file.path(outputFolder, "SccsDiagnosticsOverview.csv"))
}
