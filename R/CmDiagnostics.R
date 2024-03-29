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
generateCmDiagnostics <- function(outputFolder) {
  cmDiagnosticsFolder <- file.path(outputFolder, "cmDiagnostics")
  if (!file.exists(cmDiagnosticsFolder)) {
    dir.create(cmDiagnosticsFolder)
  }
  cmAnalysisListFile <- system.file("settings",
                                    "cmAnalysisList.json",
                                    package = "RcaCovidVaccine"
  )
  cmAnalysisList <- CohortMethod::loadCmAnalysisList(cmAnalysisListFile)
  
  # Negative control plots -------------
  cmSummaryFile <- file.path(outputFolder, "analysisSummary.csv")
  cmSummary <- readr::read_csv(cmSummaryFile, col_types = readr::cols())
  
  pathToCsv <- system.file("settings", "NegativeControls.csv", package = "RcaCovidVaccine")
  negativeControls <- readr::read_csv(pathToCsv, col_types = readr::cols())
  
  cmSummary <- cmSummary %>%
    filter(.data$outcomeId %in% !!negativeControls$outcomeId)
  
  # subset <- split(cmSummary, paste(cmSummary$analysisId, cmSummary$targetId, cmSummary$comparatorId))[[1]]
  plotCalibration <- function(subset) {
    fileName <- file.path(cmDiagnosticsFolder, sprintf("NcDistribution_a%d_t%d_c%d.png", subset$analysisId[1], subset$targetId[1], subset$comparatorId[1]))
    if (!file.exists(fileName)) {
      # title <- sprintf("%s\n%s\n%s",
      #                  gsub(".*exposure to ", "", subset$targetName[1]),
      #                  gsub(".*exposure to ", "", subset$comparatorName[1]),
      #                  subset$analysisDescription[1])
      title <- subset$analysisDescription[1]
      EmpiricalCalibration::plotCalibrationEffect(
        logRrNegatives = subset$logRr,
        seLogRrNegatives = subset$seLogRr,
        xLabel = "Hazard Ratio",
        title = title,
        showCis = TRUE,
        fileName = fileName
      )
    }
    null <- EmpiricalCalibration::fitMcmcNull(subset$logRr, subset$seLogRr)
    ease <- EmpiricalCalibration::computeExpectedAbsoluteSystematicError(null)
    return(data.frame(
      targetId = subset$targetId[1],
      comparatorId = subset$comparatorId[1],
      analysisId = subset$analysisId[1],
      ease = ease$ease
    ))
  }
  eases <- plyr::llply(split(cmSummary, paste(cmSummary$analysisId, cmSummary$targetId, cmSummary$comparatorId)), plotCalibration)
  eases <- bind_rows(eases)
  eases <- eases %>%
    mutate(easeDiagnostic = case_when(
      abs(.data$ease) < 0.1 ~ "PASS",
      abs(.data$ease) < 0.25 ~ "WARNING",
      TRUE ~ "FAIL"
    ))
  # Warning: very specific code for current analysis design matrix. Needs updating when analyses are changed:
  eases <- bind_rows(
    eases,
    eases %>%
      filter(.data$analysisId <= 19) %>%
      mutate(analysisId = .data$analysisId + 4),
    eases %>%
      filter(.data$analysisId <= 19) %>%
      mutate(analysisId = .data$analysisId + 8),
    eases %>%
      filter(.data$analysisId <= 19) %>%
      mutate(analysisId = .data$analysisId + 12)
  )
  
  
  # PS distribution plots -----------------------------------------
  subset <- readRDS(file.path(outputFolder, "cmOutput", "outcomeModelReference.rds")) %>%
    distinct(.data$targetId, .data$comparatorId, .data$sharedPsFile) %>%
    addCohortNames(IdColumnName = "targetId", nameColumnName = "targetName") %>%
    addCohortNames(IdColumnName = "comparatorId", nameColumnName = "comparatorName")
  
  # row <- split(subset, 1:nrow(subset))[[2]]
  createPsPlot <- function(row) {
    ps <- readRDS(file.path(outputFolder, "cmOutput", row$sharedPsFile))
    
    fileName <- file.path(cmDiagnosticsFolder, sprintf("Ps_t%d_c%d.png", row$targetId[1], row$comparatorId[1]))
    if (!file.exists(fileName)) {
      message("Generating plot ", fileName)
      
      
      CohortMethod::plotPs(
        data = ps,
        targetLabel = gsub(".*exposure to ", "", row$targetName),
        comparatorLabel = gsub(".*exposure to ", "", row$comparatorName),
        showAucLabel = TRUE,
        showEquiposeLabel = TRUE,
        fileName = fileName
      )
    }
    equipoiseBounds <- c(0.3, 0.7)
    equipoise <- mean(ps$preferenceScore >= equipoiseBounds[1] & ps$preferenceScore <= equipoiseBounds[2])
    return(data.frame(
      targetId = row$targetId[1],
      comparatorId = row$comparatorId[1],
      equipoise = equipoise
    ))
  }
  equipoises <- plyr::llply(split(subset, 1:nrow(subset)), createPsPlot)
  equipoises <- bind_rows(equipoises)
  equipoises <- equipoises %>%
    mutate(equipoiseDiagnostic = case_when(
      .data$equipoise >= 0.5 ~ "PASS",
      .data$equipoise >= 0.1 ~ "WARNING",
      TRUE ~ "FAIL"
    ))
  
  # Propensity model ----------------------------------------------------
  subset <- readRDS(file.path(outputFolder, "cmOutput", "outcomeModelReference.rds")) %>%
    distinct(.data$targetId, .data$comparatorId, .data$sharedPsFile, .data$cohortMethodDataFile) %>%
    addCohortNames(IdColumnName = "targetId", nameColumnName = "targetName") %>%
    addCohortNames(IdColumnName = "comparatorId", nameColumnName = "comparatorName")
  
  createPropensityModelTable <- function(row) {
    fileName <- file.path(cmDiagnosticsFolder, sprintf("PropensityModel_t%d_c%d.csv", row$targetId[1], row$comparatorId[1]))
    if (!file.exists(fileName)) {
      ps <- readRDS(file.path(outputFolder, "cmOutput", row$sharedPsFile))
      cmData <- CohortMethod::loadCohortMethodData(file.path(outputFolder, "cmOutput", row$cohortMethodDataFile))
      model <- CohortMethod::getPsModel(ps, cmData)
      readr::write_excel_csv(model, fileName)
    }
  }
  plyr::l_ply(split(subset, 1:nrow(subset)), createPropensityModelTable)
  
  # Covariate balance plots and tables 1 ------------------------------------------------------------------
  balanceFiles <- list.files(file.path(outputFolder, "balance"), "bal.*.rds", full.names = TRUE)
  pathToCsv <- system.file("settings", "CohortsToCreate.csv", package = "RcaCovidVaccine")
  cohortsToCreate <- read.csv(pathToCsv)
  
  # balanceFile = balanceFiles[1]
  plotBalance <- function(balanceFile) {
    balance <- readRDS(balanceFile)
    targetId <- gsub(".*_t", "", gsub("_c[0-9].*", "", balanceFile))
    comparatorId <- gsub(".*_c", "", gsub("_a.*", "", balanceFile))
    analysisId <- gsub(".*_a", "", gsub(".rds", "", balanceFile))
    beforeLabel <- "Before matching"
    afterLabel <- "After matching"
    
    fileName <- file.path(cmDiagnosticsFolder, sprintf("Balance_t%s_c%s_a%s.png", targetId, comparatorId, analysisId))
    if (!file.exists(fileName)) {
      CohortMethod::plotCovariateBalanceScatterPlot(
        balance = balance,
        threshold = 0.1,
        showCovariateCountLabel = TRUE,
        showMaxLabel = TRUE,
        beforeLabel = beforeLabel,
        afterLabel = afterLabel,
        fileName = fileName
      )
    }
    fileName <- file.path(cmDiagnosticsFolder, sprintf("BalanceTopVars_t%s_c%s_a%s.png", targetId, comparatorId, analysisId))
    if (!file.exists(fileName)) {
      CohortMethod::plotCovariateBalanceOfTopVariables(
        balance = balance,
        beforeLabel = beforeLabel,
        afterLabel = afterLabel,
        fileName = fileName
      )
    }
    fileName <- file.path(cmDiagnosticsFolder, sprintf("Balance_t%s_c%s_a%s.png", targetId, comparatorId, analysisId))
    if (!file.exists(fileName)) {
      table1 <- CohortMethod::createCmTable1(
        balance = balance,
        beforeLabel = beforeLabel,
        afterLabel = afterLabel
      )
      readr::write_excel_csv(table1, fileName)
    }
    fileName <- file.path(cmDiagnosticsFolder, sprintf("BalanceAllVars_t%s_c%s_a%s.csv", targetId, comparatorId, analysisId))
    if (!file.exists(fileName)) {
      readr::write_excel_csv(balance, fileName)
    }
    fileName <- file.path(cmDiagnosticsFolder, sprintf("BalancePrevVsPrev_t%s_c%s_a%s.png", targetId, comparatorId, analysisId))
    if (!file.exists(fileName)) {
      targetLabel <- cohortsToCreate$name[cohortsToCreate$cohortId == targetId]
      comparatorLabel <- cohortsToCreate$name[cohortsToCreate$cohortId == comparatorId]
      CohortMethod::plotCovariatePrevalence(
        balance = balance,
        beforeLabel = beforeLabel,
        afterLabel = afterLabel,
        targetLabel = targetLabel,
        comparatorLabel = comparatorLabel,
        fileName = fileName
      )
    }
    
    return(data.frame(
      targetId = as.numeric(targetId),
      comparatorId = as.numeric(comparatorId),
      analysisId = as.numeric(analysisId),
      maxAbsSdm = max(abs(balance$afterMatchingStdDiff), na.rm = TRUE)
    ))
  }
  balance <- plyr::llply(balanceFiles, plotBalance)
  balance <- bind_rows(balance)
  balance <- balance %>%
    mutate(balanceDiagnostic = case_when(
      .data$maxAbsSdm < 0.1 ~ "PASS",
      TRUE ~ "FAIL"
    ))
  
  # MDRR and follow up ---------------------------
  outcomeIds <- getOutcomesOfInterest()
  subset <- readRDS(file.path(outputFolder, "cmOutput", "outcomeModelReference.rds")) %>%
    filter(.data$outcomeId %in% outcomeIds) %>%
    distinct(.data$analysisId, .data$targetId, .data$comparatorId, .data$outcomeId, .data$strataFile)
  
  plotFollowUp <- function(row) {
    strataPop <- readRDS(file.path(outputFolder, "cmOutput", row$strataFile))
    fileName <- file.path(cmDiagnosticsFolder, sprintf("FollowUp_a%s_o%s.png", row$analysisId, row$outcomeId))
    if (!file.exists(fileName)) {
      row <- row %>%
        addCohortNames(IdColumnName = "targetId", nameColumnName = "targetName") %>%
        addCohortNames(IdColumnName = "comparatorId", nameColumnName = "comparatorName") %>%
        addCohortNames(IdColumnName = "outcomeId", nameColumnName = "outcomeName") %>%
        addAnalysisDescription(cmAnalysisList)
      title <- sprintf("%s\n%s", row$analysisDescription, gsub(".*\\] ", "", row$outcomeName))
      CohortMethod::plotFollowUpDistribution(
        population = strataPop,
        targetLabel = gsub(".*exposure to ", "", row$targetName),
        comparatorLabel = gsub(".*exposure to ", "", row$comparatorName),
        title = title,
        fileName = fileName
      )
    }
    # medianFollowUp <- strataPop %>%
    #   group_by(.data$treatment) %>%
    #   summarise(medianFollowUp = median(.data$survivalTime))
    
    # Wilcoxon Rank Sum Test always significant due to large sample, so not informative:
    # wilcox.test(survivalTime ~ treatment, data = strataPop)
    # return(data.frame(targetId = row$targetId,
    #                   comparatorId = row$comparatorId,
    #                   outcomeId = row$outcomeId,
    #                   analysisId = row$analysisId,
    #                   medianFollowUpT = medianFollowUp$medianFollowUp[medianFollowUp$treatment == 1],
    #                   medianFollowUpC = medianFollowUp$medianFollowUp[medianFollowUp$treatment == 0]))
  }
  plyr::l_ply(split(subset, 1:nrow(subset)), plotFollowUp)
  # followUp <- plyr::llply(split(subset, 1:nrow(subset)), plotFollowUp)
  # followUp <- bind_rows(followUp)
  # followUp <- followUp %>%
  #   mutate(followUpPass = .data$medianFollowUpT / .data$medianFollowUpC > 0.5 & .data$medianFollowUpT / .data$medianFollowUpC < 2)
  
  plotAttrition <- function(row) {
    fileName <- file.path(cmDiagnosticsFolder, sprintf("Attrition_a%s_c%s_o%s.png", row$analysisId, row$comparatorId, row$outcomeId))
    if (!file.exists(fileName)) {
      strataPop <- readRDS(file.path(outputFolder, "cmOutput", row$strataFile))
      row <- row %>%
        addCohortNames(IdColumnName = "targetId", nameColumnName = "targetName") %>%
        addCohortNames(IdColumnName = "comparatorId", nameColumnName = "comparatorName") %>%
        addCohortNames(IdColumnName = "outcomeId", nameColumnName = "outcomeName") %>%
        addAnalysisDescription(cmAnalysisList)
      CohortMethod::drawAttritionDiagram(
        object = strataPop,
        targetLabel = gsub(" .*$", "", gsub(".*exposure to ", "", row$targetName)),
        comparatorLabel = gsub(" .*$", "", gsub(".*exposure to ", "", row$comparatorName)),
        fileName = fileName
      )
    }
  }
  plyr::l_ply(split(subset, 1:nrow(subset)), plotAttrition)
  
  # row <- split(subset, 1:nrow(subset))[[1]]
  computeMdrr <- function(row) {
    strataPop <- readRDS(file.path(outputFolder, "cmOutput", row$strataFile))
    attrition <- CohortMethod::getAttritionTable(strataPop)
    # Using matching (implying ATT estimator), so only interested in attrition of target cohort when worried
    # about generalizability:
    attritionFraction <- 1 - (attrition$targetExposures[nrow(attrition)] / attrition$targetExposures[1])
    CohortMethod::computeMdrr(population = strataPop) %>%
      mutate(
        analysisId = row$analysisId,
        targetId = row$targetId,
        comparatorId = row$comparatorId,
        outcomeId = row$outcomeId,
        attritionFraction = attritionFraction
      ) %>%
      return()
  }
  results <- plyr::llply(split(subset, 1:nrow(subset)), computeMdrr)
  results <- bind_rows(results) %>%
    addCohortNames(IdColumnName = "targetId", nameColumnName = "targetName") %>%
    addCohortNames(IdColumnName = "comparatorId", nameColumnName = "comparatorName") %>%
    addCohortNames(IdColumnName = "outcomeId", nameColumnName = "outcomeName") %>%
    addAnalysisDescription(cmAnalysisList)
  readr::write_csv(results, file.path(cmDiagnosticsFolder, "MDRR.csv"))
  
  # Time to event ------------------------------------
  outcomeIds <- getOutcomesOfInterest()
  subset <- readRDS(file.path(outputFolder, "cmOutput", "outcomeModelReference.rds")) %>%
    filter(.data$outcomeId %in% outcomeIds) %>%
    distinct(.data$targetId, .data$comparatorId, .data$outcomeId, .data$cohortMethodDataFile) %>%
    addCohortNames(IdColumnName = "outcomeId", nameColumnName = "outcomeName")
  
  plotTimeToEventForOutcome <- function(row, cmData, population) {
    fileName <- file.path(cmDiagnosticsFolder, sprintf("TimeToEvent_t%s_c%s_o%s.png", row$targetId, row$comparatorId, row$outcomeId))
    if (!file.exists(fileName)) {
      title <- gsub(".*\\] ", "", row$outcomeName)
      CohortMethod::plotTimeToEvent(
        cohortMethodData = cmData,
        population = population,
        outcomeId = row$outcomeId,
        firstExposureOnly = TRUE,
        restrictToCommonPeriod = TRUE,
        washoutPeriod = 365,
        removeDuplicateSubjects = TRUE,
        minDaysAtRisk = 0,
        riskWindowStart = 0,
        startAnchor = "cohort start",
        riskWindowEnd = 28,
        endAnchor = "cohort start",
        showFittedLines = TRUE,
        includePostIndexTime = FALSE,
        title = title,
        fileName = fileName
      )
    }
  }
  
  plotTimeToEventForTcos <- function(tcosRows) {
    cmData <- CohortMethod::loadCohortMethodData(file.path(outputFolder, "cmOutput", tcosRows$cohortMethodDataFile[1]))
    population <- cmData$cohorts %>%
      collect()
    plyr::l_ply(split(tcosRows, tcosRows$outcomeId), plotTimeToEventForOutcome, cmData = cmData, population = population)
  }
  
  plyr::l_ply(split(subset, subset$cohortMethodDataFile), plotTimeToEventForTcos)
  
  # Write overview ------------------------
  overview <- results %>%
    mutate(mdrrDiagnostic = case_when(
      .data$mdrr < 2 ~ "PASS",
      .data$mdrr < 10 ~ "WARNING",
      TRUE ~ "FAIL"
    )) %>%
    mutate(attritionDiagnostic = case_when(
      .data$attritionFraction < 0.1 ~ "PASS",
      .data$attritionFraction < 0.5 ~ "WARNING",
      TRUE ~ "FAIL"
    )) %>%
    select(
      .data$analysisId,
      .data$analysisDescription,
      .data$outcomeId,
      .data$outcomeName,
      .data$targetId,
      .data$targetName,
      .data$comparatorId,
      .data$comparatorName,
      .data$mdrr,
      .data$mdrrDiagnostic,
      .data$attritionFraction,
      .data$attritionDiagnostic
    ) %>%
    inner_join(eases, by = c("analysisId", "targetId", "comparatorId")) %>%
    inner_join(equipoises, by = c("targetId", "comparatorId")) %>%
    inner_join(select(balance, -.data$analysisId), by = c("targetId", "comparatorId")) %>%
    # inner_join(followUp,  by = c("analysisId", "targetId", "comparatorId", "outcomeId")) %>%
    mutate(unblind = .data$mdrrDiagnostic != "FAIL" &
             .data$attritionDiagnostic != "FAIL" &
             .data$easeDiagnostic != "FAIL" &
             .data$equipoiseDiagnostic != "FAIL" &
             .data$balanceDiagnostic != "FAIL")
  readr::write_csv(overview, file.path(outputFolder, "CmDiagnosticsOverview.csv"))
}
