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

runCohortMethod <- function(connectionDetails,
                            cdmDatabaseSchema,
                            cohortDatabaseSchema,
                            cohortTable,
                            tempEmulationSchema,
                            outputFolder,
                            maxCores,
                            fitOutcomeModels = TRUE,
                            outcomeModelBatch = NULL,
                            outcomeModelBatchCount = NULL) {
  cmOutputFolder <- file.path(outputFolder, "cmOutput")
  if (!file.exists(cmOutputFolder)) {
    dir.create(cmOutputFolder)
  }
  cmAnalysisListFile <- system.file("settings",
    "cmAnalysisList.json",
    package = "RcaCovidVaccine"
  )
  cmAnalysisList <- CohortMethod::loadCmAnalysisList(cmAnalysisListFile)
  tcosList <- createTcos(outputFolder = outputFolder)
  outcomesOfInterest <- getOutcomesOfInterest()

  # Changing analyses and outcomes for distributed compute
  if (!fitOutcomeModels) {
    message("Not fitting outcome models (as requested)")
    for (i in 1:length(cmAnalysisList)) {
      cmAnalysisList[[i]]$fitOutcomeModel <- FALSE
      cmAnalysisList[[i]]$fitOutcomeModelArgs <- NULL
    }
  } else if (!is.null(outcomeModelBatch)) {
    if (is.null(outcomeModelBatchCount)) {
      stop("If 'outcomeModelBatch' is specified, 'outcomeModelBatchCount' must also be specified")
    }
    if (outcomeModelBatch < 1 || outcomeModelBatch > outcomeModelBatchCount) {
      stop("'outcomeModelBatch' should be >= 1 and <= 'outcomeModelBatchCount'")
    }
    message(sprintf("Fitting outcome model batch %d of %d", outcomeModelBatch, outcomeModelBatchCount))
    outcomeIds <- c()
    for (i in 1:length(tcosList)) {
      outcomeIds <- c(outcomeIds, tcosList[[i]]$outcomeIds)
    }
    outcomeIds <- unique(outcomeIds)
    outcomeIds <- outcomeIds[order(outcomeIds)]
    batchIds <- rep(1:outcomeModelBatchCount, ceiling(length(outcomeIds) / outcomeModelBatchCount))
    batchIds <- batchIds[1:length(outcomeIds)]

    batchOutcomeIds <- outcomeIds[batchIds == outcomeModelBatch]

    # Subset outcomes to batch outcomes
    for (i in 1:length(tcosList)) {
      tcosList[[i]]$outcomeIds <- tcosList[[i]]$outcomeIds[tcosList[[i]]$outcomeIds %in% batchOutcomeIds]
    }
    outcomesOfInterest <- outcomesOfInterest[outcomesOfInterest %in% batchOutcomeIds]
  }

  analysesToExclude <- getAnalysesToExclude(method = "CohortMethod")

  results <- CohortMethod::runCmAnalyses(
    connectionDetails = connectionDetails,
    cdmDatabaseSchema = cdmDatabaseSchema,
    exposureDatabaseSchema = cohortDatabaseSchema,
    exposureTable = cohortTable,
    outcomeDatabaseSchema = cohortDatabaseSchema,
    outcomeTable = cohortTable,
    outputFolder = cmOutputFolder,
    tempEmulationSchema = tempEmulationSchema,
    cmAnalysisList = cmAnalysisList,
    targetComparatorOutcomesList = tcosList,
    refitPsForEveryOutcome = FALSE,
    refitPsForEveryStudyPopulation = FALSE,
    getDbCohortMethodDataThreads = 1, # min(3, maxCores),
    createStudyPopThreads = min(3, maxCores),
    createPsThreads = max(1, round(maxCores / 10)),
    psCvThreads = min(10, maxCores),
    trimMatchStratifyThreads = min(5, maxCores),
    fitOutcomeModelThreads = max(1, round(maxCores / 4)),
    outcomeCvThreads = min(4, maxCores),
    outcomeIdsOfInterest = outcomesOfInterest,
    analysesToExclude = analysesToExclude
  )

  if (!fitOutcomeModels || !is.null(outcomeModelBatch)) {
    message("Not fitting full set of outcome models, so quitting now")
    return(NULL)
  }
  message("Summarizing results")
  analysisSummary <- CohortMethod::summarizeAnalyses(
    referenceTable = results,
    outputFolder = cmOutputFolder
  )
  analysisSummary <- addCohortNames(analysisSummary, "targetId", "targetName")
  analysisSummary <- addCohortNames(analysisSummary, "comparatorId", "comparatorName")
  analysisSummary <- addCohortNames(analysisSummary, "outcomeId", "outcomeName")
  analysisSummary <- addAnalysisDescription(analysisSummary, cmAnalysisList, "analysisId", "analysisDescription")
  write.csv(analysisSummary, file.path(outputFolder, "analysisSummary.csv"), row.names = FALSE)

  message("Extracting likelihood profiles")
  getProfile <- function(outcomeModelFile) {
    model <- readRDS(file.path(cmOutputFolder, outcomeModelFile))
    if (is.null(model$logLikelihoodProfile)) {
      return(NULL)
    } else {
      logLikelihoodProfile <- model$logLikelihoodProfile
      logLikelihoodProfile <- tibble(
        point = as.numeric(names(logLikelihoodProfile)),
        value = logLikelihoodProfile
      )
      return(logLikelihoodProfile)
    }
  }
  profiles <- lapply(results$outcomeModelFile, getProfile)
  names(profiles) <- paste(results$analysisId, results$targetId, results$comparatorId, results$outcomeId)
  saveRDS(profiles, file.path(cmOutputFolder, "profiles.rds"))

  message("Computing covariate balance")
  balanceFolder <- file.path(outputFolder, "balance")
  if (!file.exists(balanceFolder)) {
    dir.create(balanceFolder)
  }
  subset <- results[results$outcomeId %in% outcomesOfInterest & results$analysisId == 17, ]
  subset <- subset[subset$strataFile != "", ]
  # Compute balance only once per Target-Comparator-Analysis:
  subset <- subset[!duplicated(subset[, c("targetId", "comparatorId", "analysisId")]), ]

  if (nrow(subset) > 0) {
    subset <- split(subset, seq(nrow(subset)))
    cluster <- ParallelLogger::makeCluster(min(5, maxCores))
    ParallelLogger::clusterApply(cluster, subset, computeCovariateBalance, cmOutputFolder = cmOutputFolder, balanceFolder = balanceFolder)
    ParallelLogger::stopCluster(cluster)
  }
}

computeCovariateBalance <- function(row, cmOutputFolder, balanceFolder) {
  outputFileName <- file.path(
    balanceFolder,
    sprintf("bal_t%s_c%s_a%s.rds", row$targetId, row$comparatorId, row$analysisId)
  )
  if (!file.exists(outputFileName)) {
    ParallelLogger::logTrace("Creating covariate balance file ", outputFileName)

    # Not computing balance per outcome, so need to create a new matched / stratified population
    # without excluding anyone with prior outcomes:
    cohortMethodDataFile <- file.path(cmOutputFolder, row$cohortMethodDataFile)
    cohortMethodData <- CohortMethod::loadCohortMethodData(cohortMethodDataFile)
    psFile <- file.path(cmOutputFolder, row$sharedPsFile)
    ps <- readRDS(psFile)
    cmAnalysisListFile <- system.file("settings",
      "cmAnalysisList.json",
      package = "RcaCovidVaccine"
    )
    cmAnalysisList <- CohortMethod::loadCmAnalysisList(cmAnalysisListFile)
    analysisArgs <- ParallelLogger::matchInList(cmAnalysisList, list(analysisId = row$analysisId))[[1]]
    createStudyPopArgs <- analysisArgs$createStudyPopArgs
    createStudyPopArgs$population <- ps
    studyPop <- do.call(CohortMethod::createStudyPopulation, createStudyPopArgs)

    # Note: implemented matching only. Add other variants as needed:
    if (analysisArgs$matchOnPs) {
      matchOnPsArgs <- analysisArgs$matchOnPsArgs
      matchOnPsArgs$population <- studyPop
      strataPop <- do.call(CohortMethod::matchOnPs, matchOnPsArgs)
    }

    balance <- CohortMethod::computeCovariateBalance(
      population = strataPop,
      cohortMethodData = cohortMethodData,
      maxCohortSize = 250000
    )
    saveRDS(balance, outputFileName)
  }
}

createTcos <- function(outputFolder) {
  pathToCsv <- system.file("settings", "TcosOfInterest.csv", package = "RcaCovidVaccine")
  tcosOfInterest <- readr::read_csv(pathToCsv, col_types = readr::cols())
  pathToCsv <- system.file("settings", "NegativeControls.csv", package = "RcaCovidVaccine")
  negativeControls <- readr::read_csv(pathToCsv, col_types = readr::cols())

  createTco <- function(i) {
    targetId <- tcosOfInterest$targetId[i]
    comparatorId <- tcosOfInterest$comparatorId[i]
    outcomeIds <- as.character(tcosOfInterest$outcomeIds[i])
    outcomeIds <- as.numeric(strsplit(outcomeIds, split = ";")[[1]])
    outcomeIds <- c(outcomeIds, negativeControls$outcomeId)
    excludeConceptIds <- as.character(tcosOfInterest$excludedCovariateConceptIds[i])
    if (length(excludeConceptIds) == 1 && is.na(excludeConceptIds)) {
      excludeConceptIds <- c()
    } else if (length(excludeConceptIds) > 0) {
      excludeConceptIds <- as.numeric(strsplit(excludeConceptIds, split = ";")[[1]])
    }
    tco <- CohortMethod::createTargetComparatorOutcomes(
      targetId = targetId,
      comparatorId = comparatorId,
      outcomeIds = outcomeIds,
      excludedCovariateConceptIds = excludeConceptIds
    )
    return(tco)
  }
  tcosList <- lapply(1:nrow(tcosOfInterest), createTco)
  return(tcosList)
}
