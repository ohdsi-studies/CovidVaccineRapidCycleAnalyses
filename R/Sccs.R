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

runSelfControlledCaseSeries <- function(connectionDetails,
                                        cdmDatabaseSchema,
                                        cohortDatabaseSchema,
                                        cohortTable,
                                        outputFolder,
                                        maxCores) {
  start <- Sys.time()
  sccsFolder <- file.path(outputFolder, "sccsOutput")
  if (!file.exists(sccsFolder)) {
    dir.create(sccsFolder)
  }

  sccsSummaryFile <- file.path(outputFolder, "sccsSummary.csv")
  if (!file.exists(sccsSummaryFile)) {
    eoList <- createTos()
    sccsAnalysisListFile <- system.file("settings", "sccsAnalysisList.json", package = "RcaCovidVaccine")
    sccsAnalysisList <- SelfControlledCaseSeries::loadSccsAnalysisList(sccsAnalysisListFile)
    analysesToExclude <- getAnalysesToExclude(method = "SCCS")
    sccsResult <- SelfControlledCaseSeries::runSccsAnalyses(
      connectionDetails = connectionDetails,
      cdmDatabaseSchema = cdmDatabaseSchema,
      exposureDatabaseSchema = cohortDatabaseSchema,
      exposureTable = cohortTable,
      outcomeDatabaseSchema = cohortDatabaseSchema,
      outcomeTable = cohortTable,
      sccsAnalysisList = sccsAnalysisList,
      exposureOutcomeList = eoList,
      outputFolder = sccsFolder,
      combineDataFetchAcrossOutcomes = FALSE,
      getDbSccsDataThreads = 1,
      createStudyPopulationThreads = min(3, maxCores),
      createSccsIntervalDataThreads = min(3, maxCores),
      fitSccsModelThreads = max(1, floor(maxCores / 4)),
      cvThreads = min(4, maxCores),
      analysesToExclude = analysesToExclude
    )
    sccsSummary <- SelfControlledCaseSeries::summarizeSccsAnalyses(sccsResult, sccsFolder)
    sccsSummary <- addCohortNames(sccsSummary, "exposureId", "exposureName")
    sccsSummary <- addCohortNames(sccsSummary, "outcomeId", "outcomeName")
    sccsSummary <- addAnalysisDescription(sccsSummary, sccsAnalysisList, "analysisId", "analysisDescription")
    readr::write_csv(sccsSummary, sccsSummaryFile)
  }

  fileName <- file.path(sccsFolder, "profiles.rds")
  if (!file.exists(fileName)) {
    message("Extracting likelihood profiles")
    sccsResult <- readRDS(file.path(sccsFolder, "outcomeModelReference.rds"))

    getProfile <- function(sccsModelFile) {
      model <- readRDS(file.path(sccsFolder, sccsModelFile))
      if (is.null(model$logLikelihoodProfile)) {
        return(NULL)
      } else {
        return(model$logLikelihoodProfile$`1000`)
      }
    }
    profiles <- lapply(sccsResult$sccsModelFile, getProfile)
    names(profiles) <- paste(sccsResult$analysisId, sccsResult$exposureId, sccsResult$outcomeId)
    saveRDS(profiles, fileName)
  }

  delta <- Sys.time() - start
  message(paste("Completed SCCS analyses in", signif(delta, 3), attr(delta, "units")))
}

createTos <- function() {
  pathToCsv <- system.file("settings", "TosOfInterest.csv", package = "RcaCovidVaccine")
  tosOfInterest <- readr::read_csv(pathToCsv, col_types = readr::cols())

  pathToCsv <- system.file("settings", "NegativeControls.csv", package = "RcaCovidVaccine")
  negativeControls <- readr::read_csv(pathToCsv, col_types = readr::cols())

  tos <- bind_rows(
    tosOfInterest,
    tosOfInterest %>%
      distinct(.data$exposureId, .data$exposureId2) %>%
      full_join(negativeControls %>% select(.data$outcomeId), by = character(0))
  )

  createTo <- function(i) {
    exposureOutcome <- SelfControlledCaseSeries::createExposureOutcome(
      exposureId = tos$exposureId[i],
      exposureId2 = tos$exposureId2[i],
      outcomeId = tos$outcomeId[i]
    )
    return(exposureOutcome)
  }
  tosList <- lapply(1:nrow(tos), createTo)
  return(tosList)
}
