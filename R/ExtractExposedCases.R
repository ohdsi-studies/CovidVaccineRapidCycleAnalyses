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

# outputFolder <- "s:/temp/RcaCovidVaccine"
# source("R/HelperFunctions.R")
# library(dplyr)
extractExposedCasesSccs <- function(outputFolder) {
  vaccineCovariateId <- 1000
  outcomeIds <- getOutcomesOfInterest()
  omr <- readRDS(file.path(outputFolder, "sccsOutput", "outcomeModelReference.rds")) %>%
    filter(.data$outcomeId %in% outcomeIds)


  # i <- 1
  extractCases <- function(i) {
    studyPop <- readRDS(file.path(outputFolder, "sccsOutput", omr$studyPopFile[i]))
    sccsIntervalData <- SelfControlledCaseSeries::loadSccsIntervalData(file.path(
      outputFolder,
      "sccsOutput",
      omr$sccsIntervalDataFile[i]
    ))
    exposedCaseIds <- sccsIntervalData$covariates %>%
      filter(.data$covariateId == vaccineCovariateId) %>%
      inner_join(sccsIntervalData$outcomes %>% filter(.data$y > 0),
        by = c("rowId", "stratumId")
      ) %>%
      pull(.data$stratumId)
    exposedOutcomes <- studyPop$cases %>%
      inner_join(studyPop$outcomes, by = "caseId") %>%
      filter(.data$caseId %in% exposedCaseIds) %>%
      collect() %>%
      mutate(cohortStartDate = .data$startDate + .data$outcomeDay) %>%
      select(subjectId = .data$personId, .data$cohortStartDate)
    exposedOutcomes <- exposedOutcomes %>%
      mutate(
        exposureId = omr$exposureId[i],
        outcomeId = omr$outcomeId[i],
        analysisId = omr$analysisId[i]
      )
    return(exposedOutcomes)
  }
  results <- lapply(1:nrow(omr), extractCases)
  results <- bind_rows(results)
  saveRDS(results, file.path(outputFolder, "ExposedCasesSccs.rds"))
}

extractExposedCasesCohortMethod <- function(outputFolder) {
  outcomeIds <- getOutcomesOfInterest()
  omr <- readRDS(file.path(outputFolder, "cmOutput", "outcomeModelReference.rds")) %>%
    filter(.data$outcomeId %in% outcomeIds)
  cohortMethodData <- CohortMethod::loadCohortMethodData(file.path(outputFolder, "cmOutput", omr$cohortMethodDataFile[1]))
  personSeqIdToPersonId <- cohortMethodData$cohorts %>%
    select(.data$personSeqId, .data$personId) %>%
    collect()

  # i <- 1
  extractCases <- function(i) {
    strataPop <- readRDS(file.path(outputFolder, "cmOutput", omr$strataFile[i]))

    exposedOutcomes <- strataPop %>%
      filter(.data$outcomeCount > 0) %>%
      mutate(
        cohortStartDate = .data$cohortStartDate + .data$daysToEvent,
        exposureId = if_else(.data$treatment == 1, omr$targetId[i], omr$comparatorId[i])
      ) %>%
      inner_join(personSeqIdToPersonId, by = "personSeqId") %>%
      select(subjectId = .data$personId, .data$cohortStartDate, .data$exposureId)

    exposedOutcomes <- exposedOutcomes %>%
      mutate(
        outcomeId = omr$outcomeId[i],
        analysisId = omr$analysisId[i]
      )
    return(exposedOutcomes)
  }
  results <- lapply(1:nrow(omr), extractCases)
  results <- bind_rows(results)
  saveRDS(results, file.path(outputFolder, "ExposedCasesCm.rds"))
}
