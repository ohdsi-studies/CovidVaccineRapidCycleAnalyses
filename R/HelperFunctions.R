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

addCohortNames <- function(data, IdColumnName = "cohortDefinitionId", nameColumnName = "cohortName") {
  pathToCsv <- system.file("settings", "CohortsToCreate.csv", package = "RcaCovidVaccine")
  cohortsToCreate <- read.csv(pathToCsv)
  pathToCsv <- system.file("settings", "NegativeControls.csv", package = "RcaCovidVaccine")
  negativeControls <- read.csv(pathToCsv)

  idToName <- data.frame(
    cohortId = c(
      cohortsToCreate$cohortId,
      negativeControls$outcomeId
    ),
    cohortName = c(
      as.character(cohortsToCreate$atlasName),
      as.character(negativeControls$outcomeName)
    )
  )
  idToName <- idToName[order(idToName$cohortId), ]
  idToName <- idToName[!duplicated(idToName$cohortId), ]
  names(idToName)[1] <- IdColumnName
  names(idToName)[2] <- nameColumnName
  data <- merge(data, idToName, all.x = TRUE)
  # Change order of columns:
  idCol <- which(colnames(data) == IdColumnName)
  if (idCol < ncol(data) - 1) {
    data <- data[, c(1:idCol, ncol(data), (idCol + 1):(ncol(data) - 1))]
  }
  return(data)
}

getOutcomesOfInterest <- function() {
  pathToCsv <- system.file("settings", "TcosOfInterest.csv", package = "RcaCovidVaccine")
  tcosOfInterest <- readr::read_csv(pathToCsv, col_types = readr::cols())
  outcomeIds <- as.character(tcosOfInterest$outcomeIds)
  outcomeIds <- do.call("c", (strsplit(outcomeIds, split = ";")))
  outcomeIds <- unique(as.numeric(outcomeIds))
  return(outcomeIds)
}

addAnalysisDescription <- function(data, analysisSettingsList, IdColumnName = "analysisId", nameColumnName = "analysisDescription") {
  idToName <- lapply(analysisSettingsList, function(x) data.frame(analysisId = x$analysisId, description = as.character(x$description)))
  idToName <- do.call("rbind", idToName)
  names(idToName)[1] <- IdColumnName
  names(idToName)[2] <- nameColumnName
  data <- merge(data, idToName, all.x = TRUE)
  # Change order of columns:
  idCol <- which(colnames(data) == IdColumnName)
  if (idCol < ncol(data) - 1) {
    data <- data[, c(1:idCol, ncol(data), (idCol + 1):(ncol(data) - 1))]
  }
  return(data)
}

getAnalysesToExclude <- function(method = "CohortMethod") {
  pathToCsv <- system.file("settings", "NegativeControls.csv", package = "RcaCovidVaccine")
  negativeControls <- readr::read_csv(pathToCsv, col_types = readr::cols())

  specs <- openxlsx::read.xlsx(system.file("settings", "HVRCAAnalysesSpecs.xlsx", package = "RcaCovidVaccine"), colNames = FALSE)
  methodRow <- 1
  analysisIdRow <- 2
  cohortNameCol <- 1
  cohortIdCol <- 2
  startRow <- 7
  startCol <- 5
  methods <- specs[methodRow, startCol:ncol(specs)]
  methodCols <- which(methods == method) + startCol - 1

  analysesToExclude <- data.frame()
  for (rowIdx in startRow:nrow(specs)) {
    for (columnIdx in methodCols) {
      if (specs[rowIdx, columnIdx] == "N") {
        if (tolower(specs[rowIdx, cohortNameCol]) == "negative controls") {
          analysesToExclude <- rbind(
            analysesToExclude,
            data.frame(
              analysisId = as.numeric(specs[analysisIdRow, columnIdx]),
              outcomeId = negativeControls$outcomeId
            )
          )
        } else {
          analysesToExclude <- rbind(
            analysesToExclude,
            data.frame(
              analysisId = as.numeric(specs[analysisIdRow, columnIdx]),
              outcomeId = as.numeric(specs[rowIdx, cohortIdCol])
            )
          )
        }
      }
    }
  }
  return(analysesToExclude)
}

getDiscoveryAnalyses <- function() {
  specs <- openxlsx::read.xlsx(system.file("settings", "HVRCAAnalysesSpecs.xlsx", package = "RcaCovidVaccine"), colNames = FALSE)
  methodRow <- 1
  analysisIdRow <- 2
  cohortNameCol <- 1
  cohortIdCol <- 2
  startRow <- 7
  startCol <- 5
  methods <- specs[methodRow, startCol:ncol(specs)]
  methodCols <- which(methods %in% c("SCCS", "CohortMethod")) + startCol - 1
  
  discoveryAnalyses <- data.frame()
  for (rowIdx in startRow:nrow(specs)) {
    for (columnIdx in methodCols) {
      if (specs[rowIdx, columnIdx] == "D") {
        discoveryAnalyses <- rbind(
          discoveryAnalyses,
          data.frame(
            method = specs[methodRow, columnIdx],
            analysisId = as.numeric(specs[analysisIdRow, columnIdx]),
            outcomeId = as.numeric(specs[rowIdx, cohortIdCol])
          )
        )
      }
    }
  }
  return(discoveryAnalyses)
}
