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

createCohortMethodAnalysisDetails <- function(workFolder) {
  
  # Load analyses specs from overall Excell sheet. 
  # Note: makes strong assumptions on cell ranges. Must update if format is updated.
  specs <- openxlsx::read.xlsx("inst/settings/HVRCAAnalysesSpecs.xlsx", colNames = FALSE)
  analysesSpecs <- t(specs[1:7, 4:ncol(specs)])
  colnames(analysesSpecs) <- analysesSpecs[1, ]
  analysesSpecs <- as.data.frame(analysesSpecs[-1, ])
  analysesSpecs <- analysesSpecs[analysesSpecs$method == "CohortMethod", ]
  
  # Use default covariate settings but set end day to -1, excluding the index date:
  covarSettings <- FeatureExtraction::createCovariateSettings(useDemographicsGender = TRUE,
                                                              useDemographicsAgeGroup = TRUE,
                                                              useDemographicsRace = TRUE,
                                                              useDemographicsEthnicity = TRUE,
                                                              useDemographicsIndexYear = TRUE,
                                                              useDemographicsIndexMonth = TRUE,
                                                              useConditionGroupEraLongTerm = TRUE,
                                                              useConditionGroupEraShortTerm = TRUE,
                                                              useDrugGroupEraLongTerm = TRUE,
                                                              useDrugGroupEraShortTerm = TRUE,
                                                              useProcedureOccurrenceLongTerm = TRUE,
                                                              useProcedureOccurrenceShortTerm = TRUE,
                                                              useDeviceExposureLongTerm = TRUE,
                                                              useDeviceExposureShortTerm = TRUE,
                                                              useMeasurementLongTerm = TRUE,
                                                              useMeasurementShortTerm = TRUE,
                                                              useMeasurementRangeGroupLongTerm = TRUE,
                                                              useObservationLongTerm = TRUE,
                                                              useObservationShortTerm = TRUE,
                                                              useCharlsonIndex = TRUE,
                                                              useDcsi = TRUE,
                                                              useChads2 = TRUE,
                                                              useChads2Vasc = TRUE,
                                                              endDays = -1)
  
  getDbCmDataArgs <- CohortMethod::createGetDbCohortMethodDataArgs(washoutPeriod = 365,
                                                                   restrictToCommonPeriod = TRUE,
                                                                   firstExposureOnly = TRUE,
                                                                   removeDuplicateSubjects = "remove all",
                                                                   studyStartDate = "20210101",
                                                                   studyEndDate = "",
                                                                   covariateSettings = covarSettings)
  
  
  
  createPsArgs <- CohortMethod::createCreatePsArgs(maxCohortSizeForFitting = 150000,
                                                   control = Cyclops::createControl(cvType = "auto",
                                                                                    startingVariance = 0.01,
                                                                                    noiseLevel = "quiet",
                                                                                    seed = 1,
                                                                                    resetCoefficients = TRUE,
                                                                                    tolerance = 2e-07,
                                                                                    cvRepetitions = 1))
  
  matchOnPsArgs <- CohortMethod::createMatchOnPsArgs(maxRatio = 100)
  
  fitOutcomeModelArgs <- CohortMethod::createFitOutcomeModelArgs(useCovariates = FALSE,
                                                                 modelType = "cox",
                                                                 stratified = TRUE)
  
  # row <- split(analysesSpecs, 1:nrow(analysesSpecs))[[5]]
  createAnalysis <- function(row) {
    createStudyPopArgs <-  CohortMethod::createCreateStudyPopulationArgs(removeSubjectsWithPriorOutcome = TRUE,
                                                                   priorOutcomeLookback = as.numeric(row$priorOutcomeLookback),
                                                                   minDaysAtRisk = 1,
                                                                   riskWindowStart = as.numeric(row$riskWindowStart),
                                                                   startAnchor = "cohort start",
                                                                   riskWindowEnd = as.numeric(row$riskWindowEnd),
                                                                   endAnchor = "cohort start")
    description <- sprintf("Cohort method, TAR = %s-%s, lookback = %s", 
                           row$riskWindowStart, 
                           row$riskWindowEnd, 
                           row$priorOutcomeLookback)
    description <- gsub("TAR = 0-9999", "ITT", description)
    cmAnalysis <- CohortMethod::createCmAnalysis(analysisId = as.numeric(row$analysisId),
                                                 description = description,
                                                 getDbCohortMethodDataArgs = getDbCmDataArgs,
                                                 createStudyPopArgs = createStudyPopArgs,
                                                 createPs = TRUE,
                                                 createPsArgs = createPsArgs,
                                                 matchOnPs = TRUE,
                                                 matchOnPsArgs = matchOnPsArgs,
                                                 fitOutcomeModel = TRUE,
                                                 fitOutcomeModelArgs = fitOutcomeModelArgs)
    return(cmAnalysis)
    
  }
  cmAnalysisList <- lapply(split(analysesSpecs, 1:nrow(analysesSpecs)), createAnalysis)
  
  CohortMethod::saveCmAnalysisList(cmAnalysisList, file.path(workFolder, "cmAnalysisList.json"))
}
