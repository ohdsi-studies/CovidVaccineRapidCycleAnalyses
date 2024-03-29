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

createSccsAnalysisDetails <- function(workFolder) {
  
  # Load analyses specs from overall Excell sheet. 
  # Note: makes strong assumptions on cell ranges. Must update if format is updated.
  specs <- openxlsx::read.xlsx("inst/settings/HVRCAAnalysesSpecs.xlsx", colNames = FALSE)
  analysesSpecs <- t(specs[1:7, 4:ncol(specs)])
  colnames(analysesSpecs) <- analysesSpecs[1, ]
  analysesSpecs <- as.data.frame(analysesSpecs[-1, ])
  analysesSpecs <- analysesSpecs[analysesSpecs$method == "SCCS", ]
  
  
  createStudyPopulationArgs <- SelfControlledCaseSeries::createCreateStudyPopulationArgs(naivePeriod = 365,
                                                                                         firstOutcomeOnly = TRUE)
  
  
  
  covarPreExp <- SelfControlledCaseSeries::createEraCovariateSettings(label = "Pre-exposure",
                                                                      includeEraIds = "exposureId",
                                                                      start = -30,
                                                                      end = -1,
                                                                      endAnchor = "era start")
  
  calendarTimeSettings <- SelfControlledCaseSeries::createCalendarTimeCovariateSettings(calendarTimeKnots = 5, 
                                                                                        allowRegularization = TRUE,
                                                                                        computeConfidenceIntervals = FALSE)
  
  controlIntervalSettings <- SelfControlledCaseSeries::createControlIntervalSettings(
    includeEraIds = "exposureId",
    start = 0,
    startAnchor = "era start",
    end = 999999,
    endAnchor = "era start"
  )
  
  fitSccsModelArgs <- SelfControlledCaseSeries::createFitSccsModelArgs(
    control = Cyclops::createControl(cvType = "auto", 
                                     selectorType = "byPid", 
                                     startingVariance = 0.1,
                                     seed = 1, 
                                     resetCoefficients = TRUE,
                                     noiseLevel = "quiet")
  )
  
  # row <- split(analysesSpecs, 1:nrow(analysesSpecs))[[1]]
  createAnalysis <- function(row) {
    getDbSccsDataArgs <- SelfControlledCaseSeries::createGetDbSccsDataArgs(studyStartDate = row$studyStartDate,
                                                                           studyEndDate = "",
                                                                           maxCasesPerOutcome = 1e6,
                                                                           exposureIds = c("exposureId", "exposureId2"),
                                                                           deleteCovariatesSmallCount = 0)
    
    covarExposureOfInt <- SelfControlledCaseSeries::createEraCovariateSettings(label = "Main",
                                                                               includeEraIds = "exposureId",
                                                                               start = as.numeric(row$riskWindowStart),
                                                                               startAnchor = "era start",
                                                                               end = as.numeric(row$riskWindowEnd),
                                                                               endAnchor = "era start",
                                                                               profileLikelihood = TRUE)
    
    covarSecondDose <- SelfControlledCaseSeries::createEraCovariateSettings(label = "Second",
                                                                            includeEraIds = "exposureId2",
                                                                            start = as.numeric(row$riskWindowStart),
                                                                            startAnchor = "era start",
                                                                            end = as.numeric(row$riskWindowEnd),
                                                                            endAnchor = "era start",
                                                                            profileLikelihood = TRUE)
    
    if (row$postExposureOnly == "no") {
      createSccsIntervalDataArgs <- SelfControlledCaseSeries::createCreateSccsIntervalDataArgs(
        eraCovariateSettings = list(covarExposureOfInt,
                                    covarSecondDose,
                                    covarPreExp),
        calendarTimeCovariateSettings = calendarTimeSettings
      )
      
      description <- sprintf("SCCS, TAR = %s-%s", 
                             row$riskWindowStart,
                             row$riskWindowEnd)
      if (row$studyStartDate == "20180101") {
        description <- paste0(description, ", longer lookback")
      }
      
      sccsAnalysis <- SelfControlledCaseSeries::createSccsAnalysis(analysisId = as.numeric(row$analysisId),
                                                                   description = description,
                                                                   getDbSccsDataArgs = getDbSccsDataArgs,
                                                                   createStudyPopulationArgs = createStudyPopulationArgs,
                                                                   createSccsIntervalDataArgs = createSccsIntervalDataArgs,
                                                                   fitSccsModelArgs = fitSccsModelArgs)
    } else {
      # Post exposure only
      createScriIntervalDataArgsAllPost <- SelfControlledCaseSeries::createCreateScriIntervalDataArgs(
        eraCovariateSettings = list(covarExposureOfInt,
                                    covarSecondDose),
        controlIntervalSettings = controlIntervalSettings
      )
      
      description <- sprintf("SCCS, TAR = %s-%s, post-exposure only", 
                             row$riskWindowStart,
                             row$riskWindowEnd)
      sccsAnalysis <- SelfControlledCaseSeries::createSccsAnalysis(analysisId = as.numeric(row$analysisId),
                                                                   description = description,
                                                                   getDbSccsDataArgs = getDbSccsDataArgs,
                                                                   createStudyPopulationArgs = createStudyPopulationArgs,
                                                                   design = "SCRI",
                                                                   createScriIntervalDataArgs = createScriIntervalDataArgsAllPost,
                                                                   fitSccsModelArgs = fitSccsModelArgs)
    }
    return(sccsAnalysis)
  }
  sccsAnalysisList <- lapply(split(analysesSpecs, 1:nrow(analysesSpecs)), createAnalysis)
  
  SelfControlledCaseSeries::saveSccsAnalysisList(sccsAnalysisList, file.path(workFolder, "sccsAnalysisList.json"))
}
