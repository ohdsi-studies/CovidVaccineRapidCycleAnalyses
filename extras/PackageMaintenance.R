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

# Format and check code ---------------------------------------------------
remotes::install_github("ohdsi/OhdsiRTools")
styler::style_pkg()
OhdsiRTools::checkUsagePackage("RcaCovidVaccine")
OhdsiRTools::updateCopyrightYearFolder()
devtools::spell_check()

# Create manual -----------------------------------------------------------
unlink("extras/RcaCovidVaccine.pdf")
shell("R CMD Rd2pdf ./ --output=extras/RcaCovidVaccine.pdf")

# Insert cohort definitions from ATLAS into package -----------------------
remotes::install_github("ohdsi/ROhdsiWebApi")
ROhdsiWebApi::authorizeWebApi(baseUrl = Sys.getenv("baseUrl"),
                              authMethod = "windows")
ROhdsiWebApi::insertCohortDefinitionSetInPackage(fileName = "inst/settings/CohortsToCreate.csv",
                                                 baseUrl = Sys.getenv("baseUrl"),
                                                 insertTableSql = FALSE,
                                                 insertCohortCreationR = FALSE,
                                                 generateStats = FALSE,
                                                 packageName = "RcaCovidVaccine")

# Create analysis details -------------------------------------------------
source("extras/CreateCohortMethodAnalysisDetails.R")
createCohortMethodAnalysisDetails("inst/settings/")
source("extras/CreateSccsAnalysisDetails.R")
createSccsAnalysisDetails("inst/settings/")

# Capture dependencies in renv lock file --------------------------------------
OhdsiRTools::createRenvLockFile(rootPackage = "RcaCovidVaccine",
                                includeRootPackage = FALSE,
                                mode = "description",
                                additionalRequiredPackages = c("keyring", 
                                                               "OhdsiRTools",
                                                               "checkmate", 
                                                               "DT", 
                                                               "ggplot2", 
                                                               "ggiraph", 
                                                               "gtable", 
                                                               "htmltools", 
                                                               "lubridate", 
                                                               "pool", 
                                                               "purrr", 
                                                               "scales", 
                                                               "shiny", 
                                                               "shinydashboard", 
                                                               "shinyWidgets", 
                                                               "stringr", 
                                                               "tidyr",
                                                               "CirceR"))

