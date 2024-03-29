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

createCohorts <- function(connectionDetails,
                          cdmDatabaseSchema,
                          cohortDatabaseSchema,
                          cohortTable = "cohort",
                          tempEmulationSchema,
                          outputFolder) {
  if (!file.exists(outputFolder)) {
    dir.create(outputFolder)
  }

  connection <- DatabaseConnector::connect(connectionDetails)
  on.exit(DatabaseConnector::disconnect(connection))

  # Using oracleTempSchema because ROhdsiWebApi-generated code will still use that:
  .createCohorts(
    connection = connection,
    cdmDatabaseSchema = cdmDatabaseSchema,
    cohortDatabaseSchema = cohortDatabaseSchema,
    cohortTable = cohortTable,
    oracleTempSchema = tempEmulationSchema,
    outputFolder = outputFolder
  )

  pathToCsv <- system.file("settings", "NegativeControls.csv", package = "RcaCovidVaccine")
  negativeControls <- read.csv(pathToCsv)

  message("Creating negative control outcome cohorts")
  sql <- SqlRender::loadRenderTranslateSql(
    sqlFilename = "NegativeControlOutcomes.sql",
    packageName = "RcaCovidVaccine",
    dbms = connectionDetails$dbms,
    tempEmulationSchema = tempEmulationSchema,
    cdm_database_schema = cdmDatabaseSchema,
    cohort_database_schema = cohortDatabaseSchema,
    cohort_table = cohortTable,
    outcome_ids = unique(negativeControls$outcomeId)
  )
  DatabaseConnector::executeSql(connection, sql)

  # Check number of subjects per cohort:
  message("Counting cohorts")
  sql <- SqlRender::loadRenderTranslateSql(
    sqlFilename = "GetCounts.sql",
    packageName = "RcaCovidVaccine",
    dbms = connectionDetails$dbms,
    tempEmulationSchema = tempEmulationSchema,
    cdm_database_schema = cdmDatabaseSchema,
    work_database_schema = cohortDatabaseSchema,
    study_cohort_table = cohortTable
  )
  counts <- DatabaseConnector::querySql(connection, sql)
  colnames(counts) <- SqlRender::snakeCaseToCamelCase(colnames(counts))
  counts <- addCohortNames(counts)
  write.csv(counts, file.path(outputFolder, "CohortCounts.csv"), row.names = FALSE)
}
