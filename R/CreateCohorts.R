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

.createCohorts <- function(connection,
                           cdmDatabaseSchema,
                           vocabularyDatabaseSchema = cdmDatabaseSchema,
                           cohortDatabaseSchema,
                           cohortTable,
                           oracleTempSchema,
                           outputFolder) {

  # Create study cohort table structure:
  sql <- SqlRender::loadRenderTranslateSql(
    sqlFilename = "CreateCohortTable.sql",
    packageName = "RcaCovidVaccine",
    dbms = attr(connection, "dbms"),
    oracleTempSchema = oracleTempSchema,
    cohort_database_schema = cohortDatabaseSchema,
    cohort_table = cohortTable
  )
  DatabaseConnector::executeSql(connection, sql, progressBar = FALSE, reportOverallTime = FALSE)



  # Instantiate cohorts:
  pathToCsv <- system.file("settings", "CohortsToCreate.csv", package = "RcaCovidVaccine")
  cohortsToCreate <- readr::read_csv(pathToCsv, col_types = readr::cols())

  # For some reason can't create these cohorts from R
  # badCohorts <- c(3599, 610, 624, 628, 635, 1068, 1073, 1074, 1075, 1077, 1078, 1079, 1308, 2927, 3072, 366, 2796)
  badCohorts <- c()

  for (i in 1:nrow(cohortsToCreate)) {
    if (cohortsToCreate$cohortId[i] %in% badCohorts) {
      writeLines(paste("Copying cohort:", cohortsToCreate$name[i]))
      resultsDatabaseSchema <- gsub("_cdm_", "_results_", cdmDatabaseSchema)
      sql <- "INSERT INTO @cohort_database_schema.@cohort_table (subject_id, cohort_definition_id, cohort_start_date, cohort_end_date)
      SELECT subject_id,
        cohort_definition_id,
        cohort_start_date,
        cohort_end_date
      FROM @results_database_schema.cohort
      WHERE cohort_definition_id = @cohort_id;"
      DatabaseConnector::renderTranslateExecuteSql(
        connection = connection,
        sql = sql,
        results_database_schema = resultsDatabaseSchema,
        cohort_database_schema = cohortDatabaseSchema,
        cohort_table = cohortTable,
        cohort_id = cohortsToCreate$cohortId[i]
      )
    } else {
      message(paste("Creating cohort:", cohortsToCreate$name[i]))
      sql <- SqlRender::loadRenderTranslateSql(
        sqlFilename = paste0(cohortsToCreate$name[i], ".sql"),
        packageName = "RcaCovidVaccine",
        dbms = attr(connection, "dbms"),
        oracleTempSchema = oracleTempSchema,
        cdm_database_schema = cdmDatabaseSchema,
        vocabulary_database_schema = vocabularyDatabaseSchema,
        target_database_schema = cohortDatabaseSchema,
        target_cohort_table = cohortTable,
        target_cohort_id = cohortsToCreate$cohortId[i]
      )
      if (connection@dbms %in% c("redshift", "postgresql")) {
        # Some performance tuning:
        sql <- gsub("CREATE TABLE #qualified_events", "ANALYZE #Codesets;\n\nCREATE TABLE #qualified_events", sql)
        sql <- gsub("CREATE TABLE #Inclusion_0", "ANALYZE #qualified_events;\n\nCREATE TABLE #Inclusion_0", sql)
      }
      DatabaseConnector::executeSql(connection, sql)
    }
  }

  # Fetch cohort counts:
  sql <- "SELECT cohort_definition_id, COUNT(*) AS count FROM @cohort_database_schema.@cohort_table GROUP BY cohort_definition_id"
  sql <- SqlRender::render(sql,
    cohort_database_schema = cohortDatabaseSchema,
    cohort_table = cohortTable
  )
  sql <- SqlRender::translate(sql, targetDialect = attr(connection, "dbms"))
  counts <- DatabaseConnector::querySql(connection, sql)
  names(counts) <- SqlRender::snakeCaseToCamelCase(names(counts))
  counts <- merge(counts, data.frame(
    cohortDefinitionId = cohortsToCreate$cohortId,
    cohortName = cohortsToCreate$name
  ))
  readr::write_csv(x = counts, file = file.path(outputFolder, "CohortCounts.csv"))
}
