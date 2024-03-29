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

#' Execute the Study
#'
#' @details
#' This function executes the RcaCovidVaccine Study.
#'
#' @param connectionDetails    An object of type \code{connectionDetails} as created using the
#'                             \code{\link[DatabaseConnector]{createConnectionDetails}} function in the
#'                             DatabaseConnector package.
#' @param cdmDatabaseSchema    Schema name where your patient-level data in OMOP CDM format resides.
#'                             Note that for SQL Server, this should include both the database and
#'                             schema name, for example 'cdm_data.dbo'.
#' @param cohortDatabaseSchema Schema name where intermediate data can be stored. You will need to have
#'                             write privilege in this schema. Note that for SQL Server, this should
#'                             include both the database and schema name, for example 'cdm_data.dbo'.
#' @param cohortTable          The name of the table that will be created in the work database schema.
#'                             This table will hold the exposure and outcome cohorts used in this
#'                             study.
#' @param tempEmulationSchema Some database platforms like Oracle and Impala do not truly support temp tables. To
#'                            emulate temp tables, provide a schema with write privileges where temp tables
#'                            can be created.
#' @param verifyDependencies   Check whether correct package versions are installed?
#' @param outputFolder         Name of local folder to place results; make sure to use forward slashes
#'                             (/). Do not use a folder on a network drive since this greatly impacts
#'                             performance.
#' @param databaseId           A short string for identifying the database (e.g.
#'                             'Synpuf').
#' @param databaseName         The full name of the database (e.g. 'Medicare Claims
#'                             Synthetic Public Use Files (SynPUFs)').
#' @param databaseDescription  A short description (several sentences) of the database.
#' @param createCohorts        Create the cohortTable table with the exposure and outcome cohorts?
#' @param runCohortDiagnostics Perform cohort diagnostics?
#' @param runCohortMethod      Perform the cohort method analyses?
#' @param runSccs              Perform the SCCS analyses?
#' @param generateCohortMethodDiagnostics  Generate diagnostics for the cohort method analyses?
#' @param generatSccsDiagnostics  Generate diagnostics for the SCCS analyses?
#' @param maxCores             How many parallel cores should be used? If more cores are made available
#'                             this can speed up the analyses.
#' @param minCellCount         The minimum number of subjects contributing to a count before it can be included
#'                             in packaged results.
#'
#' @export
execute <- function(connectionDetails,
                    cdmDatabaseSchema,
                    cohortDatabaseSchema = cdmDatabaseSchema,
                    cohortTable = "cohort",
                    tempEmulationSchema = getOption("sqlRenderTempEmulationSchema"),
                    verifyDependencies = TRUE,
                    outputFolder,
                    databaseId = "Unknown",
                    databaseName = "Unknown",
                    databaseDescription = "Unknown",
                    createCohorts = TRUE,
                    runCohortDiagnostics = TRUE,
                    runCohortMethod = TRUE,
                    runSccs = TRUE,
                    generateCohortMethodDiagnostics = TRUE,
                    generatSccsDiagnostics = TRUE,
                    maxCores = 4,
                    minCellCount = 5) {
  if (!file.exists(outputFolder)) {
    dir.create(outputFolder, recursive = TRUE)
  }

  ParallelLogger::addDefaultFileLogger(file.path(outputFolder, "log.txt"))
  ParallelLogger::addDefaultErrorReportLogger(file.path(outputFolder, "errorReportR.txt"))
  on.exit(ParallelLogger::unregisterLogger("DEFAULT_FILE_LOGGER", silent = TRUE))
  on.exit(ParallelLogger::unregisterLogger("DEFAULT_ERRORREPORT_LOGGER", silent = TRUE), add = TRUE)

  if (verifyDependencies) {
    message("Checking whether correct package versions are installed")
    verifyDependencies()
  }

  if (createCohorts) {
    message("Creating exposure and outcome cohorts")
    createCohorts(
      connectionDetails = connectionDetails,
      cdmDatabaseSchema = cdmDatabaseSchema,
      cohortDatabaseSchema = cohortDatabaseSchema,
      cohortTable = cohortTable,
      tempEmulationSchema = tempEmulationSchema,
      outputFolder = outputFolder
    )
  }

  if (runCohortDiagnostics) {
    message("Running cohort diagnostics")
    exportFolder <- file.path(outputFolder, "cohortDiagnostics")
    # if you want - we could do the switch to executeDiagnostics
    # https://github.com/OHDSI/SkeletonCohortDiagnosticsStudy/blob/main/R/CohortDiagnostics.R#L127-L168
    # if yes, it will need cohortDefinitionSet object etc as in this script

    CohortDiagnostics::runCohortDiagnostics(
      packageName = "RcaCovidVaccine",
      cohortToCreateFile = "settings/CohortsToCreate.csv",
      exportFolder = exportFolder,
      connectionDetails = connectionDetails,
      cdmDatabaseSchema = cdmDatabaseSchema,
      cohortDatabaseSchema = cohortDatabaseSchema,
      cohortTable = cohortTable,
      tempEmulationSchema = tempEmulationSchema,
      databaseId = databaseId,
      databaseName = databaseName,
      databaseDescription = databaseDescription,
      runInclusionStatistics = FALSE,
      runIncludedSourceConcepts = TRUE,
      runOrphanConcepts = TRUE,
      runTimeDistributions = TRUE,
      runVisitContext = TRUE,
      runBreakdownIndexEvents = TRUE,
      runIncidenceRate = TRUE,
      runTimeSeries = FALSE,
      runCohortOverlap = TRUE,
      runCohortCharacterization = TRUE,
      minCellCount = minCellCount
    )
    CohortDiagnostics::preMergeDiagnosticsFiles(exportFolder)
  }

  if (runCohortMethod) {
    message("Running CohortMethod analyses")
    runCohortMethod(
      connectionDetails = connectionDetails,
      cdmDatabaseSchema = cdmDatabaseSchema,
      cohortDatabaseSchema = cohortDatabaseSchema,
      cohortTable = cohortTable,
      tempEmulationSchema = tempEmulationSchema,
      outputFolder = outputFolder,
      maxCores = maxCores
    )
  }

  if (runSccs) {
    message("Running SCCS analyses")
    runSelfControlledCaseSeries(
      connectionDetails = connectionDetails,
      cdmDatabaseSchema = cdmDatabaseSchema,
      cohortDatabaseSchema = cohortDatabaseSchema,
      cohortTable = cohortTable,
      outputFolder = outputFolder,
      maxCores = maxCores
    )
  }

  if (generateCohortMethodDiagnostics) {
    message("Generating CohortMethod diagnostics")
    generateCmDiagnostics(outputFolder = outputFolder)
  }

  if (generatSccsDiagnostics) {
    message("Generating SCCS diagnostics")
    generateSccsDiagnostics(outputFolder = outputFolder, maxCores)
  }

  invisible(NULL)
}
