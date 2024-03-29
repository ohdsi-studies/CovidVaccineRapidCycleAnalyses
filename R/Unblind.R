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

#' Unblind results
#'
#' @details
#' Unblind the results of analyses that have passed diagnostics.
#'
#' @param verifyDependencies   Check whether correct package versions are installed?
#' @param outputFolder         Name of local folder containing the results; make sure to use forward slashes
#'                             (/). Do not use a folder on a network drive since this greatly impacts
#'                             performance.
#' @param databaseId           A short string for identifying the database (e.g.
#'                             'Synpuf').
#' @param unblindCohortMethod      Unblind the cohort method analyses?
#' @param unblindSccs              Unblind the SCCS analyses?
#'
#' @export
unblind <- function(verifyDependencies = TRUE,
                    outputFolder,
                    databaseId,
                    unblindCohortMethod = TRUE,
                    unblindSccs = TRUE) {
  ParallelLogger::addDefaultFileLogger(file.path(outputFolder, "log.txt"))
  ParallelLogger::addDefaultErrorReportLogger(file.path(outputFolder, "errorReportR.txt"))
  on.exit(ParallelLogger::unregisterLogger("DEFAULT_FILE_LOGGER", silent = TRUE))
  on.exit(ParallelLogger::unregisterLogger("DEFAULT_ERRORREPORT_LOGGER", silent = TRUE), add = TRUE)

  if (verifyDependencies) {
    message("Checking whether correct package versions are installed")
    verifyDependencies()
  }

  if (unblindCohortMethod) {
    message("Unblinding cohort method results")
    exportCmResults(outputFolder, databaseId)
  }

  if (unblindSccs) {
    message("Unblinding SCCS")
    exportSccsResults(outputFolder, databaseId)
  }
  invisible(NULL)
}
