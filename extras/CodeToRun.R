library(RcaCovidVaccine)

# Optional: specify where the temporary files (used by the Andromeda package) will be created:
options(andromedaTempFolder = "d:/andromedaTemp")

# Maximum number of cores to be used:
maxCores <- parallel::detectCores()

# For some database platforms (e.g. Oracle): define a schema that can be used to emulate temp tables:
options(sqlRenderTempEmulationSchema = NULL)

# Settings for Pharmetrics 
outputFolder <- "d:/RcaCovidVaccine_Pharmetrics"
connectionDetails <- DatabaseConnector::createConnectionDetails(dbms = "redshift",
                                                                connectionString = keyring::key_get("redShiftConnectionStringOhdaPharmetrics"),
                                                                user = keyring::key_get("redShiftUserName"),
                                                                password = keyring::key_get("redShiftPassword"))
cdmDatabaseSchema <- "cdm_iqvia_pharmetrics_plus_v2421"
cohortDatabaseSchema <- "scratch_mschuemi"
cohortTable <- "cohort_epi_1071_pharmetrics"
databaseId <- "Pharmetrics"
databaseName <- "Pharmetrics Plus"
databaseDescription <- "Data is from 2006 - 2022 and comprises of fully adjudicated medical and pharmacy claims. It contains a longitudinal view of inpatient and outpatient services, prescription and office/outpatient administered drugs, costs and enrollment information. With IQVIA Adjudicated Health Plan Claims, an enrolled patient can be tracked across all sites of care: hospital, specialist, emergency room, pharmacy, primary care, and more."

# Settings for CCAE
outputFolder <- "d:/RcaCovidVaccine_ccae"
connectionDetails <- DatabaseConnector::createConnectionDetails(dbms = "redshift",
                                                                connectionString = keyring::key_get("redShiftConnectionStringOhdaCcae"),
                                                                user = keyring::key_get("redShiftUserName"),
                                                                password = keyring::key_get("redShiftPassword"))
cohortDatabaseSchema <- "scratch_mschuemi"
cohortTable <- "cohort_epi_1071_ccae"
cdmDatabaseSchema <- "cdm_truven_ccae_v2435"
databaseId <- "CCAE"
databaseName <- "IBM MarketScan Commercial Claims and Encounters Database"
databaseDescription <- "IBM MarketScan® Commercial Claims and Encounters Database (CCAE) represent data from individuals enrolled in United States employer-sponsored insurance health plans. The data includes adjudicated health insurance claims (e.g. inpatient, outpatient, and outpatient pharmacy) as well as enrollment data from large employers and health plans who provide private healthcare coverage to employees, their spouses, and dependents. Additionally, it captures laboratory tests for a subset of the covered lives. This administrative claims database includes a variety of fee-for-service, preferred provider organizations, and capitated health plans."

# Settings for Optum
outputFolder <- "d:/RcaCovidVaccine_optum"
connectionDetails <- DatabaseConnector::createConnectionDetails(dbms = "redshift",
                                                                connectionString = keyring::key_get("redShiftConnectionStringOhdaOptumSes"),
                                                                user = keyring::key_get("redShiftUserName"),
                                                                password = keyring::key_get("redShiftPassword"))
cohortDatabaseSchema <- "scratch_mschuemi"
cohortTable <- "cohort_epi_1071_optum"
cdmDatabaseSchema <- "cdm_optum_extended_ses_v2437"
databaseId <- "Optum"
databaseName <- "Optum Clinformatics Extended Data Mart - Social-economic Status (SeS)"
databaseDescription <- "Optum Clinformatics Extended DataMart is an adjudicated US administrative health claims database for members of private health insurance, who are fully insured in commercial plans or in administrative services only (ASOs), Legacy Medicare Choice Lives (prior to January 2006), and Medicare Advantage (Medicare Advantage Prescription Drug coverage starting January 2006).  The population is primarily representative of commercial claims patients (0-65 years old) with some Medicare (65+ years old) however ages are capped at 90 years.  It includes data captured from administrative claims processed from inpatient and outpatient medical services and prescriptions as dispensed, as well as results for outpatient lab tests processed by large national lab vendors who participate in data exchange with Optum."


execute(connectionDetails = connectionDetails,
        cdmDatabaseSchema = cdmDatabaseSchema,
        cohortDatabaseSchema = cohortDatabaseSchema,
        cohortTable = cohortTable,
        outputFolder = outputFolder,
        databaseId = databaseId,
        databaseName = databaseName,
        databaseDescription = databaseDescription,
        verifyDependencies = TRUE,
        createCohorts = TRUE,
        runCohortDiagnostics = F,
        runCohortMethod = TRUE,
        runSccs = TRUE,
        generateCohortMethodDiagnostics = TRUE,
        generatSccsDiagnostics = TRUE,
        maxCores = maxCores)

gitResultsFolder <- sprintf("../Results/%s", databaseId)
dir.create(gitResultsFolder, recursive = TRUE)
file.copy(file.path(outputFolder, "CmDiagnosticsOverview.csv"), gitResultsFolder, overwrite = TRUE)
file.copy(file.path(outputFolder, "cmDiagnostics"), gitResultsFolder, recursive = TRUE, overwrite = TRUE)
file.copy(file.path(outputFolder, "SccsDiagnosticsOverview.csv"), gitResultsFolder, overwrite = TRUE)
file.copy(file.path(outputFolder, "sccsDiagnostics"), gitResultsFolder, recursive = TRUE, overwrite = TRUE)
file.copy(file.path(outputFolder, "CohortCounts.csv"), gitResultsFolder, overwrite = TRUE)

# Make sure to check unblind decisions in 'CmDiagnosticsOverview.csv' and 'SccsDiagnosticsOverview.csv' 
# before proceeding.

unblind(verifyDependencies = TRUE,
        outputFolder = outputFolder,
        databaseId = databaseId,
        unblindCohortMethod = TRUE,
        unblindSccs = TRUE)

gitResultsFolder <- sprintf("../Results/%s", databaseId)
file.copy(file.path(outputFolder, "cmResults"), gitResultsFolder, recursive = TRUE, overwrite = TRUE)
file.copy(file.path(outputFolder, "sccsResults"), gitResultsFolder, recursive = TRUE, overwrite = TRUE)
