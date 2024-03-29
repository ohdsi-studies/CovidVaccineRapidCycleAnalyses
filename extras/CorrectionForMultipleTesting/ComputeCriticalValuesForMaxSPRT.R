library(dplyr)

familyWiseAlpha <- 0.5
maxCores <- 10

resultFolder <- "../Results"

dbFolders <- list.files(resultFolder, include.dirs = TRUE)
dbFolders <- dbFolders[dir.exists(file.path(resultFolder, dbFolders))]

# Bonferroni -----------------------------------------
specs <- openxlsx::read.xlsx("inst/settings/HVRCAAnalysesSpecs.xlsx", colNames = FALSE)
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
discoveryAnalyses<- discoveryAnalyses %>%
  mutate(nAnalyses = ifelse(method == "CohortMethod", 2, 1))
numberOfTests <- length(dbFolders) * sum(discoveryAnalyses$nAnalyses) 

alpha <- familyWiseAlpha / numberOfTests

writeLines(sprintf("Number of test in a look: %d. alpha per test: %0.5f", numberOfTests, alpha))

# Cohort method ------------------------------------------
allCmSampleSizes <- data.frame()
for (dbFolder in dbFolders) {
  mdrr <- read.csv(file.path(resultFolder, dbFolder, "cmDiagnostics", "mdrr.csv"))
  sampleSizes <- mdrr %>%
    inner_join(select(discoveryAnalyses, analysisId, outcomeId), by = c("analysisId", "outcomeId")) %>%
    mutate(z = comparatorDays / targetDays,
           events = totalOutcomes) %>%
    select(analysisId, outcomeId, targetId, comparatorId, z, events) %>%
    mutate(database = basename(dbFolder),
           look = 1,
           method = "Cohort method")
  allCmSampleSizes <- rbind(allCmSampleSizes, sampleSizes)
}

# 2% more subjects at each look:
allCmSampleSizes <- bind_rows(
  allCmSampleSizes,
  allCmSampleSizes %>%
    mutate(look = 2,
           events = events * 1.02),
  allCmSampleSizes %>%
    mutate(look = 3,
           events = events * 1.04),
  allCmSampleSizes %>%
    mutate(look = 4,
           events = events * 1.06)
)
write.csv(allCmSampleSizes, "extras/CorrectionForMultipleTesting/expectedPowerCohortMethod.csv", row.names = FALSE)

## Remove Health Verity, add Pharmetrics
allCmSampleSizes <- readr::read_csv("extras/CorrectionForMultipleTesting/expectedPowerCohortMethod.csv", show_col_types = FALSE)
allCmSampleSizes <- allCmSampleSizes %>%
  filter(!database %in% c("HealthVerity", "Pharmetrics"))

mdrr <- read.csv(file.path(resultFolder, "Pharmetrics", "cmDiagnostics", "mdrr.csv"))
sampleSizes <- mdrr %>%
  inner_join(select(discoveryAnalyses, analysisId, outcomeId), by = c("analysisId", "outcomeId")) %>%
  mutate(z = comparatorDays / targetDays,
         events = totalOutcomes) %>%
  select(analysisId, outcomeId, targetId, comparatorId, z, events) %>%
  mutate(database = "Pharmetrics",
         look = 1,
         method = "Cohort method")
allCmSampleSizes <- rbind(allCmSampleSizes, sampleSizes)

write.csv(allCmSampleSizes, "extras/CorrectionForMultipleTesting/expectedPowerCohortMethod.csv", row.names = FALSE)

# SCCS ------------------------------------------
allSccsSampleSizes <- data.frame()
for (dbFolder in dbFolders) {
  mdrr <- read.csv(file.path(resultFolder, dbFolder, "sccsDiagnostics", "mdrr.csv"))
  sampleSizes <- mdrr %>%
    filter(exposureId == 1082) %>%
    inner_join(select(discoveryAnalyses, analysisId, outcomeId), by = c("analysisId", "outcomeId")) %>%
    mutate(z = 1/propTimeExposed - 1) %>%
    select(analysisId, outcomeId, exposureId, z, events) %>%
    mutate(database = basename(dbFolder),
           look = 1,
           method = "SCCS")
  allSccsSampleSizes <- rbind(allSccsSampleSizes, sampleSizes)
}
# 2% more subjects, and 3 months more time per subject, at each look: 
allSccsSampleSizes <- bind_rows(
  allSccsSampleSizes,
  allSccsSampleSizes %>%
    mutate(look = 2,
           events = case_when(
             analysisId %in% 1:5 ~ events * (15/12) * 1.02,
             analysisId %in% 6:10 ~ events * (39/36) * 1.02,
             TRUE ~ events * (12/9) * 1.02)),
  allSccsSampleSizes %>%
    mutate(look = 3,
           events = case_when(
             analysisId %in% 1:5 ~ events * (18/12) * 1.04,
             analysisId %in% 6:10 ~ events * (42/36) * 1.04,
             TRUE ~ events * (15/9) * 1.04)),
  allSccsSampleSizes %>%
    mutate(look = 4,
           events = case_when(
             analysisId %in% 1:5 ~ events * (21/12) * 1.06,
             analysisId %in% 6:10 ~ events * (45/36) * 1.06,
             TRUE ~ events * (18/9) * 1.06))
)
write.csv(allSccsSampleSizes, "extras/CorrectionForMultipleTesting/expectedPowerSccs.csv", row.names = FALSE)

## Remove Health Verity, add Pharmetrics
allSccsSampleSizes <- readr::read_csv("extras/CorrectionForMultipleTesting/expectedPowerSccs.csv", show_col_types = FALSE)
allSccsSampleSizes <- allSccsSampleSizes %>%
  filter(!database %in% c("HealthVerity", "Pharmetrics"))

mdrr <- read.csv(file.path(resultFolder, "Pharmetrics", "sccsDiagnostics", "mdrr.csv"))
sampleSizes <- mdrr %>%
  filter(exposureId == 1082) %>%
  inner_join(select(discoveryAnalyses, analysisId, outcomeId), by = c("analysisId", "outcomeId")) %>%
  mutate(z = 1/propTimeExposed - 1) %>%
  select(analysisId, outcomeId, exposureId, z, events) %>%
  mutate(database = "Pharmetrics",
         look = 1,
         method = "SCCS")
allSccsSampleSizes <- rbind(allSccsSampleSizes, sampleSizes)

write.csv(allSccsSampleSizes, "extras/CorrectionForMultipleTesting/expectedPowerSccs.csv", row.names = FALSE)

# Compute critical values ------------------------------------------------------

# subset <- subset[[1]]
# subset <- allSccsSampleSizes[allSccsSampleSizes$analysisId == 13 &
#                                allSccsSampleSizes$outcomeId == 5831 &
#                                allSccsSampleSizes$database == "CCAE", ]
computeCv <- function(subset, alpha) {
  if (all(is.na(subset$events))) {
    sampleSizeUpperLimit <- 0
  } else {
    sampleSizeUpperLimit <- max(subset$events, na.rm = TRUE)
  }
  events <- subset$events
  if (length(events) > 1) {
    events[2:length(events)] <- events[2:length(events)] - events[1:(length(events)-1)]
    events <- events[events != 0]
  }
  if (sampleSizeUpperLimit == 0 || length(events) == 0) {
    cv <- Inf
    cvAlpha <- 0
  } else {
    suppressMessages(
      cv <- EmpiricalCalibration::computeCvBinomial(groupSizes = events,
                                                    z = mean(subset$z),
                                                    minimumEvents = 1,
                                                    sampleSize = 1e7,
                                                    alpha = alpha)
    )
    cvAlpha <- attr(cv, "alpha")
  }
  row <- subset %>%
    head(1) %>%
    select(-z, -events, -look) %>%
    mutate(totalExpectedEvents = sampleSizeUpperLimit,
           alpha = alpha,
           cv = cv,
           cvAlpha = cvAlpha)
  return(row)  
}

cluster <- ParallelLogger::makeCluster(maxCores)
ParallelLogger::clusterRequire(cluster, "dplyr")
subsets <- split(allCmSampleSizes, 
                 paste(allCmSampleSizes$analysisId, allCmSampleSizes$targetId, allCmSampleSizes$comparatorId, allCmSampleSizes$outcomeId, allCmSampleSizes$database))
cvsCohortMethod <- ParallelLogger::clusterApply(cluster, subsets, computeCv, alpha = alpha)
cvsCohortMethod <- bind_rows(cvsCohortMethod)

subsets <- split(allSccsSampleSizes, 
                 paste(allSccsSampleSizes$analysisId, allSccsSampleSizes$exposureId, allSccsSampleSizes$outcomeId, allSccsSampleSizes$database))
cvsSccs <- ParallelLogger::clusterApply(cluster, subsets, computeCv, alpha = alpha)
cvsSccs <- bind_rows(cvsSccs)

ParallelLogger::stopCluster(cluster)

write.csv(cvsCohortMethod, "inst/settings/criticalValuesCohortMethod.csv", row.names = FALSE)
write.csv(cvsSccs, "inst/settings/criticalValuesSccs.csv", row.names = FALSE)

writeLines(sprintf("Effective family-wise alpha: %0.3f", sum(cvsCohortMethod$cvAlpha) + sum(cvsSccs$cvAlpha)))

