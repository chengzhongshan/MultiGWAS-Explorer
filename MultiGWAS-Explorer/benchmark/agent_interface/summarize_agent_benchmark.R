d <- read.delim("agent_interface_runs.tsv", stringsAsFactors = FALSE)
rows <- lapply(split(d, d$path), function(x) {
  q <- quantile(x$elapsed_seconds, c(0, .25, .5, .75, 1), names = FALSE, type = 7)
  data.frame(
    PATH = x$path[1], RUNS = nrow(x), SUCCESS_RATE = mean(x$success),
    PARITY_RATE = mean(x$parity_with_cli), UNIQUE_ARTIFACT_HASHES = length(unique(x$artifact_set_sha256)),
    ELAPSED_MIN = q[1], ELAPSED_Q1 = q[2], ELAPSED_MEDIAN = q[3],
    ELAPSED_Q3 = q[4], ELAPSED_MAX = q[5]
  )
})
out <- do.call(rbind, rows)
write.table(out, "agent_interface_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

r <- read.delim("failure_recovery.tsv", stringsAsFactors = FALSE)
recovery <- data.frame(
  INJECTED_FAILURE_RUNS = nrow(r),
  FAILURE_DETECTION_RATE = mean(r$injected_failure_detected),
  CORRECTED_RETRY_RECOVERY_RATE = mean(r$corrected_retry_success),
  RECOVERY_ELAPSED_MEDIAN = median(r$recovery_seconds),
  RECOVERY_ELAPSED_Q1 = unname(quantile(r$recovery_seconds, .25)),
  RECOVERY_ELAPSED_Q3 = unname(quantile(r$recovery_seconds, .75))
)
write.table(recovery, "failure_recovery_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
