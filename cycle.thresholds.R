source("packages.R")

well.pattern <- paste0(
  "(?<letter>[A-P])",
  "(?<number>[0-9]+)")

csv.path.vec <- Sys.glob("cycle.thresholds/*.csv")
cycle.thresholds.list <- list()
for(csv.path.i in seq_along(csv.path.vec)){
  csv.path <- csv.path.vec[[csv.path.i]]
  csv.line.vec <- readLines(csv.path)
  print(header.row <- grep("^Well", csv.line.vec))
  ct.tab <- fread(
    csv.path, skip=header.row-1, header=TRUE,
    na.strings="Undetermined",
    select=c("CT", "Well Position"))
  well.mat <- str_match_named(ct.tab$Well, well.pattern, list(
    number=as.integer))
  cycle.thresholds.list[[csv.path]] <- data.table(
    csv.path.i, csv.path, well.mat, ct.tab,
    number.fac=factor(well.mat[, "number"]))
}
cycle.thresholds <- do.call(rbind, cycle.thresholds.list)

save(cycle.thresholds, file="cycle.thresholds.RData")
