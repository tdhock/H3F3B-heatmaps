source("packages.R")

well.pattern <- paste0(
  "(?<letter>[A-P])",
  "(?<number>[0-9]+)")

txt.path.vec <- Sys.glob("cycle.thresholds/*/*.txt")
cycle.thresholds.list <- list()
for(txt.path.i in seq_along(txt.path.vec)){
  txt.path <- txt.path.vec[[txt.path.i]]
  csv.line.vec <- readLines(txt.path)
  print(header.rows <- grep("^Well", csv.line.vec))
  data.line.vec <- grep("^[0-9]+\t", csv.line.vec, value=TRUE)
  first.header.row <- header.rows[1]
  header.data.line.vec <- c(csv.line.vec[first.header.row], data.line.vec)
  cleaned.path <- sub("txt$", "cleaned.csv", txt.path)
  writeLines(header.data.line.vec, cleaned.path)
  ct.tab <- fread(
    cleaned.path, header=TRUE,
    ##nrows=last.table.row-header.row,
    na.strings="Undetermined",
    select=c("Ct", "Sample Name"))
  well.mat <- str_match_named(ct.tab[["Sample Name"]], well.pattern, list(
    number=as.integer))
  only.data <- data.table(
    txt.path.i, txt.path, well.mat, ct.tab,
    number.fac=factor(well.mat[, "number"], 1:24),
    letter.fac=factor(well.mat[, "letter"], LETTERS[1:18]))
  data.dir <- dirname(txt.path)
  row.csv <- Sys.glob(file.path(data.dir, "*_Row.csv"))
  row.info <- fread(row.csv)
  table(row.info[["Plate Row"]])
  table(only.data[["letter"]])
  setkey(row.info, `Plate Row`)
  setkey(only.data, letter)
  data.row.info <- only.data[row.info]
  col.csv <- Sys.glob(file.path(data.dir, "*_Column.csv"))
  col.info <- fread(col.csv)
  table(data.row.info[["number"]])
  table(col.info[["Column"]])
  setkey(col.info, Column)
  setkey(data.row.info, number)
  data.row.col.info <- data.row.info[col.info]
  rbind(
    data=nrow(ct.tab),
    rows=nrow(row.info),
    data.rows=nrow(data.row.info),
    cols=nrow(col.info),
    data.rows.cols=nrow(data.row.col.info))
  cycle.thresholds.list[[paste("regular", txt.path)]] <- data.row.col.info
  irregular.csv <- Sys.glob(file.path(data.dir, "*_Irregular.csv"))
  if(length(irregular.csv)){
    irregular.info <- fread(irregular.csv)
    setkey(irregular.info, `Plate Row`, Column)
    setkey(only.data, letter, number)
    cycle.thresholds.list[[paste("irregular", txt.path)]] <-
      only.data[irregular.info]
  }
}

(names.list <- lapply(cycle.thresholds.list, names))
name.counts <- table(unlist(names.list))
stopifnot(all(name.counts==length(names.list)))
name.counts[order(name.counts)]

ordered.columns.list <- list()
for(data.name in names(cycle.thresholds.list)){
  one.dt <- cycle.thresholds.list[[data.name]]
  ordered.columns.list[[data.name]] <-
    one.dt[, names(name.counts), with=FALSE]
}
no.stats <- do.call(rbind, ordered.columns.list)

gapdh.means <- no.stats[Gene=="gapdh", list(
  gapdh.mean.Ct=mean(Ct)),
  by=.(Genotype, Time, Treatment)]
setkey(gapdh.means, Genotype, Time, Treatment)
setkey(no.stats, Genotype, Time, Treatment)
no.gapdh <- no.stats[Gene!="gapdh",]
cycle.thresholds <- no.gapdh[gapdh.means]
cycle.thresholds[, Ct.diff := Ct - gapdh.mean.Ct]
cycle.thresholds[, relative.expression := gapdh.mean.Ct - Ct]
cycle.thresholds[, fold.difference := 2^(-Ct.diff)]
cycle.thresholds[, log.fold.difference := log(fold.difference)]

save(cycle.thresholds, file="cycle.thresholds.RData")
