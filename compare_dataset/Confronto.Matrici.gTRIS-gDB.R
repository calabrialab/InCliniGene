### nuova FE Luca - 1 giugno 2021
library(reshape2)
library(ggplot2)
library(gplots)
library(openxlsx)
library(psych)

###############################################################
## functions
###############################################################
# source("/Users/calabria.andrea/Dropbox (FONDAZIONE TELETHON)/TIGET/Workbench/isatk/script/R/isa_utils_functions.R")
reannotate_Repeat_IS <- function(df, id_cols = c("Chr", "Pos", "Strand")) {
  # slice matrix with only repeats
  slice_R <- df[which(df[id_cols[1]] == "R" & df[id_cols[2]] == -1),]
  slice_notR <- df[which( !(df[id_cols[1]] == "R" & df[id_cols[2]] == -1) ),]
  slice_R[id_cols[2]] <- seq(1,nrow(slice_R))
  return( rbind(slice_notR, slice_R) )
}
###############################################################

id_cols <- c("chr", "integration_locus", "strand", "GeneName", "GeneStrand")
id_cols <- c("Chr", "Pos", "Strand")

#### ------------- on VA -------------- #####
# gtris
gtrisf <- paste0("data/gtris.clu_ClonalExpansion_extended_repeat.csv.gz", sep = "") 
gtris <- read.csv(gtrisf, header=TRUE, fill=T, sep=',', check.names = FALSE, na.strings = c("NONE", "NA", "NULL", "NaN", "ND", ""), row.names = 1)
names(gtris) <- gsub("ClonalExpansion_unique_r1_", "", colnames(gtris))
rownames(gtris) <- apply(gtris[id_cols], 1, function(x) { paste(x[1],"_",as.numeric(x[2]),"_",x[3], sep="") })

gdbf <- paste0("data/gdb.TableFile_clu.csv", sep = "") 
gdb <- read.csv(gdbf, header=TRUE, fill=T, sep=',', check.names = FALSE, na.strings = c("NONE", "NA", "NULL", "NaN", "ND", ""), row.names = 1)

gdb_reannotated <- reannotate_Repeat_IS(df = gdb, id_cols = id_cols)
# gdb <- read.xlsx(xlsxFile = gdbf, check.names = F)
colnames(gdb)
gdb_reannotated$Chr <- gsub("R", "chrR", gdb_reannotated$Chr)
gdb_reannotated$Strand <- gsub("R", "+", gdb_reannotated$Strand)
rownames(gdb_reannotated) <- apply(gdb_reannotated[id_cols], 1, function(x) { paste(x[1],"_",as.numeric(x[2]),"_",x[3], sep="") })
# gdb$Pos <- ifelse(gdb$Pos == -1, gdb$Pos == gdb$ID, gdb$Pos)

gdb <- gdb_reannotated
# intersect on cols -> stats
common_samples <- intersect(colnames(gtris), colnames(gdb))
venn(list("gTRIS" = rownames(gtris[which(!(gtris$Chr %in% "chrR")),]), "gDB"= rownames(gdb[which(!(gdb$Chr %in% "chrR")),])))
venn(list("gTRIS" = colnames(gtris), "gDB"= colnames(gdb)))

id_vars_melt <- c("Chr", "Pos", "Strand")

gdb[gdb==0] <- NA
gdb_stats <- as.data.frame(describe(gdb[common_samples]))
gdb_melt <- melt(data = gdb[common_samples], id.vars = id_vars_melt, variable.name = "Sample", value.name = "Abundance", na.rm = T)
gdb_melt_allchr <- as.data.frame(table(gdb_melt$Chr, gdb_melt$Sample))
names(gdb_melt_allchr) <- c("chr", "Sample", "nIS_byChr")
gdb_melt_chrR <- gdb_melt_allchr[which(gdb_melt_allchr$chr == "chrR"),]
rownames(gdb_melt_chrR) <- gdb_melt_chrR$Sample
gdb_stats <- merge(x = gdb_stats, y = gdb_melt_chrR, by = 0)
rownames(gdb_stats) <- gdb_stats$Row.names
gdb_stats <- gdb_stats[, !(names(gdb_stats) %in% c("Row.names"))]
gdb_sums <- data.frame(
  "NumIS" = apply(gdb[setdiff(common_samples, id_vars_melt)], 2, function(x) {length(x[!is.na(x)])}),
  "SeqCountSum" = apply(gdb[setdiff(common_samples, id_vars_melt)], 2, function(x) {sum(x, na.rm = T)})
)
gdb_stats <- merge(x = gdb_stats, y = gdb_sums, by = 0)
rownames(gdb_stats) <- gdb_stats$Row.names
gdb_stats <- gdb_stats[, !(names(gdb_stats) %in% c("Row.names"))]
names(gdb_stats) <- paste("gDB", colnames(gdb_stats))

gtris[gtris==0] <- NA
gtris_stats <- as.data.frame(describe(gtris[common_samples]))
gtris_melt <- melt(data = gtris[common_samples], id.vars = id_vars_melt, variable.name = "Sample", value.name = "Abundance", na.rm = T)
gtris_melt_allchr <- as.data.frame(table(gtris_melt$Chr, gtris_melt$Sample))
names(gtris_melt_allchr) <- c("chr", "Sample", "nIS_byChr")
gtris_melt_chrR <- gtris_melt_allchr[which(gtris_melt_allchr$chr == "chrR"),]
rownames(gtris_melt_chrR) <- gtris_melt_chrR$Sample
gtris_stats <- merge(x = gtris_stats, y = gtris_melt_chrR, by = 0)
rownames(gtris_stats) <- gtris_stats$Row.names
gtris_stats <- gtris_stats[, !(names(gtris_stats) %in% c("Row.names"))]
gtris_sums <- data.frame(
  "NumIS" = apply(gtris[setdiff(common_samples, id_vars_melt)], 2, function(x) {length(x[!is.na(x)])}),
  "SeqCountSum" = apply(gtris[setdiff(common_samples, id_vars_melt)], 2, function(x) {sum(x, na.rm = T)})
)
gtris_stats <- merge(x = gtris_stats, y = gtris_sums, by = 0)
rownames(gtris_stats) <- gtris_stats$Row.names
gtris_stats <- gtris_stats[, !(names(gtris_stats) %in% c("Row.names"))]
names(gtris_stats) <- paste("gTRIS", colnames(gtris_stats))

## merge
merge_and_compare <- merge(x = gdb_stats, y = gtris_stats, by = 0, all.x = T)
rownames(merge_and_compare) <- merge_and_compare$Row.names
merge_and_compare <- merge_and_compare[, !(names(merge_and_compare) %in% c("Row.names"))]
# add some results
merge_and_compare$delta_nIS <- merge_and_compare$`gDB n` - merge_and_compare$`gTRIS n`
merge_and_compare$delta_nIS_on_gDB <- merge_and_compare$delta_nIS/merge_and_compare$`gDB n`

merge_and_compare$delta_SeqCountSUM <- merge_and_compare$`gDB SeqCountSum` - merge_and_compare$`gTRIS SeqCountSum`
merge_and_compare$delta_SeqCountSUM_on_gTRIS <- merge_and_compare$delta_SeqCountSUM / merge_and_compare$`gTRIS SeqCountSum`

merge_and_compare$delta_nIS_chrR <- merge_and_compare$`gDB nIS_byChr` - merge_and_compare$`gTRIS nIS_byChr`
merge_and_compare$delta_nIS_chrR_on_gDB <- merge_and_compare$delta_nIS_chrR / merge_and_compare$`gDB nIS_byChr`
merge_and_compare$delta_ab_max <- merge_and_compare$`gDB max` - merge_and_compare$`gTRIS max`
merge_and_compare$ratio_skew <- merge_and_compare$`gDB skew` / merge_and_compare$`gTRIS skew`
merge_and_compare$ratio_kurt <- merge_and_compare$`gDB kurtosis` / merge_and_compare$`gTRIS kurtosis`

write.xlsx(x = merge_and_compare, file = "data/Comparison.xlsx")

# specific tests on samples and chrR
View(gdb[c(id_vars_melt, "ID00000000000000027038", "SUBG")])
View(gtris[c(id_vars_melt, "ID00000000000000027038")])

View(gdb[c(id_vars_melt, "ID00000000000000027002")])
View(gtris[c(id_vars_melt, "ID00000000000000027002")])

View(gdb[c(id_vars_melt, "ID00000000000000026998")])
View(gtris[c(id_vars_melt, "ID00000000000000026998")])

