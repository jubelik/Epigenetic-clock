cran_mirror <- Sys.getenv("R_CRAN_MIRROR", "http://cran.rstudio.com/")
options(repos = cran_mirror)

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("methylKit")

library(BiocManager)
library(stringi)
library(pacman)
library(methylKit)
library(ade4)
library(FactoMineR)
library(devtools)
library(factoextra)
library(tibble)
library(caret)

file.list <- list(
"./21-105_AHCY7CDSXC_S101_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./21-111_AHCY7CDSXC_S20_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./21-112_AHCY7CDSXC_S58_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./21-113_AHCY7CDSXC_S47_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./21-118_AHCY7CDSXC_S91_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./21-142_AHCY7CDSXC_S32_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./21-143_AHCY7CDSXC_S102_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./21-146_AHCY7CDSXC_S97_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./21-147_AHCY7CDSXC_S95_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./21-176_AHCY7CDSXC_S92_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./21-178_AHCY7CDSXC_S21_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./21-180_AHCY7CDSXC_S98_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./21-182_AHCY7CDSXC_S36_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./21-183_AHCY7CDSXC_S31_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./21-184_AHCY7CDSXC_S69_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./21-27_AHCY7CDSXC_S82_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./21-28_AHCY7CDSXC_S89_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./21-68_AHCY7CDSXC_S110_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./21-71_AHCY7CDSXC_S93_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./22-01_AHCY7CDSXC_S49_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./22-03_AHCY7CDSXC_S104_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./22-05_AHCY7CDSXC_S100_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./22-10_AHCY7CDSXC_S71_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./22-11_AHCY7CDSXC_S60_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./22-12_AHCY7CDSXC_S63_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./22-13_AHCY7CDSXC_S50_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./22-14_AHCY7CDSXC_S61_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./22-18_AHCY7CDSXC_S80_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./22-19_AHCY7CDSXC_S96_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./22-20_AHCY7CDSXC_S39_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./22-28_AHCY7CDSXC_S35_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./22-37_AHCY7CDSXC_S72_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./22-38_AHCY7CDSXC_S90_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./22-39_AHCY7CDSXC_S34_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./22-43_AHCY7CDSXC_S52_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./22-44_AHCY7CDSXC_S73_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./22-45_AHCY7CDSXC_S19_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./22-49_AHCY7CDSXC_S87_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./22-50_AHCY7CDSXC_S48_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./22-51_AHCY7CDSXC_S26_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./22-52_AHCY7CDSXC_S55_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./23-27_AHCY7CDSXC_S79_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./23-28_AHCY7CDSXC_S74_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./23-29_AHCY7CDSXC_S33_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./23-30_AHCY7CDSXC_S27_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./23-32_AHCY7CDSXC_S78_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./23-33_AHCY7CDSXC_S77_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./23-35_AHCY7CDSXC_S53_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./23-36_AHCY7CDSXC_S83_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./23-37_AHCY7CDSXC_S76_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./23-38_AHCY7CDSXC_S37_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./23-39_AHCY7CDSXC_S38_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./23-40_AHCY7CDSXC_S59_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./23-41_AHCY7CDSXC_S30_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./23-42_AHCY7CDSXC_S56_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./23-43_AHCY7CDSXC_S54_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./23-44_AHCY7CDSXC_S43_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./23-45_AHCY7CDSXC_S86_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./23-46_AHCY7CDSXC_S45_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./23-47_AHCY7CDSXC_S99_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./23-48_AHCY7CDSXC_S17_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./23-49_AHCY7CDSXC_S81_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./23-50_AHCY7CDSXC_S85_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./23-52_AHCY7CDSXC_S103_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./23-53_AHCY7CDSXC_S68_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./23-54_AHCY7CDSXC_S105_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./23-55_AHCY7CDSXC_S22_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./23-56_AHCY7CDSXC_S46_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./23-57_AHCY7CDSXC_S75_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./23-58_AHCY7CDSXC_S40_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./23-59_AHCY7CDSXC_S51_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./23-60_AHCY7CDSXC_S44_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./23-61_AHCY7CDSXC_S57_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./23-62_AHCY7CDSXC_S66_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./23-63_AHCY7CDSXC_S94_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam",
"./23-64_AHCY7CDSXC_S62_L001_R1_001_val_1_bismark_bt2_pe.sorted.bam")

myobj <- processBismarkAln(file.list,
                           sample.id=list("21-105", "21-111", "21-112", "21-113", "21-118", "21-142", "21-143", "21-146", "21-147","21-176", "21-178", "21-180", "21-182", "21-183", "21-184", "21-27", "21-28", "21-68", "21-71", "22-01", "22-03", "22-05", "22-10","22-11", "22-12", "22-13", "22-14", "22-18", "22-19", "22-20", "22-28", "22-37", "22-38","22-39", "22-43", "22-44", "22-45", "22-49", "22-50", "22-51", "22-52", "23-27", "23-28","23-29", "23-30", "23-32", "23-33", "23-35", "23-36", "23-37", "23-38", "23-39", "23-40","23-41", "23-42", "23-43", "23-44", "23-45", "23-46", "23-47", "23-48", "23-49", "23-50","23-52", "23-53", "23-54", "23-55", "23-56", "23-57", "23-58", "23-59", "23-60", "23-61","23-62", "23-63", "23-64"),
                           assembly="ASM164957v2",
                           treatment=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
                           read.context="CpG",
                           save.folder = getwd(),
                           mincov = 10
)

saveRDS(myobj, file="/home/jbelik/RRBS/Clock/Data/Bismark_results/Clock-new-young.RDS")

