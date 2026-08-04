library(data.table)
library(stringr)
library(purrr)
library(dplyr)

# source("/work/bio-tanyj/scripts/TanYongjun_code.R", echo = F)
# options(scipen = 200) # not to use scienrific notation in bed files.
setwd("/scratch/2025-09-01/bio-shenw/Cardiovascular_disease_STARR-seq/workspace_tanyj_pipeline/workspace_nonrisk/mapping/UMI_dedup")

## rmdup and count (20210817/20210910)--------------
TrimmedLength <- 0 # length of sequence trimmed from 5'end of reads in `truncat_reads.sh`："-b 25"。
NrowsForTest <- c(1000, Inf)[2] # set Inf to use all data.

for(i in list.files(path = ".",
                    pattern = ".+_UMI.bed$")){
    # i <- "si562_C25_S1_UMI.bed" # test
    id <- stringr::str_replace_all(i, "(.+)_UMI.+", "\\1")

    SummaryDf <- data.frame(ID = id,
                            Num_All_Raw_Frag = NA, Num_dT_Raw_Frag = NA, Num_nodT_Raw_Frag = NA,
                            Num_Old_All_Uniq_Frag = NA, Num_Old_nodT_Uniq_Frag = NA, Num_Old_dT_Uniq_Frag = NA,
                            Num_All_UMI_uniq_Frag = NA, Num_nodT_UMI_Uniq_Frag = NA, Num_dT_UMI_Uniq_Frag = NA,
                            Total_Uniq_UMI_Used = NA, Total_N_in_all_UMI = NA, Rate_of_N_in_all_UMI = NA,
                            Median_Len_Frag_dT = NA, Median_Len_Frag_nodT = NA)
    if(!file.exists(str_c(i, ".rmduping"))){
        write.table(NA, file = str_c(i, ".rmduping"))
        cat(paste("\n>>>>>> sample ID: ", id))

        # load files
            # fragments
            cat(paste("\n       reading BED files which contain all fragments....\n     "))
            timestamp()
            df <- fread(i, nrows = NrowsForTest, header = F)
            cat(paste("\n       add ", TrimmedLength, "bp (trimmed in truncate step) in each side of fragments"))
            df$V2 <- df$V2 - TrimmedLength
            df$V2[df$V2 < 0] <- 1
            df$V3 <- df$V3 + TrimmedLength
            SummaryDf$Num_All_Raw_Frag <- nrow(df)

            # readID
            # cat(paste("\n       reading ID of reads which contain 18dT....\n        "))
            # timestamp()
            # reads18dT <- fread(file = str_c("Trimed_", id, "_18dT_reads.readsIDtab"),
                               # nrows = NrowsForTest, header = F)

        # mark fragments with dT
            cat(paste("\n       Mark fragments with dT....\n        "))
            timestamp()
            df$dT <- "N"
            # df$dT[df$V5 %in% reads18dT$V1] <- "Y"
            df$V5 <- NULL
            cat(paste("\n       Count raw fragments with or without dT....\n        "))
            timestamp()
            SummaryDf$Num_dT_Raw_Frag <- sum(df$dT == "Y")
            SummaryDf$Num_nodT_Raw_Frag <- sum(df$dT == "N")
            # rm(reads18dT)

        # UMI stats
            cat(paste("\n       Summary of UMI....\n        "))
            timestamp()
            SummaryDf$Total_Uniq_UMI_Used <- length(unique(df$V6))
            SummaryDf$Total_N_in_all_UMI <- sum(str_count(df$V6, "N"))
            SummaryDf$Rate_of_N_in_all_UMI <- SummaryDf$Total_N_in_all_UMI / (SummaryDf$Num_All_Raw_Frag * 8)
            UMIsummary <- df$V6 %>%
                table() %>%
                table() %>%
                as.data.frame()
            colnames(UMIsummary) <- c("Times_of_individual_UMI_used_in_raw_fragments", "Number_of_UMI")
            UMIsummary %>%
                dplyr::arrange(desc(Number_of_UMI)) %>%
                fwrite(file = str_c("Summary_of_used_UMI_", id, ".tsv", sep = ""),
                      quote = F, sep = "\t")

        # rmdup
            dfLs <- split(df, list(df$V1, df$V4, df$dT))
            rm(df)

            # old strategy (dT/nodT, strand-/+ in one file)
            cat(paste("\n       rmdup with old strategy (chr, start, end, strand)....\n     "))
            timestamp()
            UniqOldDf <- map_dfr(dfLs, function(x){
                x %>%
                    dplyr::select(V1, V2, V3, V4, dT) %>%
                    dplyr::distinct() %>%
                    return()
            }) %>%
                dplyr::select(1:5) %>%
                dplyr::arrange(V1, V2, V3)
            SummaryDf$Num_Old_All_Uniq_Frag <- nrow(UniqOldDf)
            SummaryDf$Num_Old_nodT_Uniq_Frag <- sum(UniqOldDf$dT == "N")
            SummaryDf$Num_Old_dT_Uniq_Frag <- sum(UniqOldDf$dT == "Y")
            UniqOldDf %>%
                dplyr::select(1:4) %>%
                fwrite(file = str_c("Uniq_Old_strategy_", id, ".bed", sep = ""),
                        quote = F, sep = "\t", col.names = F)
            rm(UniqOldDf)

        # UMI method
            cat(paste("\n       rmdup with UMI (chr, start, end, strand, UMI)....\n     "))
            timestamp()
            UniqUMIDf <- map_dfr(dfLs, function(x){
                x %>%
                    dplyr::distinct()
            })
            SummaryDf$Num_All_UMI_uniq_Frag <- nrow(UniqUMIDf)
            UniqUMIDfdT <- UniqUMIDf %>% 
                dplyr::filter(dT == "Y") %>% 
                dplyr::select(1:4) %>% 
                dplyr::arrange(V1, V2, V3)
            UniqUMIDfNodT <- UniqUMIDf %>%
                dplyr::filter(dT == "N") %>%
                dplyr::select(1:4) %>%
                dplyr::arrange(V1, V2, V3)
            SummaryDf$Num_nodT_UMI_Uniq_Frag <- nrow(UniqUMIDfNodT)
            SummaryDf$Num_dT_UMI_Uniq_Frag <- nrow(UniqUMIDfdT)
            UniqUMIDfdT %>%
                fwrite(file = str_c("Uniq_UMI_dT_notExtend_", id, ".bed", sep = ""),
                                 quote = F, sep = "\t", col.names = F)
            UniqUMIDfNodT %>%
                dplyr::arrange(V1, V2, V3) %>%
                fwrite(file = str_c("Uniq_UMI_NodT_", id, ".bed", sep = ""),
                                 quote = F, sep = "\t", col.names = F)

            # length distribution of fragments (dT before extention)
            cat(paste("\n       Summary length of unique fragments (UMI)....\n      "))
            timestamp()
            FragLenDist <- function(x){
                x %>%
                    dplyr::mutate(FragmentLength = V3 - V2) %>%
                    group_by(FragmentLength) %>%
                    dplyr::summarise(Num = n()) %>%
                    return()
            }
            UniqUMIDfdT %>%
                FragLenDist %>%
                fwrite(str_c("Length_distribution_UniqUMI_dT_", id, ".tab", sep = ""),quote = F, sep = "\t")
            UniqUMIDfNodT %>%
                FragLenDist %>%
                fwrite(str_c("Length_distribution_UniqUMI_nodT_", id, ".tab", sep = ""),quote = F, sep = "\t")
            UniqUMIDf %>%
                FragLenDist %>%
                fwrite(str_c("Length_distribution_UniqUMI_all_", id, ".tab", sep = ""),quote = F, sep = "\t")
            rm(UniqUMIDf)

        # extend fragments with dT for some part of sequences were truncated.(To median length of fragments without dT)
            MedianFragLen <- median(UniqUMIDfNodT$V3 - UniqUMIDfNodT$V2)
            cat(str_c("\n     The median length of fragments with or without dT: ",
                      MedianFragLen, median(UniqUMIDfdT$V3 - UniqUMIDfdT$V2), "\n       ", 
                      sep = " "))
            UniqUMIDfdT %>%
                dplyr::mutate(Len = V3 - V2,
                              V3 = ifelse(V4 == "+" & Len < MedianFragLen,
                                          V2 + MedianFragLen,
                                          V3),
                              V2 = ifelse(V4 == "-" & Len < MedianFragLen,
                                          V3 - MedianFragLen,
                                          V2)) %>%
                dplyr::mutate(V2 = ifelse(V2 < 1, 1, V2)) %>%
                dplyr::select(1:4) %>%
                fwrite(file = str_c("Uniq_UMI_dT_Extend_", id, ".bed", sep = ""),
                                 quote = F, sep = "\t", col.names = F)

        # Summary of fragments
            SummaryDf$Median_Len_Frag_dT <- median(UniqUMIDfdT$V3 - UniqUMIDfdT$V2)
            SummaryDf$Median_Len_Frag_nodT <- MedianFragLen
            SummaryDf %>%
                fwrite(file = str_c("Summary_of_fragments_", id, ".tab", sep = ""),
                      sep = "\t", quote = F)
    }
}

# ## rmdup with UMI (202105)-----------------
# TrimmedLength <- 5 # length of sequence trimmed from 5'end of reads in `truncat_reads.sh`
# for(i in list.files(path = "./",
#                     pattern = ".+UMI.bed$")){
#     # i <- "Test_si116_P1_UMI.bed"
#     if(!file.exists(str_replace_all(i, ".bed", ".rmduping"))){
#         write.table(NA, file = str_replace_all(i, ".bed", ".rmduping"))

#         # load
#             cat(paste("\n>>>>>> Start reading", i))
#             timestamp()
#             df <- fread(i)
#             cat(paste("\n       extend ", TrimmedLength, "bp in each side of fragments"))
#             df$V2 <- df$V2 - TrimmedLength
#             df$V3 <- df$V3 + TrimmedLength
#             NumTotalFrag <- nrow(df)
#             df <- split(df, df$V1)

#         # rmdup based on chr,start, end, strand. (the previous method)
#             cat("\n       load complete, start rmdup with old strategy....    ")
#             timestamp()
#             cat("         ")
#             NumUniqFragOldStrategy <- 0
#             dfOld <- list()
#             for(j in names(df)){
#                 cat(paste(j, " "))
#                 dfOld[[j]] <- df[[j]] %>%
#                     dplyr::select(V1, V2, V3, V4) %>%
#                     unique()
#                 NumUniqFragOldStrategy <- NumUniqFragOldStrategy + nrow(dfOld[[j]])
#             }

#             # write
#             cat("\n       write to bed file using old strategy to bed file    ")
#             timestamp()
#             dfOld %>%
#                 purrr::reduce(., rbind) %>%
#                 arrange(V1, V2, V3, V4) %>%
#                 fwrite(file = str_replace_all(i, "UMI.bed", "OldMethod_rmdup.bed"), 
#                        sep = "\t", row.names = F, quote = F, col.names = F)
#             rm(dfOld)

#         # split based on strand
#             cat("\n       splite into two files....    ")
#             timestamp()
#             dfF <- list()
#             dfR <- list()
#             for(j in names(df)){
#                 dfF[[j]] <- df[[j]] %>%
#                     dplyr::filter(V4 == "+")
#                 dfR[[j]] <- df[[j]] %>%
#                     dplyr::filter(V4 == "-")
#             }
#             NumRawFor <- map_int(dfF, nrow) %>% sum()
#             NumRawRev <- map_int(dfR, nrow) %>% sum()
#             rm(df)

#         # rmdup based on chr,start (based on direction),strand, and UMI.
#             cat("\n       rmdup based on chr, start, strand, and UMI    ")
#             timestamp()
#             cat("         ")
#             for(j in names(dfF)){
#                 cat(paste(j, " "))
#                 dfF[[j]] <- dfF[[j]] %>%
#                     dplyr::group_by(V1, V2, V4, V5) %>%
#                     dplyr::summarise(V3 = max(V3)) %>%
#                     dplyr::select(-V5) %>%
#                     dplyr::ungroup()
#                 dfR[[j]] <- dfR[[j]] %>%
#                     dplyr::group_by(V1, V3, V4, V5) %>%
#                     dplyr::summarise(V2 = min(V2)) %>%
#                     dplyr::select(-V5) %>%
#                     dplyr::ungroup()
#             }
#             NumForFragUni <- map_int(dfF, nrow) %>% sum()
#             NumRevFragUni <- map_int(dfR, nrow) %>% sum()

#             # merge
#             cat("\n       merging    ")
#             timestamp()
#             dfF <- purrr::reduce(dfF, rbind)
#             dfR <- purrr::reduce(dfR, rbind)
#             df <- rbind(dfF, dfR)
#             rm(dfF, dfR)

#             # sort
#             df <- dplyr::arrange(df, V1, V2, V3, V4)
#             NumToalFragUni <- nrow(df)

#             # write to file
#             cat("\n       write to file    ")
#             timestamp()
#             data.table::fwrite(df[c("V1", "V2", "V3", "V4")],
#                             file = str_replace_all(i, ".bed", "_rmdup.bed"),
#                             sep = "\t", row.names = F, col.names = F)

#         # Report Num of fragments
#             write.table(data.frame(NumRawFragments = NumTotalFrag,
#                                 NumUniqFragOldStrategy = NumUniqFragOldStrategy,
#                                 NumRawFor = NumRawFor,
#                                 NumRawRev = NumRawRev,
#                                 NumUMIUniqFragments = NumToalFragUni,
#                                 NumUMIUniqFor = NumForFragUni,
#                                 NumUMIUniqRev = NumRevFragUni),
#                         file = str_replace_all(i, ".bed", "_NumOfFragments.tsv"),
#                         sep = "\t", row.names = F, quote = F, col.names = F)
#     }
# }

# ## Summary of UMI (20210712)---------------
# #Count: Total number of UMI used, Number of PCR dup for the same fragments.
# for(i in list.files(path = "./",
#                     pattern = ".+UMI.bed$")){
#     # i <- "Test_si116_P1_UMI.bed"
#     if(!file.exists(str_replace_all(i, ".bed", ".Counting"))){
#         write.table(NA, file = str_replace_all(i, ".bed", ".Counting"))

#         # load
#         cat(paste("\n>>>>>> Start reading", i))
#         timestamp()
#         df <- fread(i)
#         NumTotalFrag <- nrow(df)

#         # count
#         NumUMI <- table(df$V5)
#         TotalNumberUMI <- length(NumUMI)
#         QuantileAllUMI <- table(NumUMI) %>% as.data.frame()

#         DupFragmentsNoUMI <- df %>%
#            dplyr::group_by(V1, V2, V3, V4) %>%
#            dplyr::summarise(Num = n()) %>%
#            dplyr::ungroup() %>%
#            dplyr::select(Num) %>%
#            unlist() %>%
#            table() %>%
#            as.data.frame()

#         DupFragmentsUMI <- df %>%
#            dplyr::group_by(V1, V2, V3, V4, V5) %>%
#            dplyr::summarise(Num = n()) %>%
#            dplyr::ungroup() %>%
#            dplyr::select(Num) %>%
#            unlist() %>%
#            table() %>%
#            as.data.frame()

#         # write to file
#         fileName <- str_replace_all(i, "_UMI.bed", ".tsv")
#         write.table(str_c("NumTotalFragments:", NumTotalFrag), file = fileName)
#         write.table("Table of All UMI's in all fragments (not take fragments into consideration)",
#                     file = fileName, append = T, sep = "\t")
#         fwrite(QuantileAllUMI, file = fileName, append = T, sep = "\t", row.names = T, col.names = T)
#         write.table("Table of dups of the same fragments (not use UMI info)",
#                     file = fileName, append = T, sep = "\t")
#         fwrite(DupFragmentsNoUMI, file = fileName, append = T, sep = "\t", row.names = T, col.names = T)
#         write.table("Table of dups of the same fragments+UMI",
#                     file = fileName, append = T, sep = "\t")
#         fwrite(DupFragmentsUMI, file = fileName, append = T, sep = "\t", row.names = T, col.names = T)
#     }
# }




