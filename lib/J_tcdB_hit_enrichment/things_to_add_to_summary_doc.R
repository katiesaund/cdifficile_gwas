sum(hogwash_continuous$hit_pvals$fdr_corrected_pvals > -log(0.15)) # 220

sum(hogwash_continuous$hit_pvals$fdr_corrected_pvals > -log(0.15) & 
      hogwash_continuous$convergence$epsilon > (0.15)) # 8


sum(grepl("tcdB", row.names(hogwash_continuous$hit_pvals)[hogwash_continuous$hit_pvals$fdr_corrected_pvals > -log(0.15)]))
sum(grepl("tcdE", row.names(hogwash_continuous$hit_pvals)[hogwash_continuous$hit_pvals$fdr_corrected_pvals > -log(0.15)]))
sum(grepl("tcdC", row.names(hogwash_continuous$hit_pvals)[hogwash_continuous$hit_pvals$fdr_corrected_pvals > -log(0.15)]))
sum(grepl("tcdR", row.names(hogwash_continuous$hit_pvals)[hogwash_continuous$hit_pvals$fdr_corrected_pvals > -log(0.15)]))
sum(grepl("tcdA", row.names(hogwash_continuous$hit_pvals)[hogwash_continuous$hit_pvals$fdr_corrected_pvals > -log(0.15)]))


load("../../data")

sum(grepl("[|]tcdR-tcdB[|]", row.names(hogwash_continuous$hit_pvals)[hogwash_continuous$hit_pvals$fdr_corrected_pvals > -log(0.15)])) #2

sum(grepl("[|]tcdB[|]", row.names(hogwash_continuous$hit_pvals)[hogwash_continuous$hit_pvals$fdr_corrected_pvals > -log(0.15)])) # 75

tcdB_names <- row.names(hogwash_continuous$hit_pvals)[hogwash_continuous$hit_pvals$fdr_corrected_pvals > -log(0.15)][grep("[|]tcdB[|]", row.names(hogwash_continuous$hit_pvals)[hogwash_continuous$hit_pvals$fdr_corrected_pvals > -log(0.15)])]

tcdB_names[grepl("1967", tcdB_names)]
tcdB_names[grepl("1231", tcdB_names)]


sum(grepl("[|]tcdR-tcdB[|]", row.names(hogwash_continuous$hit_pvals))) #2


# grouped
# gene="tcdB"
# /locus_tag="CD630_06600"
load("../../data/5_hogwash/log_toxin/hogwash_continuous_grouped_continuous_log_toxin_grouped_post-ar.rda")
hogwash_continuous$hit_pvals[grepl("CD630_06600", row.names(hogwash_continuous$hit_pvals)), , drop = FALSE]
