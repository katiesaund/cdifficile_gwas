# Goal: svae cdu1-tcdR intergenic region as a fasta file

library(ape)

cdu1_tcdR_intergenic_start_position <- 785614 # I just looked at the CD630 genome https://www.ncbi.nlm.nih.gov/nuccore/AM180355.1?&withparts=on&expand-gaps=on
cdu1_tcdR_intergenic_stop_position <- 786520


cd630_fa <- read.FASTA(file = "../../data/16_reference_genome/cdiff_630.fasta", type = "DNA")
cdu1_tcdR_fasta <- cd630_fa
cdu1_tcdR_fasta$`gi|126697566|ref|NC_009089.1| Peptoclostridium difficile 630, complete genome, plasmid` <- cdu1_tcdR_fasta$`gi|126697566|ref|NC_009089.1| Peptoclostridium difficile 630, complete genome, plasmid`[cdu1_tcdR_intergenic_start_position:cdu1_tcdR_intergenic_stop_position]
write.FASTA(cdu1_tcdR_fasta, file ="../../data/10_paloc/cdu1_tcdR.fna", header = "cdu1-tcdR")
# rename header on command line now

# foo <- seqinr::read.fasta(file = "../../../../Sequence_data/Project_Hanna_collections/reference_genome/cdiff_630.fasta", seqtype = "DNA")
# foo$`gi|126697566|ref|NC_009089.1|`[cdu1_tcdR_intergenic_start_position:cdu1_tcdR_intergenic_stop_position]
# foo$`gi|126697566|ref|NC_009089.1|`[(cdu1_tcdR_intergenic_stop_position + 1):(cdu1_tcdR_intergenic_stop_position + 100)]

cdu1_start <- 785233
hyp_protein_stop <- 805679
paloc_plus_fasta <- cd630_fa
paloc_plus_fasta$`gi|126697566|ref|NC_009089.1| Peptoclostridium difficile 630, complete genome, plasmid` <- 
  paloc_plus_fasta$`gi|126697566|ref|NC_009089.1| Peptoclostridium difficile 630, complete genome, plasmid`[cdu1_start:hyp_protein_stop]
write.FASTA(paloc_plus_fasta, file ="../../data/10_paloc/paloc_plus.fna", header = "cdu1_thru_hypo_paloc_plus")
