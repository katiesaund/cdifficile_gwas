module load Bioinformatics
module load ncbi-blast/2.9.0
/nfs/esnitkin/bin_group/blast_suite/bin/blastx -db /scratch/esnitkin_root/esnitkin/apirani/Testing_pipelines/Variant_calling/2020_08_06_Test_Provean/nr_database/v4_db/nr -query ../../data/17_blast_overlap_genes/toxin_severity_sig_overlap_indiv_hits.fna  -evalue 1e-4 -outfmt 10 -max_target_seqs 1 -max_hsps 1 -out ../../data/17_blast_overlap_genes/blast_toxin_sev_overlap_accessory_hits_vs_nr_db.csv -query_gencode 11

