# Analyzing 16S Microbiome Insight data - added phyloseq

source("package.R")

dada2analysis <- drake_plan(
  import_metadata = read.csv("data/meta.csv"),
  sample_name_list = as.character(levels(import_metadata$MBI.Sample.ID)),
  metadata = import_metadata %>% 
    set_rownames(sample_name_list) %>% 
    select(-c(MBI.Sample.ID)),
  set_path = here("raw_sequences/"),
  files = list.files(set_path, pattern = "fastq.gz"),
  fnFs = sort(list.files(set_path, pattern = "R1_001.fastq.gz", full.names = TRUE)),
  fnRs = sort(list.files(set_path, pattern = "R2_001.fastq.gz", full.names = TRUE)),
  sample.names = sapply(strsplit(basename(fnFs), "_"), `[`, 1),
  filtFs = file.path(set_path, "filtered", paste0(sample.names, "_F_filt.fastq.gz")),
  filtRs = file.path(set_path, "filtered", paste0(sample.names, "_R_filt.fastq.gz")),
  out = filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                      maxN=0, maxEE=c(2,30), truncQ=2, rm.phix=TRUE,
                      compress=TRUE, multithread=FALSE),
  errF = learnErrors(filtFs, multithread=TRUE),
  errR = learnErrors(filtRs, multithread=TRUE),
  derepFs = derepFastq(filtFs, verbose=TRUE) %>% 
    setNames(sample.names),
  derepRs = derepFastq(filtRs, verbose=TRUE) %>% 
    setNames(sample.names),
  dadaFs = dada(derepFs, err=errF, multithread=TRUE),
  dadaRs = dada(derepRs, err=errR, multithread=TRUE),
  mergers = mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE),
  seqtab = makeSequenceTable(mergers),
  seqtab.nochim = removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE),
  getN = function(x) sum(getUniques(x)),
  track = cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim)) %>% 
    set_colnames(c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")) %>% 
    set_rownames(sample.names),
  taxa = assignTaxonomy(seqtab.nochim, here("taxa/silva_nr_v138_train_set.fa.gz"), multithread=TRUE),
  taxa_sp = addSpecies(taxa, here("taxa/silva_species_assignment_v138.fa.gz")),
  taxa.print = set_rownames(taxa_sp, NULL),
  ps = phyloseq(otu_table(as.matrix(seqtab.nochim), taxa_are_rows = FALSE),
                sample_data(metadata),
                tax_table(as.matrix(taxa_sp))),
  dna = Biostrings::DNAStringSet(taxa_names(ps)),
  dna_name = set_names(dna, taxa_names(ps)),
  ps_object = merge_phyloseq(ps, dna_name),
  add_asv = paste0("ASV", seq(ntaxa(ps))),
  taxa_name_ps = taxa_names(ps_object)
)

make(dada2analysis)

loadd(set_path, files, fnFs, fnRs, sample.names, filtFs, filtRs, out,
      errF, errR, derepFs, derepRs, dadaFs, dadaRs, mergers, seqtab,
      seqtab.nochim, getN, track, taxa, taxa_sp, taxa.print, import_metadata,
      sample_name_list, metadata, ps_object, add_asv, taxa_name_ps,
      ps, dna, dna_name, ps_object)

vis_drake_graph(dada2analysis)
