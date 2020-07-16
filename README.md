# DADA2 workflow using drake

## About

The code in this repository is for a workflow for analyzing microbial community data using the `DADA2` and `phloseq` R package. We have built a workflow with the `drake` R package, which enables the creation of reproducible workflows for automating analyses.

We were inspired by other projects such as `ampliseq` for automating QIIME2 runs, but in this case, we are looking into providing an analysis tool that can be run entirely within R.

## Dependencies

Install R & Download R Studio. (Version used 3.6.3)
In order to run this workflow, the following packages need to be installed in R Studio

1. `drake`		    (version 7.12.0)
2. `DADA2`		    (version 1.14.1)
3. `phyloseq`     (version 1.30.0)
4. `tidyverse`	  (version 1.3.0)
5. `here`		      (version 0.1)
6. `magrittr`    	(version 1.5)

## Folder Structure

![](folder_structure_380.png)

* set the working directory `to 16s_rrna/`
* enter input data to the `16s_rrna/raw_sequences/` folder 

## Input
Before running samples through the workflow, ensure reads are:
* Demultiplexed
* Adapter and primer sequences are removed
* And forward and reverse reads are matched with file names (if pair-end sequencing) 

## Targets

| Description        | Function          |
| ------------------ | ------------------|
| Quality filtering  | `filterAndTrim()` |
| Dereplication      | `derepFastq()` |
| Learn error rates  | `learnErrors()` |
| Sample inference	 | `dada()` |
| Merge of paired reads	| `mergePairs()` |
| Construct a sequence table |	`makeSequenceTable()` |
| Chimera removal	| `removeBimeraDenovo()` |
| Taxonomic Classification	| `assignTaxonomy()`, `addSpecies()`|

*Documentation is reproduced from `DADA2` package*

Two taxonomic tables are created, with and without species classification and can be omitted from the output, but is added as an option. `assignSpecies()` can be used instead of `assignTaxonomy()` to get a species classification. 

## Analysis Functions

Inspect read quality profiles (forward and reverse reads)
> `plotQualityProfile(fnFs[<index_start>:<index_end>])`

*Output:*
* Gray scale heat map : frequency of each quality score at each base position
* Green line : mean quality score at each position is shown
* Orange lines : quartiles of the quality score distribution
* Red line : scaled proportion of reads that extend to at least that position

Plot observed and estimated error rate for forward and reverse reads
> `plotErrors(errF)`

## Applying modifications based on personalized data
|  Function       |   Description    |
|-------- |-------|
|`filterAndTrim(truncLen = c(<>,<>))` | Default 0 (no truncation). Truncate reads after truncLen bases. Reads shorter than this are discarded.|
|`removePrimers()` |	Removes primers and orients reads in a consistent direction. Intended use for PacBio CSS data|
|`seqComplexity()` | Determine if input sequence(s) are low complexity.|
|`show(<dada2 object>)`	| Method extensions to show for dada2 objects.|
|`tperr1`	| An empirical error matrix|

*Documentation is reproduced from drake package*

## Usage
| Function  | Description  |
|---|---|
|`make(<drake plan>)`	| Run the drake plan |
|`loadd(<drake target>)` |	Run and load drake target |
|`readd(<drake target>)` |	Run and return drake target |
|`vis_drake_graph(<drake plan>)` |	Show an interactive visual network representation of drake plan |
|`clean()` |	Remove Targets/Imports From The Cache |

## Taxonomic Reference Data 

DADA2 maintains [reference fastas](https://benjjneb.github.io/dada2/training.html) for the three most common 16S databases:  
* Silva  
* RDP 
* GreenGenes 
The workflow uses Silva (version 138) but can be modified integrate any of the 3 databases above. 

## Visualize drake plan
* quickly view the status of targets
* `drake` tracks changes in targets
* overview of target dependencies

> Example: Up to date `vis_drake_graph(dada2analysis)`

![](/16s_rrna/vis_drake_graph.png)


> Example: Outdated `vis_drake_graph(dada2analysis`

![](/16s_rrna/vis_drake_graph_outdated.png)

## Analysis of phyloseq Object 

| Function | Description |
|----------|-------------|
| `sample_variables(<phyloseq object>)` | Get The Sample Variables Present In Sample Data |
| `sample_data(<phyloseq object>)`      | Build Or Access Sample Data |
| `refseq(<phyloseq object>)`           | Retrieve Reference Sequences (XStringSet-Class) From Object |
| `taxa_sums(<phyloseq object>)`        | Returns The Total Number Of Individuals Observed From Each Species/Taxa/OTU |
| `ntaxa(<phyloseq object>)`            | Number of Taxa |
| `rank_names(<phyloseq object>)`       | Taxonomic Ranks |
| `otu_table(<phyloseq object>)`        | Access or build the OTU Table |
| `tax_table(<phyloseq object>)`        | Access or build the Taxonomy Table |

## Modifications to phyloseq Object 

| Function | Description |
|----------|-------------|
| `prune_samples(<phyloseq object>)` | Define A Subset Of Samples To Keep In A Phyloseq Object |
| `prune_taxa(<phyloseq object>)` | Prune Unwanted OTUs / Taxa From A Phylogenetic Object |
| `merge_samples(<phyloseq object>)` | Merge Samples Based On A Sample Variable Or Factor |

## Notes
This is intended for data from Illumina sequencing, but can be modified for other sequencing data formats such as 454 and Ion Torrent sequencing data with some modification in the DADA2 parameters. For more information, please refer to the [DADA2 documentation](https://www.bioconductor.org/packages/release/bioc/manuals/dada2/man/dada2.pdf). 

## Reference
* DADA2 tutorial - https://benjjneb.github.io/dada2/tutorial.html
* phyloseq tutorial - https://joey711.github.io/phyloseq/
* The drake R Package User Manual -	https://books.ropensci.org/drake/
* ampliseq - https://github.com/nf-core/ampliseq
