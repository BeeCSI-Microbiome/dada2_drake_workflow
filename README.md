# DADA2 workflow using drake

## About

The code in this repository is for a workflow for analyzing microbial community data using the `DADA2` R package. We have built a workflow with the `drake` R package, which enables the creation of reproducible workflows for automating analyses.

We were inspired by other projects such as `ampliseq` for automating QIIME2 runs, but in this case we are looking into providing an analysis tool that can be run entirely within R.

## Dependencies

Install R & Download R Studio. (Version used 3.6.3)
In order to run this workflow, the following packages need to be installed in R Studio
1. `tidyverse`	  (version 1.3.0)
2. `here`		      (version 0.1)
3. `DADA2`		    (version 1.14.1)
4. `drake`		    (version 7.12.0)
5. `magrittr`    	(version 1.5)

## Targets

| Description        | Function          |
| ------------------ | ------------------|
| Quality filtering  | `filterAndTrim()` |
| Dereplication      | `derepFastq()` |
| Learn error rates  | `learnErrors()` |
|Sample inference	   | `dada()` |
|Merge of paired reads	| `mergePairs()` |
|Construct a sequence table |	`makeSequenceTable()` |
|Chimera removal	| `removeBimeraDenovo()` |
|Taxonomic Classification	| `assignTaxonomy()`, `addSpecies()`|

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

## Visualize drake plan
> `vis_drake_graph(dada2analysis)`

![](/16s_rrna/vis_drake_graph.png)

## Reference
* DADA2 tutorial - https://benjjneb.github.io/dada2/tutorial.html
* The drake R Package User Manual -	https://books.ropensci.org/drake/
* ampliseq - https://github.com/nf-core/ampliseq
