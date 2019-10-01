## Course schedule

![](diagram.png)

The figure was generated by [workflow.R](workflow.R).


## Data

+ [test.db](data/test.db)
+ parsing genbank
    - [AB115403.gb](data/AB115403.gb)    
    - [AB115403.fasta](data/AB115403.fasta)
    - [some gb files](data/gb)
+ [GEO](data/GEO)
+ Flu
    - [flu_seq.fas](data/flu_seq.fas)
    - [flu_seq_v2.fas](data/flu_seq.fas)
+ Beta Coronavirus
    - [accession number](data/betaCov/betaCov_acc.txt)
    - [genbank files](data/betaCov/gb)
    - [fasta sequence](data/betaCov/betaCov.fas)
    - [partial rdrp fasta sequence](data/betaCov/rdrp.fas)
    - [aligned partial rdrp fasta sequence](data/betaCov/rdrp_aln.fas)
    - [rdrp newick tree](data/betaCov/rdrp_aln.fas.treefile)


## R package

+ [bioinfo.practice](bioinfo.practice)
    - demo functions you need to implement in the class
+ [buffalo](buffalo)
    - demo R package development


## Links and resources

+ SQLite
    - <https://www.sqlite.org/aff_short.html>
    - [How to Use SQLite with R](https://www.bioconductor.org/help/course-materials/2006/rforbioinformatics/labs/thurs/SQLite-R-howto.pdf) 
    - [RSQLite Tutorial](https://github.com/ysquared2/RSQLiteTutorial)
    - [SQLite in R](https://www.datacamp.com/community/tutorials/sqlite-in-r)
    - [A Gentle Introduction to SQL Using SQLite](https://a-gentle-introduction-to-sql.readthedocs.io/en/latest/)
+ GenBank 
    - <https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html>
    - <https://www.ncbi.nlm.nih.gov/nuccore/AB115403>
    - <https://www.ncbi.nlm.nih.gov/books/NBK25501/>
+ Alignment
    - <https://www.nature.com/articles/nbt0704-909>
        - [global.c](global.c)
    - <https://developer.ibm.com/articles/j-seqalign/>
    - <https://www.geeksforgeeks.org/sequence-alignment-problem/>
    - <https://mp.weixin.qq.com/s/4DL-pJVItOhkYJAbPEg7BQ>
+ Make
    - [Make命令教程](http://www.ruanyifeng.com/blog/2015/02/make.html)
    - [A simple makefile tutorial](http://www.cs.colby.edu/maxwell/courses/tutorials/maketutor/)
+ git/github
    - <https://help.github.com/en/articles/set-up-git>
    - <http://jaimeiniesta.github.io/learn.github.com/p/diff.html>
    - <https://git-scm.com/book/en/v2/Git-Branching-Basic-Branching-and-Merging>
+ Reproducible research
    - <http://www.nature.com/news/reproducibility-1.17552>    
    - [Repeatability of published microarray gene expression analyses](https://www.ncbi.nlm.nih.gov/pubmed/19174838)
    - <https://github.com/GuangchuangYu/plotting_tree_with_data>
    - <https://bookdown.org/yihui/rmarkdown/>
+ Differential gene analysis
    - Marioni, J. C., Mason, C. E., Mane, S. M., Stephens, M. & Gilad, Y. RNA-seq: an assessment of technical reproducibility and comparison with gene expression arrays. Genome Res. 18, 1509-1517 (2008).
    - Robinson, M. D. & Oshlack, A. A scaling normalization method for differential expression analysis of RNA-seq data. Genome Biol. 11, R25 (2010).
    - Anders, S. et al. Count-based differential expression analysis of RNA sequencing data using R and Bioconductor. Nat Protoc 8, 1765-1786 (2013).
    - Lund, S. P., Nettleton, D., McCarthy, D. J. & Smyth, G. K. Detecting differential expression in RNA-sequence data using quasilikelihood with shrunken dispersion estimates. Stat Appl Genet Mol Biol 11, (2012).
    - Anders, S. & Huber, W. Differential expression analysis for sequence count data. Genome Biol. 11, R106 (2010).
    - Law, C. W., Chen, Y., Shi, W. & Smyth, G. K. Voom: precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biol. 15, R29 (2014).
+ GSEA
    - Subramanian, A. et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. Proc. Natl. Acad. Sci. U.S.A. 102, 15545-15550 (2005).
    - <https://yulab-smu.github.io/clusterProfiler-book>
    - <https://bioconductor.org/packages/release/BiocViews.html#___GeneSetEnrichment>
+ Multiple Sequence Alignment
    - <https://www.bioconductor.org/packages/muscle/>
    - <https://www.bioconductor.org/packages/msa>
    - <https://cran.r-project.org/package=ape>
    - <http://www.mbio.ncsu.edu/BioEdit/bioedit.html>
    - <https://mp.weixin.qq.com/s/F0TxMk1CrYVNvxWqZgv-Aw>
+ Phylogeny
    - <https://yulab-smu.github.io/treedata-book/>
    - [Evolution 101](https://evolution.berkeley.edu/evolibrary/article/evo_01)
    - [How to read a phylogenetic tree](https://artic.network/how-to-read-a-tree.html)
    - [Reading a Phylogenetic Tree: The Meaning of Monophyletic Groups](https://www.nature.com/scitable/topicpage/reading-a-phylogenetic-tree-the-meaning-of-41956/)