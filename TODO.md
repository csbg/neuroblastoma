# TODO

* doublet removal (see chapter 16 in OSCA)
* densMAP (https://github.com/hhcho/densvis/tree/master/densmap)
* compare number of NB cells (esp. in group III) to tumour infiltration rates
  - also correlate fraction of HSCs and total hematopoietic fraction
* think about a way to include several layers of information (samples, groups,
  clusters, cell types ...) when visualizing gene expression
* use SingleR cell type classification to split clusters into subclusters automatically; e.g., 
  - for each cluster, remove low-abundant cell types (like <10%)
  - split remainder of cluster by cell type -> subclusters like `5_Neuron` or `5_pro-B`
  - this should yield clearer gene signatures
* DGE analysis on these subclusters



# References

OSCA: Orchestrating single-cell analysis with Bioconductor (http://bioconductor.org/books/release/OSCA/)