# RNA velocity for NB patient bone marrow samples

## Prepare files

Download 10x Human reference (GRCh38) dataset (version 2020-A from 2020-07-07).

```bash
curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
tar xfvz refdata-gex-GRCh38-2020-A.tar.gz
```

Download expressed repeat annotation from [here](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=611454127_NtvlaW6xBSIRYJEBI0iRDEWisITa&clade=mammal&org=Human&db=0&hgta_group=allTracks&hgta_track=rmsk&hgta_table=rmsk&hgta_regionType=genome&position=&hgta_outputType=gff&hgta_outFileName=hg38_rmsk.gtf).

```bash
gunzip hg38_rmsk.gtf.gz
```

Prepare directories.

```bash
mkdir metadata
mv refdata-gex-GRCh38-2020-A/genes/genes.gtf metadata/
mv hg38_rmsk.gtf metadata/
rm -rf refdata-gex*
```

If not available, assemble a singularity container with `velocyto`, `scvelo`, and their dependencies. Use the following definition file: [csbg/containers/scvelo.def](https://github.com/csbg/containers/blob/main/scvelo.def).


## Generate spliced/unspliced counts

For each sample, the following files are required:
- `possorted_genome_bam.bam`
- `possorted_genome_bam.bam.bai`
- `filtered_feature_bc_matrix/barcodes.tsv.gz`

In order to process these files for a single sample, execute `run_velocyto.sh` with the sample name as first argument. This will store spliced/unspliced counts in a loom file `output_[sample name]/possorted_genome_bam.loom`.

Several samples can be processed, e.g., via an SGE array job. To this end, prepare a file `samples.txt` with each sample name on a line, and start the array job via

```bash
qsub sge_run_velocyto.sh
```


## Calculate RNA velocity

Run `scvelo` via `run_scvelo.sh`, which executes `analyse_velocity.py` in the singularity container.


## Visualization

Results of the RNA velocity analysis are plotted in the Jupyter notebook `plot_velocity.ipynb`.


## Links

- https://github.com/basilkhuder/Seurat-to-RNA-Velocity
- http://velocyto.org/velocyto.py/tutorial/cli.html
