# Somersaulting Pigeon

This repository contains the analysis scripts associated with the study of the
genetic, molecular, cellular, and behavioral basis of somersaulting behavior in
domestic pigeons.

## Repository structure

### `01_bulk_RNAseq`

Scripts for bulk RNA-seq analyses.

- `01.DEGs`: Differential gene-expression analysis
- `02.WGCNA`: Weighted gene co-expression network analysis

### `02_population_genomics`

Scripts for population genomic analyses.

- `01.SNPcalling`: SNP calling and variant processing
- `02.phylogeneticTree`: Phylogenetic tree construction
- `03.ADMIXTURE`: Population ancestry and structure analysis
- `04.PCA`: Principal component analysis
- `05.LSBL`: Locus-specific branch length analysis
- `06.Pi_ratio`: Nucleotide diversity ratio analysis
- `07.XP-CLR`: Cross-population composite likelihood ratio analysis
- `08.GO_semantic_similarity`: Gene Ontology semantic similarity analysis

### `03_single_nucleus_RNAseq`

Scripts for single-nucleus RNA-seq analyses.

- `01.expression_matrix`: Expression-matrix generation and preprocessing
- `02.scrublet`: Doublet detection using Scrublet
- `03.quality_control`: Quality control and cell filtering
- `04.annotation`: Data integration, clustering, and cell-type annotation
- `05.scCODA`: Differential cell-composition analysis
- `06.UCell`: Single-cell gene-signature scoring using UCell
- `07.GSEA`: Gene set enrichment analysis
- `08.scTenifoldKnk`: Cell-type-specific virtual-knockout analysis
- `09.cellChat`: Cell–cell communication analysis using CellChat
- `10.DiseaseGeneSetEnrichment`: Disease-associated gene-set enrichment analysis

### `04_zebrafish_tracking`

Python-based zebrafish behavioral tracking and analysis pipeline.

- `configs`: Configuration files
- `data`: Input data and example data structure
- `docs`: Documentation
- `scripts`: Executable analysis scripts
- `src`: Source code for tracking, trajectory processing, and behavioral analysis
- `tests`: Test scripts
- `checkpoints`: Model or analysis checkpoints
- `outputs`: Generated analysis outputs
- `legacy`: Earlier or deprecated versions of scripts
- `pyproject.toml`: Python project configuration
- `requirements.txt`: Python package dependencies
- `README.md`: Detailed usage instructions for the zebrafish tracking pipeline

## Usage

Scripts within each analysis directory are organized according to the
corresponding workflow.

Before running the analyses, users should:

1. Update input and output paths.
2. Check software and package dependencies.
3. Review user-defined parameters in each script.
4. Follow the execution order indicated by the directory and script numbering.

For the zebrafish tracking pipeline, please refer to:

```text
04_zebrafish_tracking/README.md

## Citation

Citation information will be added after publication.

## License

License information will be added.
