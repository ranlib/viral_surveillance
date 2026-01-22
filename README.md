# viral_surveillance
Viral surveillance analysis geared towards the Illumina viral surveillance panel.

# Overview of Pipeline

* Read QC & trimming

* Human/host depletion (alignment to human)

* Taxonomic classification (virus detection)

* Reference alignment for detected viruses

* Variant calling

* Consensus genome generation

* Phylogenetic placement

Each step is a task with a Docker container specified.

# Notes

* This workflow assumes pre-computed human and viral reference indices
  (BWA/Bowtie2) are available.

* The Viral Surveillance Panel targets >200 viruses including RNA and
  DNA viruses, useful for surveillance and evolutionary studies.

* You can expand taxonomic detection to include abundance estimates
  (e.g., Bracken for Kraken2 output) or add dragen or similar tools if
  you have access.

* Depth and allele frequency thresholds (min_depth / min_var_af) are
  configurable inputs for variant calling.
  