# micPDP
This repository contains a tool to identifiy micropeptides from genome wide multiple species alignments. The micpdp tool has been used in

[Mackowiak at al, Extensive identification and analysis of conserved small ORFs in animals. Genome Biology 2015](https://pubmed.ncbi.nlm.nih.gov/26364619)

[Bazzini et al, Identification of small ORFs in vertebrates using ribosome footprinting and evolutionary conservation. Embo J. 2014](https://pubmed.ncbi.nlm.nih.gov/24705786/)

## Authorship and legal information
The scripts, in the form provided in this archive, were implemented by Sebastian Mackowiak, being a PhD student at this time in the lab of Nikolaus Rajewsky at the [Max-Delbrueck-Center for Molecular Biology](https://www.mdc-berlin.de/1151037/en/research/research_teams/systems_biology_of_gene_regulatory_elements). The code is hereby released is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License](LICENSE).
You should have received a copy of the license along with this
work. If not, see [http://creativecommons.org/licenses/by-nc-sa/4.0/](https://creativecommons.org/licenses/by-nc-sa/4.0/)

## Installation

```
git clone https://github.com/Drmirdeep/micpdp.git
```

The tool requires mutliple species alignments for regions of interest to call micropeptides in them. For a whole genome wide analysis maf blocks were downloaded from [UCSC](https://hgdownload.soe.ucsc.edu/downloads.html) and stitched together with the scripts from the [stitch_maf package](https://github.com/Drmirdeep/stitch_maf).
Detailed instructions on how to build these alignments can be found in the file build_genome_maf_plus_index.txt.

Default scoring of micpdp is done by using Ka-Ks ratios. As an alternative scoring procedure the [PhyloCSF](https://pubmed.ncbi.nlm.nih.gov/21685081/) tool can be used.

## Example use case:
```bash
micpdp.pl -m chrI/chrI.maf.stitched.cmpl  -c chrI -l 27 -s species_ce10_reordered -t project_cel -b chrI_refseq.bed -T 0.5 -S  -F stitched_mafs
```
Would identify micropeptides in transcripts of chromosome I with a minimum length of 9 AA (27 nt => -l 27) and keeps sequences with less than 50% of gaps.
