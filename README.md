# busco2snap.pl v0.2

## Description
__Creation of a SNAP model from [BUSCO](https://gitlab.com/ezlab/busco) output.__

Using [MAKER's](https://www.yandell-lab.org/software/maker.html) `cegma2zff` as archetype, `busco2snap.pl` creates a SNAP model from complete, single-copy BUSCOs by predicting the gene structure with the retrained Augustus model and converting to zff format. The tools `augustus` and `blastp` need to be in your `$PATH`. If `fathom`, `forge` and `hmm-assembler.pl` are in your `$PATH` the SNAP model will be created automatically.   
If Augustus predicts more than one gene at a BUSCO locus, the predicted proteins are blasted against the ancestral variants of the corresponding BUSCO lineage. The gene(s) where the best hit (defined by highest bit score) is on the searched BUSCO are used only.

#### Tested with BUSCO versions 3.0.2, 4.1.4 and 5.2.2

## Dependencies

`busco2snap.pl` needs the following Perl modules and will search for executables in your `$PATH`:

Mandatory:
- [Parallel::Loops](https://metacpan.org/pod/Parallel::Loops)
- [Augustus](https://github.com/Gaius-Augustus/Augustus): `augustus`
- [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download): `blastp`

Optional:
- [SNAP](https://github.com/KorfLab/SNAP): `fathom`, `forge`, `hmm-assembler.pl`

## Usage

```
busco2snap.pl -b <busco_dir>

Mandatory:
	-b STR		BUSCO output directory (run_*)
			(augustus_output/retraining_parameters/)

Options: [default]
	-dup		Include duplicated BUSCOs [off]
	-frag		Include fragmented BUSCOs [off]
	-o STR		Output directory [.]
	-pre STR	Prefix for file name of the SNAP model file [species model name]
	-t INT		Number of parallel processed BUSCOs [1]
	-kt		Keep temporary files [off]
	-v		Print executed commands to STDERR [off]
	-dry-run	Only print commands to STDERR instead of executing [off]

	-h or -help	Print this help and exit
	-version	Print version number and exit
```

## Citation

__If you use this tool please cite the dependencies:__

- Augustus:  
Stanke M, Steinkamp R, Waack S, Morgenstern B (2004). AUGUSTUS: a web server for gene finding in eukaryotes. _Nucleic Acids Research_, 32(suppl_2):W309–W312, <https://doi.org/10.1093/nar/gkh379>
- BLAST+:  
Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J et al. (2009). BLAST+: architecture and applications. _BMC Bioinformatics_, 10(1):421, <https://dx.doi.org/10.1186%2F1471-2105-10-421>
- SNAP:  
Korf I (2004). Gene finding in novel genomes. _BMC Bioinformatics_, 5(1):59, <https://doi.org/10.1186/1471-2105-5-59> 
