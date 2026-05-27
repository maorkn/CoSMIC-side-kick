# Methods

For this comparison, Agelas and Oscarella MAG collections were unpacked and reorganized into MAG-only directories using `setup_agelas_oscarella.sh`. Full-length and V3V4 16S amplicon (CoSMIC genrated output) FASTA files were then converted into CoSMIC-compatible CSV inputs with `fasta_to_csv.py`, This produced four input tables: Agelas full-length (100 ASVs), Agelas V3V4 (100 ASVs), Oscarella full-length (50 ASVs), and Oscarella V3V4 (50 ASVs).

The four CoSMIC Sidekick analyses were executed with `run_all.sh`, which called `cosmic_sidekick.py run` and `cosmic_sidekick.py report` using the run-specific configuration files `config_agelas_V3V4.yaml`, `config_agelas_full.yaml`, `config_oscarella_V3V4.yaml`, and `config_oscarella_full.yaml`. In each run, `cosmic_sidekick.py` first used Barrnap to detect 16S and 18S rRNA loci in the MAG assemblies and to generate MAG-linked rRNA reference files (`barrnap_rrna_mapping.csv` and `barrnap_rrna_sequences.fasta`). Amplicon sequences were then recruited to these rRNA references with `vsearch --usearch_global` at a 95% identity threshold and both strand orientations, as specified in all four configuration files. Only MAGs with at least one retained CoSMIC hit were passed forward for annotation with Prokka. The resulting MAG-level and ASV-level summaries, including Prokka-derived product, EC, and COG annotations, were compiled with `cosmic_sidekick.py report`.
Run-specific differences between V3V4 and full-length recruitment were evaluated separately for Agelas and Oscarella using `compare_mags.py`, which parsed the corresponding `cosmic_llm_report.md` files and identified MAGs unique to each run type. Functional discrepancies associated with V3V4-only recruits were then summarized with `extract_functions_v2.py`, which extracted the dominant Prokka product annotations from those MAG-specific report sections. These post-processing steps were used to generate the comparative interpretation summarized in `CoSMIC_V3V4_vs_FullLength_Comparison.md`.

## References

- Rognes T, Flouri T, Nichols B, Quince C, Mahe F. 2016. VSEARCH: a versatile open source tool for metagenomics. *PeerJ* 4:e2584.
- Seemann T. 2014. Prokka: rapid prokaryotic genome annotation. *Bioinformatics* 30(14):2068-2069.
- Seemann T. Barrnap. Software repository: `https://github.com/tseemann/barrnap`.
