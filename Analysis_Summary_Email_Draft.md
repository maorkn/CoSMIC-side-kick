
**Subject**: CoSMIC Analysis: Comparison of V3V4 vs. Full-Length 16S MAG Recruitment for Agelas and Oscarella

Hi [Name],

I've completed the CoSMIC pipeline analysis for both *Agelas* and *Oscarella*, comparing the results when using **V3V4 extracted amplicons** versus **full-length 16S sequences** for MAG recruitment.

**Attached:**
1.  `CoSMIC_V3V4_vs_FullLength_Comparison.md`: A detailed report highlighting the specific functional discrepancies.
2.  `Agelas_Oscarella_Runs.zip`: Contains the full output folders for all 4 pipeline runs (including the LLM-ready functional reports for each).

### **Methods Summary (This Run)**
For the February 12, 2026 CoSMIC Sidekick comparison run, the Agelas and Oscarella assemblies were first organized into MAG-only directories with `setup_agelas_oscarella.sh`, and the four amplicon FASTA inputs (Agelas/Oscarella × V3V4/full-length) were converted into CoSMIC-compatible CSV tables with `fasta_to_csv.py`, using a uniform `dummy_abundance` column because this comparison was based on sequence recruitment rather than sample-resolved abundance. The four analyses were then executed with `run_all.sh`, which called `cosmic_sidekick.py run` and `cosmic_sidekick.py report` with `config_agelas_V3V4.yaml`, `config_agelas_full.yaml`, `config_oscarella_V3V4.yaml`, and `config_oscarella_full.yaml`. Within CoSMIC, rRNA loci were extracted from MAGs with Barrnap, metabarcoding sequences were mapped to the extracted rRNA references with `vsearch --usearch_global` at a **95% identity threshold** (Rognes et al., 2016), and only MAGs with at least one CoSMIC hit were annotated with Prokka (Seemann, 2014). Final V3V4-versus-full-length differences were derived from the resulting `cosmic_llm_report.md` outputs using `compare_mags.py` to identify MAGs unique to each run and `extract_functions_v2.py` to extract the dominant Prokka-derived product summaries from V3V4-only MAGs. The optional `minimap2`, eggNOG, and SILVA enrichment paths in CoSMIC were not used for this comparison.

### **Key Findings**

**1. Agelas: Full-Length Provides Superior Coverage**
*   **Mapping Success**: The full-length sequences recruited significantly more ASVs (**78%**) compared to V3V4 (**62%**).
*   **Erroneous V3V4 Functions**: Using only V3V4 resulted in the "erroneous" detection of 10 MAGs not supported by the full-length analysis. This would have led us to incorrectly assume the community had specific capabilities for **xenobiotic degradation** (e.g., *Caprolactamase*, *Caffeine dehydrogenase*) and **fatty acid metabolism**.
*   **Missed Functions**: The V3V4 run completely missed key sponge symbionts like **Poribacteria** and **Candidatus Binatia**, representing a major loss of functional insight.

**2. Oscarella: V3V4 Introduces False Positives**
*   **Erroneous Recruitment**: The V3V4 run recruited an additional MAG (*CATOTD01*) that was rejected by the full-length analysis.
*   **Functional Artifacts**: This false positive introduced specific stress response functions (**Chaperonin GroEL/GroES**) and **sialic acid biosynthesis** pathways that are likely not present in the actual amplicon-associated population.

**Conclusion**: The analysis strongly suggests that using **full-length 16S sequences** yields a more robust and accurate functional profile for these sponge communities, avoiding the specific functional artifacts and false positives seen with the V3V4 region.

### **Software / Method References**
*   Rognes T, Flouri T, Nichols B, Quince C, Mahé F. 2016. **VSEARCH: a versatile open source tool for metagenomics.** *PeerJ* 4:e2584.
*   Seemann T. 2014. **Prokka: rapid prokaryotic genome annotation.** *Bioinformatics* 30(14):2068-2069.
*   Barrnap was used for rRNA prediction via the `barrnap` software by Torsten Seemann; unlike `vsearch` and Prokka, this draft does not cite a primary journal article for Barrnap.

Best,

[Your Name]
