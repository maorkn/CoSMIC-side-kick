# CoSMIC Sidekick: V3V4 vs. Full-Length 16S Comparison

This report summarizes the differences in MAG recruitment and functional potential when using **V3V4 amplicons** versus **full-length 16S sequences** for the *Agelas* and *Oscarella* datasets.

## 1. Summary Statistics

### Agelas Dataset
| Metric | V3V4 Run | Full-Length Run | Difference |
| :--- | :--- | :--- | :--- |
| **Total ASVs** | 100 | 100 | - |
| **Mapped ASVs** | **62** (62%) | **78** (78%) | **+16** (Full-length recruited more) |
| **Recruited MAGs** | 49 | 50 | Similar count, but different composition |

### Oscarella Dataset
| Metric | V3V4 Run | Full-Length Run | Difference |
| :--- | :--- | :--- | :--- |
| **Total ASVs** | 50 | 50 | - |
| **Mapped ASVs** | 2 (4%) | 2 (4%) | Identical count |
| **Recruited MAGs** | **2** | **1** | **V3V4 recruited one extra MAG** |

---

## 2. Erroneously Assumed Functions (V3V4 Artifacts)

Using only the V3V4 region resulted in the recruitment of several MAGs that were **not** supported by the full-length analysis. If relyed upon, these would lead to erroneous assumptions about the community's functional potential.

### Oscarella: Erroneous Functional Assumptions
The V3V4 analysis uniquely recruited **MAG CATOTD01**. This MAG was absent in the full-length analysis.

**Functions erroneously assumed present:**
*   **Restriction-Modification Systems**: *Type II methyltransferase M.EcoRI*, *Histidine N-alpha-methyltransferase*.
*   **Stress Response**: *Chaperonin GroEL* and *GroES* (indicating potential heat shock or stress adaptation).
*   **Sialic Acid Biosynthesis**: *N,N'-diacetyllegionaminic acid synthase* (often involved in cell surface modification or host interaction).

### Agelas: Erroneous Functional Assumptions
The V3V4 analysis uniquely recruited **10 MAGs** that were not found in the full-length run.

**Key erroneous MAGs and their functions:**
*   **MAG CATQHL01**:
    *   **Xenobiotic Degradation**: *Caprolactamase* (degradation of synthetic compounds) and *Caffeine dehydrogenase*.
    *   **Fatty Acid Metabolism**: *3-oxoacyl-[acyl-carrier-protein] reductase FabG*.
*   **MAG odAgeOroi1 Series (Gammaproteobacteria, Alphaproteobacteria, etc.)**:
    *   **Osmoregulation**: *Ectoine hydrolase* (salt tolerance mechanism).
    *   **Antibiotic Resistance**: *Macrolide export ATP-binding/permease protein MacB*.
    *   **Metabolism**: *Carnitine monooxygenase* and *Succinate-glutarate CoA-transferase*.

**Impact:** Relying on V3V4 would falsely suggest the community has specific capabilities for **degrading caffeine/caprolactam**, **synthesizing sialic acids**, and specific **stress/salt tolerance mechanisms** that are likely not present (or at least not associated with these specific taxa) when the full gene context is considered.

---

## 3. Missed Functions (False Negatives in V3V4)

The V3V4 analysis failed to detect several MAGs that *were* present in the full-length analysis. This represents a loss of functional insight.

### Agelas: Missed Functions
The V3V4 analysis missed **11 MAGs**, including:
*   **CATQHS01** & **CATQHW01**
*   **GCA_958295495.1** (Large genome, >5000 CDS)
*   Several **odAgeOroi1** members (*Aggregatilineales*, *Candidatus Binatia*, *Poribacteria*, *Vicinamibacteria*)

**Functional Loss:**
*   **Candidatus Binatia**: Often associated with complex degradation capabilities in sponges.
*   **Poribacteria**: A key sponge symbiont group often involved in secondary metabolite production. Missing this group significantly alters the understanding of the sponge holobiont.
*   **Aggregatilineales**: Likely involved in filament formation or sequence-specific behaviors.

---

## 4. Conclusion

*   **V3V4 Analysis Risks**: It introduces **false positives** (e.g., *CATOTD01* in Oscarella, *CATQHL01* in Agelas) leading to unsupported functional claims (e.g., caffeine degradation, sialic acid synthesis). It also suffers from **false negatives**, missing key symbionts like *Poribacteria*.
*   **Full-Length Superiority**: The full-length analysis provided broader coverage (78% vs 62% mapped in Agelas) and recruited a different, likely more accurate set of MAGs, avoiding the specific functional artifacts highlighted above.
