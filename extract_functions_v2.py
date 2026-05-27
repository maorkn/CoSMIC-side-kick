
import re

def extract_mag_functions(report_path, target_mags):
    mag_funcs = {}
    current_mag = None
    capture_products = False
    
    with open(report_path, 'r', encoding='utf-8') as f:
        # Read entire file into memory to avoid seeking issues
        lines = f.readlines()
        
    for i, line in enumerate(lines):
        line = line.strip()
        if line.startswith("### MAG "):
            mag_name = line.replace("### MAG ", "")
            if mag_name in target_mags:
                current_mag = mag_name
            else:
                current_mag = None
        
        if current_mag and "- Most frequent annotated products:" in line:
            funcs = []
            # Look at next few lines
            for j in range(1, 10): # Look ahead up to 10 lines
                if i + j >= len(lines): break
                next_line = lines[i+j].strip()
                if next_line.startswith("-"):
                    cleaned = next_line.lstrip("- ").strip()
                    if cleaned and not cleaned.startswith("Most frequent") and not cleaned.startswith("Top EC"):
                        funcs.append(cleaned)
                    if len(funcs) >= 5: break
                elif next_line.startswith("#") or next_line.strip() == "":
                    break
            mag_funcs[current_mag] = funcs

    return mag_funcs

if __name__ == "__main__":
    # Agelas - Unique to V3V4
    agelas_mags = [
        "CATQHL01", 
        "GCA_958295645.1", 
        "odAgeOroi1.Acidimicrobiales_4.1.primary",
        "odAgeOroi1.Alphaproteobacteria_4.1.primary",
        "odAgeOroi1.Caldilineaceae_1.1.primary",
        "odAgeOroi1.Gammaproteobacteria_26.1.primary",
        "odAgeOroi1.Gammaproteobacteria_38.1.primary",
        "odAgeOroi1.Gammaproteobacteria_41.1.primary",
        "odAgeOroi1.Longimicrobiales_1.1.primary",
        "odAgeOroi1.Pseudomonadales_2.1.primary"
    ]
    
    print("=== AGELAS - ERRONEOUSLY DETECTED FUNCTIONS (V3V4 ONLY) ===")
    funcs = extract_mag_functions("run_Agelas_V3V4/cosmic_llm_report.md", set(agelas_mags))
    for mag in sorted(funcs.keys()):
        print(f"\nMAG: {mag}")
        for p in funcs[mag]:
            print(f"  - {p}")

    # Oscarella - Unique to V3V4
    oscarella_mags = ["CATOTD01"]
    print("\n\n=== OSCARELLA - ERRONEOUSLY DETECTED FUNCTIONS (V3V4 ONLY) ===")
    funcs_osc = extract_mag_functions("run_Oscarella_V3V4/cosmic_llm_report.md", set(oscarella_mags))
    for mag in sorted(funcs_osc.keys()):
        print(f"\nMAG: {mag}")
        for p in funcs_osc[mag]:
            print(f"  - {p}")
