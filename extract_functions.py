
import re

def extract_mag_functions(report_path, target_mags):
    mag_funcs = {}
    current_mag = None
    capture = False
    
    with open(report_path, 'r', encoding='utf-8') as f:
        # Create an iterator to manually advance lines
        lines = iter(f)
        for line in lines:
            line = line.strip()
            if line.startswith("### MAG "):
                mag_name = line.replace("### MAG ", "")
                if mag_name in target_mags:
                    current_mag = mag_name
                    capture = True
                else:
                    capture = False
                    current_mag = None
            
            if capture and current_mag and line.startswith("- Most frequent annotated products:"):
                funcs = []
                for _ in range(5): # get top 5
                    try:
                        func_line = next(lines).strip()
                        if func_line.startswith("-"):
                             funcs.append(func_line.lstrip("- "))
                    except StopIteration:
                        break
                mag_funcs[current_mag] = funcs
                capture = False # Reset
                
    return mag_funcs

if __name__ == "__main__":
    # Agelas
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
    
    print("=== Agelas V3V4-Only Functions ===")
    funcs = extract_mag_functions("run_Agelas_V3V4/cosmic_llm_report.md", set(agelas_mags))
    for mag in sorted(funcs.keys()):
        print(f"\nMAG: {mag}")
        for p in funcs[mag]:
            print(f"  - {p}")

    # Oscarella
    oscarella_mags = ["CATOTD01"]
    print("\n=== Oscarella V3V4-Only Functions ===")
    funcs_osc = extract_mag_functions("run_Oscarella_V3V4/cosmic_llm_report.md", set(oscarella_mags))
    for mag in sorted(funcs_osc.keys()):
        print(f"\nMAG: {mag}")
        for p in funcs_osc[mag]:
            print(f"  - {p}")
