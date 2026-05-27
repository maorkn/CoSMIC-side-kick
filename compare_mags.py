
import os
import re

def get_mags_from_report(report_path):
    mags = set()
    try:
        with open(report_path, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith("### MAG "):
                    mags.add(line.strip().replace("### MAG ", ""))
    except Exception as e:
        print(f"Error reading {report_path}: {e}")
    return mags

def compare(organism, v3v4_path, full_path):
    v3v4 = get_mags_from_report(v3v4_path)
    full = get_mags_from_report(full_path)
    
    print(f"=== {organism} Analysis ===")
    print(f"V3V4 MAG count: {len(v3v4)}")
    print(f"Full-Length MAG count: {len(full)}")
    
    # In V3V4 but NOT in Full (Potential False Positives)
    v3v4_only = v3v4 - full
    print(f"\nMAGs ONLY in V3V4 Run (Potential 'Erroneous' Hits):")
    if not v3v4_only:
        print("  <None>")
    for mag in sorted(v3v4_only):
        print(f"  - {mag}")
        
    # In Full but NOT in V3V4 (Potential False Negatives)
    full_only = full - v3v4
    print(f"\nMAGs ONLY in Full-Length Run (Missed by V3V4):")
    if not full_only:
        print("  <None>")
    for mag in sorted(full_only):
        print(f"  - {mag}")
    print("\n")

if __name__ == "__main__":
    compare("Agelas", "run_Agelas_V3V4/cosmic_llm_report.md", "run_Agelas_full/cosmic_llm_report.md")
    compare("Oscarella", "run_Oscarella_V3V4/cosmic_llm_report.md", "run_Oscarella_full/cosmic_llm_report.md")
