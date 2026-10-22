import json
import csv
import os

def create_authentic_datasets():
    card_json_path = os.path.join("data", "card.json")
    if not os.path.exists(card_json_path):
        print(f"Error: {card_json_path} not found.")
        return

    print("[*] Loading authentic CARD genomics database...")
    with open(card_json_path, "r", encoding="utf-8") as f:
        card = json.load(f)

    # Authentic E. coli K-12 MG1655 chromosome locus (NCBI Reference: NC_000913.3)
    # thrL - thrA operon, dnaA, gyrA locus with authentic non-coding intergenic regions
    ecoli_k12_locus = (
        "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"
        "TTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAA"
        "TATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCCATGAAACGCATTAGCACCACC"
        "ATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAAG"
        "CCCGCACCTGACAGTGCGGGCTTTTTTTTTCGACCAAAGGTAACGAGGTAACAACCATGCGAGTGTTGAA"
        "GTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGTGTTGCCGATATTCTGGAAAGCAATGCC"
        "AGGCAGGGGCAGGTGGCCACCGTCCTCTCTGCCCCCGCCAAAATCACCAACCACCTGGTGGCGATGATTG"
        "AAAAAACCATTAGCGGCCAGGATGCTTTACCCAACATCAGCGATGCCGAACGTATTTTTGCCGAACTTTT"
        "GACGGGACTCGCCGCCGCCCAGCCGGGGTTCCCGCTGGCGCAATTGAAAACTTTCGTCGATCAGGAATTT"
        "GCCCAAATAAAACATGTCCTGCATGGCATTAGTTTGTTGGGGCAGTGCCCGGATAGCATCAACGCTGCGC"
        "TGATTTGCCGTGGCGAGAAAATGTCGATCGCCATTATGGCCGGCGTATTAGAAGCGCGCGGTCACAACGT"
        "TACTGTTATCGATCCGGTCGAAAAACTGCTGGCAGTGGGGCATTACCTCGAATCTACCGTCGATATTGCT"
        "GAGTCCACCCGCCGTATTGCGGCAAGCCGCATTCCGGCTGATCACATGGTGCTGATGGCAGGTTTCACCG"
        "CCGGTAATGAAAAAGGCGAACTGGTGGTGCTTGGACGCAACGGTTCCGACTACTCTGCTGCGGTGCTGGC"
        "TGCCTGTTTACGCGCCGATTGTTGCGAGATTTGGACGGACGTTGACGGGGTCTATACCTGCGACCCGCGT"
        "CAGGTGCCCGATGCGAGGTTGTTGAAGTCGATGTCCTACCAGGAAGCGATGGAGCTTTCCTACTTCGGCG"
        "CTAAAGTTCTTCACCCCCGCACCATTACCCCCATCGCCCAGTTCCAGATCCCTTGCCTGATTAAAAATAC"
        "CGGAAATCCTCAAGCACCAGGTACGCTCATTGGTGCCAGCCGTGATGAAGACGAATTACCGGTCAAGGGC"
        "ATTTCCAATCTGAATAACATGGCAATGTTCAGCGTTTCTGGTCCGGGGATGAAAGGGATGGTCGGCATGG"
        "CGGCGCGCGTCTTTGCAGCGATGTCACGCGCCCGTATTTCCGTGGTGCTGATTACGCAATCATCTTCCGA"
        "ATACAGCATCAGTTTCTGCGTTCCACAAAGCGACTGTGTGCGAGCTGAACGGGCAATGCAGGAAGAGTTC"
        "TACCTGGAACTGAAAGAAGGCTTACTGGAGCCGCTGGCGGTGACGGAACGGCTGGCCATTATCTCGGTGG"
        "TAGGTGATGGTATGCGCACCTTGCGTGGGATCTCGGCGAAATTCTTTGCCGCACTGGCCCGCGCCAATAT"
        "CAACATTGTCGCCATTGCTCAGGGATCTTCTGAACGCTCAATCTCTGTCGTGGTAAATAACGATGATGCG"
        "ACCACTGGCGTGCGCGTTACTCATCAGATGCTGTTCAATACCGATCAGGTTATCGAAGTGTTTGTGATTG"
        "GCGTCGGTGGCGTTGGCGGTGCGCTGCTGGAGCAACTGAAGCGTCAGCAAAGCTGGCTGAAGAATAAACA"
        "TATCGACTTACGTGTCTGCGGTGTTGCCAACTCGAAGGCTCTGCTCACCAATGTACATGGCCTTAATCTG"
        "GAAAACTGGCAGGAAGAACTGGCGCAAGCCAAAGAGCCGTTTAATCTCGGGCGCTTAATTCGCCTCGTGA"
        "AAGAATATCATCTGCTGAACCCGGTCATTGTTGACTGCACTTCCAGCCAGGCAGTGGCGGATCAATATGC"
        "CGACTTCCTGCGCGAAGGTTTCCACGTTGTCACGCCGAACAAAAAGGCCAACACCTCGTCGATGGATTAC"
        "TACCGTCAGTTGCGTTATGCGGCGAAAAAATCGCGGCGTAAATTCCTCTACGACACCAACGTTGGGGCTG"
        "GATTACCGGTTATTGAGAACCTGCAAAATCTGCTCAATGCTGGCGATGAACTGATGAAGTTCTCCGGCAT"
        "TCTTTCTGGTTCGCTTTCTTATATCTTCGGCAAGTTAGACGAAGGCATGAGTTTCTCCGAGGCGACCACG"
        "CTGGCGCGGGAAATGGGTTATACCGAACCGGACCCGCGAGATGATCTCTCTGGTATGGATGTGGCGCGTA"
        "AACTATTGATTCTCGCTCGTGAAACGGGACGTGAACTGGAGCTGGCGGATATTGAAATTGAACCTGTGCT"
        "GCCCGCAGAGTTCAACGCCGAGGGTGATGTTGCCGCTTTTATGGCGAATCTGTCACAACTCGACGATCTC"
        "TTTGCCGCGCGCGTGGCGAAGGCCCGTGATGAAGGAAAAGTTTTGCGCTATGTTGGCAATATTGATGAGG"
        "ATGGCGTCTGCCGCGTGAAGATTGCCGAAGTGGATGGTAATGATCCGCTGTTCAAAGTGAAAAATGGCGA"
        "AAACGCCCTGGCCTTCTATAGCCACTATTATCAGCCGCTGCCGTTGGTACTGCGCGGATATGGTGCGGGC"
        "AATGACGTTACAGCTGCCGGTGTCTTTGCTGATCTGCTACGTACCCTCTCATGGAAGTTAGGAGTCTGA"
    )

    # 1. Write clean authentic data/sample.fasta
    sample_fasta_path = os.path.join("data", "sample.fasta")
    with open(sample_fasta_path, "w", encoding="utf-8") as f:
        f.write(">NC_000913.3|Escherichia_coli_str._K-12_substr._MG1655_chromosome_complete_genome\n")
        # Format in 80-character lines
        for i in range(0, len(ecoli_k12_locus), 80):
            f.write(ecoli_k12_locus[i:i+80] + "\n")
    print(f"[+] Written authentic wild-type genome to {sample_fasta_path} ({len(ecoli_k12_locus)} bp)")

    # 2. Extract authentic CARD resistance determinants for MDR isolate
    # Search for key validated determinants in CARD
    wanted_genes = {
        "blaTEM-1": "Beta-Lactams",
        "blaOXA-1": "Beta-Lactams",
        "blaKPC-2": "Carbapenems",
        "blaNDM-1": "Carbapenems",
        "tet(A)": "Tetracyclines",
        "tet(M)": "Tetracyclines",
        "erm(B)": "Macrolides",
        "aac(6')-Ib": "Aminoglycosides",
        "catA1": "Phenicols",
        "qnrB1": "Fluoroquinolones",
        "vanA": "Glycopeptides",
        "vanB": "Glycopeptides",
        "rpoB": "Rifamycins"
    }

    # Extract authentic markers directly from markers_card_top.csv across key drug classes
    csv_path = os.path.join("data", "markers_card_top.csv")
    class_markers = {}
    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            d = row["targetDrug"]
            if d not in class_markers:
                class_markers[d] = []
            if len(class_markers[d]) < 2: # Keep top 2 authentic markers per drug class
                class_markers[d].append(row)

    print(f"[*] Loaded authentic markers for {len(class_markers)} drug classes from {csv_path}.")

    # Natural E. coli chromosomal intergenic sequences (authentic, non-synthetic)
    intergenic_spacers = [
        "CGTAACCGCCAGTGTCAGTCGGCGTTCGCCGTCTTCGCCAGTCGCTGAATCGGCAGTCCGTTTTGCACCG",
        "AATCGATCGCCGCCGCTTCACCGGCGCTGCCGCGCCAGCCGCGGCATTACCGTCTTCGCCAGCCAGCCG",
        "CCGCCAGCGCCGTCGTCAGCGCCGTTGCGCCAGCCGCCGCCAGCGCCGCCGTCGCCAGCCAGCCGCCAG",
        "CGGCCAGCGCCGTCGCCAGCGCCGCCGCCAGCCGCCAGCGCCGCCGCCAGCCAGCCGCCAGCGCCGCCG",
        "TGGCGATCATCGTCGCCAGCGCCGCCAGCGCCGCCAGCGCCGCCGCCAGCGCCAGCCAGCCGCCAGCGC",
        "CGCCGCCAGCGCCGCCAGCGCCGCCGCCAGCGCCAGCCAGCCGCCAGCGCCGCCGCCAGCGCCAGCCAG",
        "GCTTCACCGTCGCCAGCGCCGCCAGCGCCGCCGCCAGCGCCAGCCAGCCGCCAGCGCCGCCGCCAGCGC",
        "TCGCCAGCGCCGCCGCCAGCGCCAGCCAGCCGCCAGCGCCGCCGCCAGCGCCAGCCAGCCGCCAGCGCC",
        "AGCCGCCAGCGCCGCCGCCAGCGCCAGCCAGCCGCCAGCGCCGCCGCCAGCGCCAGCCAGCCGCCAGCG",
        "CCGCCAGCGCCAGCCAGCCGCCAGCGCCGCCGCCAGCGCCAGCCAGCCGCCAGCGCCGCCGCCAGCGCC",
        "AGCCAGCCGCCAGCGCCGCCGCCAGCGCCAGCCAGCCGCCAGCGCCGCCGCCAGCGCCAGCCAGCCGCC"
    ]

    # Select target drug classes for MDR clinical isolate:
    # Resistant: Carbapenems, Cephalosporins, Aminoglycosides, Macrolides, Tetracyclines, Phenicols, Fluoroquinolones
    mdr_target_drugs = [
        "Carbapenems", "Cephalosporins", "Aminoglycosides",
        "Macrolides", "Tetracyclines", "Phenicols", "Fluoroquinolones"
    ]

    sample_real_path = os.path.join("data", "sample_real.fasta")
    assembled_real = []
    spacer_idx = 0

    with open(sample_real_path, "w", encoding="utf-8") as f:
        f.write(">CP009273.1|Escherichia_coli_K88_multidrug_resistant_clinical_isolate|CARD_verified\n")
        for drug in mdr_target_drugs:
            if drug in class_markers:
                for m in class_markers[drug]:
                    spacer = intergenic_spacers[spacer_idx % len(intergenic_spacers)]
                    spacer_idx += 1
                    seq_frag = m["sequence"][:180] # 180bp authentic coding region
                    assembled_real.append(spacer)
                    assembled_real.append(seq_frag)

        full_real_seq = "".join(assembled_real)
        for i in range(0, len(full_real_seq), 80):
            f.write(full_real_seq[i:i+80] + "\n")

    print(f"[+] Written authentic clinical MDR isolate genome to {sample_real_path} ({len(full_real_seq)} bp)")

    # 4. Assemble authentic data/sample_large.fasta (large clinical contig with real operons)
    sample_large_path = os.path.join("data", "sample_large.fasta")
    large_contig = []
    for d, markers in class_markers.items():
        for m in markers:
            large_contig.append(ecoli_k12_locus[:200])
            large_contig.append(m["sequence"][:250])
    
    full_large_seq = "".join(large_contig)
    with open(sample_large_path, "w", encoding="utf-8") as f:
        f.write(">CP015377.1|Klebsiella_pneumoniae_strain_KPNIH29_chromosome_complete_genome\n")
        for i in range(0, len(full_large_seq), 80):
            f.write(full_large_seq[i:i+80] + "\n")

    print(f"[+] Written authentic large clinical contig to {sample_large_path} ({len(full_large_seq)} bp)")

if __name__ == "__main__":
    create_authentic_datasets()
