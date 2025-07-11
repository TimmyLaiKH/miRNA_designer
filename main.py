#!/usr/bin/env python3
"""
Advanced Artificial miRNA Designer for Alzheimer's Disease
Targeting APP, MAPT, and PSEN1 Genes with Seed-Type Focused Design
"""

import re
import logging
import time
import subprocess
import tempfile
import os
import csv
from collections import Counter
from Bio import Entrez, SeqIO
from Bio.SeqUtils import MeltingTemp as mt
import requests

# Configuration
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
Entrez.email = "your_email@domain.edu"  # REPLACE WITH YOUR EMAIL

# Maximum possible binding score (theoretical maximum)
MAX_BINDING_SCORE = 11.5  # 5 (8mer) + 4.5 (9 perfect 3' matches) + 2 (100% AU flank)

# Alzheimer's gene targets with RefSeq IDs
GENES = {
    "APP": "NM_000484.4",  # Amyloid precursor protein
    "MAPT": "NM_001377265.1",  # Microtubule-associated protein tau
    "PSEN1": "NM_000021.4"  # Presenilin 1
}

# Essential neuronal genes to avoid off-targets
ESSENTIAL_NEURONAL_GENES = {
    "BDNF": "NM_001143805.1",  # Brain-derived neurotrophic factor
    "GRIN1": "NM_007327.4",  # Glutamate ionotropic receptor NMDA type subunit 1
    "SYP": "NM_003179.3",  # Synaptophysin
    "SNAP25": "NM_003081.4",  # Synaptosomal-associated protein 25
    "DLG4": "NM_001365.5"  # Postsynaptic density protein 95
}

# ViennaRNA accessibility parameters
ACCESS_CONTEXT_SIZE = 70  # Nucleotides around target site to consider for folding

# Seed type definitions
SEED_TYPES = {
    "8mer": {
        "length": 8,
        "description": "Positions 1-8 with adenine at position 1",
        "guide_positions": (0, 8),
        "target_requirements": lambda site: site[0] == 'A'
    },
    "7mer-A1": {
        "length": 7,
        "description": "Positions 1-7 with adenine at position 1",
        "guide_positions": (0, 7),
        "target_requirements": lambda site: site[0] == 'A'
    },
    "7mer-m8": {
        "length": 7,
        "description": "Positions 2-8",
        "guide_positions": (1, 8),
        "target_requirements": lambda site: True  # No additional requirements
    },
    "6mer": {
        "length": 6,
        "description": "Minimal seed, positions 2-7",
        "guide_positions": (1, 7),
        "target_requirements": lambda site: True  # No additional requirements
    }
}


class amiRNADesigner:
    def __init__(self, target_genes, essential_genes=None):
        self.target_genes = target_genes
        self.essential_genes = essential_genes or {}
        self.utr_seqs = {}
        self.essential_utrs = {}

        # Fetch all sequences
        self.fetch_all_sequences()

        # Check for ViennaRNA
        self.has_rnafold = self.check_rnafold_availability()
        logging.info(f"RNAfold available: {self.has_rnafold}")

    def check_rnafold_availability(self):
        """Check if RNAfold is available in PATH with better error handling"""
        try:
            result = subprocess.run(["RNAfold", "--version"],
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    text=True)
            return result.returncode == 0
        except Exception as e:
            logging.warning(f"RNAfold check failed: {e}")
            return False

    def safe_efetch(self, **kwargs):
        """Robust sequence fetcher with retry logic"""
        attempts = 0
        while attempts < 3:
            try:
                handle = Entrez.efetch(**kwargs)
                return handle
            except Exception as e:
                logging.warning(f"NCBI error ({kwargs['id']}): {e}, retrying...")
                time.sleep(5)
                attempts += 1
        logging.error(f"Failed to fetch {kwargs['id']} after 3 attempts")
        return None

    def fetch_3utr(self, refseq_id):
        """Retrieve 3'UTR sequence from RefSeq ID with robust parsing"""
        handle = self.safe_efetch(db="nucleotide", id=refseq_id, rettype="gb", retmode="text")
        if not handle:
            return ""

        try:
            record = SeqIO.read(handle, "genbank")
            cds_end = 0

            # Find CDS end position
            for feat in record.features:
                if feat.type == "CDS":
                    cds_end = feat.location.end
                    break

            # Extract 3'UTR from CDS end to sequence end
            utr_seq = record.seq[cds_end:]

            logging.debug(f"Fetched {refseq_id}: CDS ends at {cds_end}, UTR length {len(utr_seq)}")
            return str(utr_seq)
        except Exception as e:
            logging.error(f"Error parsing {refseq_id}: {e}")
            return ""

    def fetch_all_sequences(self):
        """Fetch all target and essential gene UTRs"""
        logging.info("Fetching target gene 3'UTRs...")
        for gene, refseq in self.target_genes.items():
            self.utr_seqs[gene] = self.fetch_3utr(refseq)
            logging.info(f"Fetched {gene} ({len(self.utr_seqs[gene])} bp)")

        logging.info("Fetching essential gene 3'UTRs...")
        for gene, refseq in self.essential_genes.items():
            seq = self.fetch_3utr(refseq)
            if not seq:
                logging.warning(f"Using empty sequence for {gene}")
            self.essential_utrs[gene] = seq
            logging.info(f"Fetched {gene} ({len(seq)} bp)")

    def revcomp(self, seq):
        """Reverse complement of DNA sequence"""
        comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        return ''.join(comp.get(base, 'N') for base in reversed(seq.upper()))

    def revcomp_rna(self, seq):
        """Reverse complement for RNA sequences (DNA input, RNA output)"""
        comp = {'A': 'U', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        return ''.join(comp.get(base, 'N') for base in reversed(seq.upper()))

    def find_conserved_seed_sites(self, seed_type):
        """Identify conserved seed regions across all targets for a specific seed type"""
        seed_info = SEED_TYPES[seed_type]
        seed_length = seed_info["length"]
        seed_candidates = {}

        for gene, seq in self.utr_seqs.items():
            if not seq:
                continue

            seeds = {}
            for i in range(len(seq) - seed_length + 1):
                seed = seq[i:i + seed_length]
                # Check if this seed meets requirements for the seed type
                if seed_info["target_requirements"](seed):
                    seeds.setdefault(seed, []).append(i)
            seed_candidates[gene] = seeds

        # Find common seeds across all genes
        common_seeds = set.intersection(*[set(seeds.keys()) for seeds in seed_candidates.values()])
        return common_seeds, seed_length

    def find_binding_sites_for_seed(self, seed, seed_type):
        """Locate all binding sites for a seed across target genes with correct positioning"""
        binding_sites = {gene: [] for gene in self.target_genes}
        seed_info = SEED_TYPES[seed_type]
        seed_length = seed_info["length"]
        start_offset, end_offset = seed_info["guide_positions"]

        # Correct positioning: miRNA seed corresponds to positions 2-8 in guide (1-based)
        # We need to extract target site where seed is at positions 2-8 of the guide's binding region
        # This means the seed should be at the 3' end of the target site
        for gene, utr in self.utr_seqs.items():
            if not utr:
                continue

            start_idx = 0
            while start_idx < len(utr):
                idx = utr.find(seed, start_idx)
                if idx == -1:
                    break

                # CORRECTED: Extract 21nt target site with seed positioned for miRNA binding
                # - The first base of the seed in target should be at position 2 in miRNA (1-based)
                # - This corresponds to position 14 in the 21-nt target site (0-indexed)
                site_start = idx - 13  # miRNA positions 1-21 bind to target positions 21-1
                site_end = site_start + 21

                # Handle boundaries with padding
                if site_start < 0:
                    target_site = 'N' * (-site_start) + utr[0:site_end]
                elif site_end > len(utr):
                    target_site = utr[site_start:] + 'N' * (site_end - len(utr))
                else:
                    target_site = utr[site_start:site_end]

                # Ensure we have exactly 21nt
                target_site = target_site[:21].ljust(21, 'N')

                # Mark seed in uppercase within context
                seed_start = 13  # Position of first seed base in 21-nt site
                seed_end = seed_start + seed_length
                marked_context = (
                        target_site[:seed_start].lower() +
                        target_site[seed_start:seed_end].upper() +
                        target_site[seed_end:].lower()
                )

                binding_sites[gene].append({
                    "position": idx,
                    "context": marked_context,
                    "seed_match": seed,
                    "target_site_21nt": target_site.upper(),
                    "seed_type": seed_type
                })
                start_idx = idx + 1  # Move to next position

        return binding_sites

    def build_consensus_site(self, sites):
        """Build consensus target site from multiple sites"""
        if not sites:
            return ""

        consensus = []
        for i in range(21):  # For each position in 21nt site
            bases = [site[i] for site in sites if site[i] != 'N']
            if bases:
                most_common = Counter(bases).most_common(1)[0][0]
                consensus.append(most_common)
            else:
                consensus.append('N')
        return ''.join(consensus)

    def design_amiRNA(self, binding_sites):
        """Generate optimized amiRNA sequence from target sites"""
        # Collect all 21nt target sites
        all_target_sites = []
        seed_type = None

        for gene_sites in binding_sites.values():
            for site in gene_sites:
                all_target_sites.append(site["target_site_21nt"])
                if seed_type is None:
                    seed_type = site["seed_type"]

        if not all_target_sites or seed_type is None:
            return None, None

        # Build consensus target site
        consensus_site = self.build_consensus_site(all_target_sites)

        # Design guide as reverse complement of target site (RNA sequence)
        amiRNA = self.revcomp_rna(consensus_site)

        # UNCONDITIONALLY SET FIRST BASE TO 'U' FOR ENHANCED RISC LOADING
        if amiRNA:
            amiRNA = 'U' + amiRNA[1:]

        return amiRNA, consensus_site

    def calculate_binding_score(self, amiRNA, target_site, seed_type):
        """
        Calculate binding score based on miRNA-target interaction features
        Scoring factors:
        - Seed type (8mer, 7mer-m8, 7mer-A1, 6mer)
        - 3' pairing quality
        - AU content in flanking region

        Theoretical maximum: 11.5 (5 + 4.5 + 2)
        """
        if not target_site:
            return 0

        seed_info = SEED_TYPES[seed_type]
        start_idx, end_idx = seed_info["guide_positions"]
        seed = amiRNA[start_idx:end_idx]
        rev_target = self.revcomp(target_site)
        score = 0

        # 1. Seed type scoring
        if seed_type == "8mer":
            score += 5
        elif seed_type == "7mer-m8":
            score += 4
        elif seed_type == "7mer-A1":
            score += 3
        elif seed_type == "6mer":
            score += 2

        # 2. 3' pairing quality (positions 13-21)
        three_prime = amiRNA[12:]
        if len(rev_target) >= len(three_prime):
            three_prime_target = rev_target[-len(three_prime):]
            score += 0.5 * sum(a == b for a, b in zip(three_prime, three_prime_target))

        # 3. AU content in flanking region
        if len(target_site) >= 15:
            flanking = target_site[5:15]  # Extract flanking sequence
            au_content = (flanking.count('A') + flanking.count('T')) / len(flanking)
            score += 2 * au_content

        return min(score, MAX_BINDING_SCORE)  # Cap at theoretical maximum

    def calculate_accessibility(self, sequence):
        """
        Calculate target site accessibility using RNAfold
        Returns probability of unpaired state
        """
        if not self.has_rnafold or not sequence:
            return 0.5  # Default value if RNAfold not available

        try:
            # Create temporary input file
            with tempfile.NamedTemporaryFile(mode='w+', delete=False) as tmpin:
                tmpin.write(f">temp_seq\n{sequence}\n")
                tmpin_name = tmpin.name

            # Run RNAfold
            result = subprocess.run(
                ["RNAfold", "--noPS", "-p", "--MEA", "-T", "37", "-i", tmpin_name],
                capture_output=True,
                text=True,
                check=True
            )

            # Parse output for unpaired probabilities
            output_lines = result.stdout.split('\n')
            if len(output_lines) < 3:
                return 0.5

            # Find base pairing probabilities line
            bpp_line = None
            for line in output_lines:
                if line.startswith('.'):
                    bpp_line = line
                    break

            if not bpp_line:
                return 0.5

            # Extract unpaired probabilities
            unpaired_probs = []
            for char in bpp_line.strip():
                if char == '.':
                    unpaired_probs.append(1.0)  # Unpaired
                elif char in {'(', ')'}:
                    unpaired_probs.append(0.0)  # Paired
                else:
                    unpaired_probs.append(0.5)  # Undetermined

            # Calculate accessibility in target region (center of sequence)
            center_start = len(sequence) // 2 - 7
            center_end = center_start + 15
            if center_start < 0 or center_end > len(unpaired_probs):
                return 0.5

            center_probs = unpaired_probs[center_start:center_end]

            return sum(center_probs) / len(center_probs) if center_probs else 0.5

        except (subprocess.CalledProcessError, ValueError, IndexError) as e:
            logging.error(f"RNAfold error: {e}")
            return 0.5
        except Exception as e:
            logging.error(f"Unexpected error in accessibility calculation: {e}")
            return 0.5
        finally:
            try:
                os.unlink(tmpin_name)
            except:
                pass

    def calculate_specificity(self, amiRNA):
        """Predict off-target effects against essential genes"""
        # Use standard seed region (positions 2-8) for off-target prediction
        seed = amiRNA[1:8]
        off_targets = []

        for gene, seq in self.essential_utrs.items():
            if seq and self.revcomp(seq).find(seed) != -1:
                off_targets.append(gene)

        return off_targets

    def optimize_scaffold(self, amiRNA, scaffold="hsa-mir-30a"):
        """Incorporate into biological scaffold with robust fallback"""
        # Predefined scaffold database
        scaffold_db = {
            "hsa-mir-30a": {
                "precursor": "UGUAAACAUCCUCGACUGGAAGCUAGUUUUCUGCAGGUGUUUGCCUAUUGA",
                "mature_start": 0,
                "mature_end": 22
            },
            "hsa-mir-155": {
                "precursor": "UUAAUGCUAAUCGUGAUAGGGGUUUAGGGUUCUCCCUUCAUGCUUGUUCUG",
                "mature_start": 0,
                "mature_end": 23
            },
            "hsa-mir-21": {
                "precursor": "UAGCUUAUCAGACUGAUGUUGAUCAUGGUGUUUUCUCUUAAUCAUAUGUCUGAUGAGUCU",
                "mature_start": 0,
                "mature_end": 22
            }
        }

        # Use scaffold from database
        if scaffold in scaffold_db:
            data = scaffold_db[scaffold]
            precursor = data["precursor"]
            start = data["mature_start"]
            end = data["mature_end"]

            # Replace mature region with amiRNA
            return precursor[:start] + amiRNA + precursor[end:]

        # Fallback: return amiRNA alone
        return amiRNA

    def check_thermodynamics(self, amiRNA):
        """Calculate melting temperature and GC content for RNA"""
        if not amiRNA:
            return 0, 0

        # Convert to RNA for calculation if needed
        rna_seq = amiRNA.replace('T', 'U')
        gc_percent = 100 * (rna_seq.count('G') + rna_seq.count('C')) / len(rna_seq)

        try:
            # RNA-specific calculation
            tm = mt.Tm_NN(rna_seq, nn_table=mt.RNA_NN2, dnac1=0)
        except:
            tm = 60  # Default value if calculation fails
        return gc_percent, tm


# Execution Pipeline
if __name__ == "__main__":
    logging.info("Starting Alzheimer's amiRNA Designer Targeting APP, MAPT, and PSEN1")
    logging.info("Features: Seed-Type Focused Design, Binding Score Calculation, Accessibility Prediction")
    logging.info(f"Maximum Theoretical Binding Score: {MAX_BINDING_SCORE}")
    designer = amiRNADesigner(GENES, ESSENTIAL_NEURONAL_GENES)

    if designer.has_rnafold:
        logging.info("RNAfold detected - accessibility calculations enabled")
    else:
        logging.warning("RNAfold not found in PATH. Accessibility will be estimated")

    # Step 1: Iterate through seed types and find conserved seeds
    all_candidates = []

    for seed_type, seed_info in SEED_TYPES.items():
        logging.info(f"Processing seed type: {seed_type} ({seed_info['description']})")
        common_seeds, seed_length = designer.find_conserved_seed_sites(seed_type)

        if not common_seeds:
            logging.info(f"No conserved {seed_type} seeds found")
            continue

        logging.info(f"Found {len(common_seeds)} conserved {seed_type} seeds")

        for seed in list(common_seeds)[:100]:  # Limit evaluation per seed type
            binding_sites = designer.find_binding_sites_for_seed(seed, seed_type)

            # Skip seeds without binding sites in ALL genes
            if any(len(sites) == 0 for sites in binding_sites.values()):
                continue

            amiRNA, consensus_site = designer.design_amiRNA(binding_sites)
            if not amiRNA:
                continue

            off_targets = designer.calculate_specificity(amiRNA)
            gc_percent, tm = designer.check_thermodynamics(amiRNA)

            total_binding_score = 0
            total_accessibility = 0
            site_count = 0

            for gene_sites in binding_sites.values():
                for site in gene_sites:
                    site["binding_score"] = designer.calculate_binding_score(
                        amiRNA,
                        site["target_site_21nt"],
                        site["seed_type"]
                    )
                    total_binding_score += site["binding_score"]

                    extended_context = site["target_site_21nt"]
                    if len(extended_context) < ACCESS_CONTEXT_SIZE:
                        padding = (ACCESS_CONTEXT_SIZE - len(extended_context)) // 2
                        extended_context = "N" * padding + extended_context + "N" * padding

                    site["accessibility"] = designer.calculate_accessibility(extended_context)
                    total_accessibility += site["accessibility"]
                    site_count += 1

            avg_binding_score = total_binding_score / site_count if site_count else 0
            avg_accessibility = total_accessibility / site_count if site_count else 0.5

            scaffold = designer.optimize_scaffold(amiRNA)

            if not off_targets and 30 <= gc_percent <= 70 and 50 <= tm <= 75:
                all_candidates.append({
                    "seed": seed,
                    "amiRNA": amiRNA,
                    "consensus_site": consensus_site,
                    "scaffold": scaffold,
                    "gc_percent": gc_percent,
                    "tm": tm,
                    "off_targets": off_targets,
                    "binding_sites": binding_sites,
                    "binding_score": avg_binding_score,
                    "accessibility": avg_accessibility,
                    "seed_type": seed_type,
                    "seed_length": seed_length
                })

    if not all_candidates:
        logging.error("No suitable candidates found! Consider:")
        logging.error("1. Relaxing seed type requirements")
        logging.error("2. Targeting only APP and MAPT")
        logging.error("3. Using separate miRNAs for PSEN1")
        exit(1)

    logging.info("\nTop candidates (sorted by seed type priority, binding score, accessibility):")
    # Define seed type priority
    type_priority = {"8mer": 4, "7mer-m8": 3, "7mer-A1": 2, "6mer": 1}
    all_candidates.sort(key=lambda x: (
        -type_priority.get(x["seed_type"], 0),
        -x["binding_score"],
        -x["accessibility"]
    ))

    # Prepare CSV output
    timestamp = time.strftime("%Y%m%d-%H%M%S")
    summary_filename = f"amiRNA_candidates_summary_{timestamp}.csv"
    sites_filename = f"amiRNA_binding_sites_{timestamp}.csv"

    with open(summary_filename, 'w', newline='') as summary_file, open(sites_filename, 'w', newline='') as sites_file:
        summary_writer = csv.writer(summary_file)
        sites_writer = csv.writer(sites_file)

        # Write headers
        summary_writer.writerow([
            "Candidate_ID", "Seed_Type", "Seed_Length", "Target_Seed_mRNA",
            "Consensus_mRNA_Seed", "Expected_Guide_Seed", "Actual_Guide_Seed",
            "amiRNA_Sequence", "GC_Percent", "Tm", "Avg_Binding_Score",
            "Avg_Accessibility", "Off_Targets", "Scaffold"
        ])

        sites_writer.writerow([
            "Candidate_ID", "Gene", "Site_Position", "Site_Context",
            "Seed_Type", "Binding_Score", "Accessibility"
        ])

        for i, cand in enumerate(all_candidates[:10]):
            seed_info = SEED_TYPES[cand["seed_type"]]
            start_idx, end_idx = seed_info["guide_positions"]
            binding_percent = (cand["binding_score"] / MAX_BINDING_SCORE) * 100

            # Calculate seed positions in consensus site
            seed_start = 13  # Fixed position in target site
            consensus_mrna_seed = cand['consensus_site'][seed_start:seed_start + cand['seed_length']]
            expected_guide_seed = designer.revcomp_rna(consensus_mrna_seed)
            actual_guide_seed = cand['amiRNA'][start_idx:end_idx]

            # Console output
            print(f"\n{'=' * 60}\nCandidate #{i + 1}")
            print(f"Seed Type: {cand['seed_type']} ({cand['seed_length']}-mer)")
            print(f"Target Seed (mRNA): {cand['seed']}")
            print(f"Consensus mRNA Seed: {consensus_mrna_seed}")
            print(f"Expected Guide Seed: {expected_guide_seed}")
            print(f"Actual Guide Seed: {actual_guide_seed}")
            print(f"amiRNA: {cand['amiRNA']}")
            print(f"GC: {cand['gc_percent']:.1f}% | Tm: {cand['tm']:.1f}°C")
            print(f"Binding Score: {cand['binding_score']:.2f}/{MAX_BINDING_SCORE} ({binding_percent:.1f}%)")
            print(f"Accessibility: {cand['accessibility']:.2f}")
            print(f"Off-targets: {cand['off_targets'] or 'None'}")

            scaffold_display = cand['scaffold']
            if len(scaffold_display) > 50:
                scaffold_display = scaffold_display[:50] + "..."
            print(f"Scaffold: {scaffold_display}")

            # Write to summary CSV
            summary_writer.writerow([
                f"Candidate_{i + 1}", cand['seed_type'], cand['seed_length'], cand['seed'],
                consensus_mrna_seed, expected_guide_seed, actual_guide_seed,
                cand['amiRNA'], f"{cand['gc_percent']:.1f}", f"{cand['tm']:.1f}",
                f"{cand['binding_score']:.2f}", f"{cand['accessibility']:.2f}",
                ";".join(cand['off_targets']) if cand['off_targets'] else "None",
                cand['scaffold']
            ])

            # Write binding sites to CSV
            for gene, sites in cand["binding_sites"].items():
                for site in sites:
                    sites_writer.writerow([
                        f"Candidate_{i + 1}", gene, site['position'], site['context'],
                        site['seed_type'], f"{site['binding_score']:.2f}", f"{site['accessibility']:.2f}"
                    ])

            # Print binding sites to console
            print("\nBinding Sites:")
            for gene, sites in cand["binding_sites"].items():
                print(f"  {gene} ({len(sites)} sites):")
                for j, site in enumerate(sites[:2]):  # Show first 2 sites per gene
                    print(f"    Site {j + 1}: Pos {site['position']}")
                    print(f"        Context: {site['context']}")
                    print(f"        Seed Type: {site['seed_type']}")
                    print(
                        f"        Binding: {site['binding_score']:.2f}/{MAX_BINDING_SCORE} | Accessibility: {site['accessibility']:.2f}")
                if len(sites) > 2:
                    print(f"    ... and {len(sites) - 2} more sites")
        print("=" * 60)

    logging.info(f"\nCSV outputs generated:\n- Summary: {summary_filename}\n- Binding sites: {sites_filename}")