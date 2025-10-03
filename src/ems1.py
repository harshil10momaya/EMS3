from motif_finder import MotifFinderBase
from utils import Params, decode_motif
from typing import Dict, List, Set, Any
import itertools

class Ems1(MotifFinderBase):
    """
    Implements EMS Version 1: Brute-force counting of motif candidates 
    found in ALL sequences (replaces C++ Ems1).
    """

    def __init__(self, input_path: str, l: int, d: int, params: Params):
        super().__init__("ems1", input_path, l, d, params)
        
        # Equivalent of A2 (motif counts) and A1 (last sequence ID matched)
        self.motif_counts: Dict[str, int] = {} 
        self.last_seq_match: Dict[str, int] = {}

    def _generate_neighborhood(self, kmer: str, errors_remaining: int, length: int, current_motif_set: Set[str]):
        """
        Recursively generates the (l, d) neighborhood (Deletions, Substitutions, Insertions).
        """
        current_len = len(kmer)

        if errors_remaining == 0:
            if current_len == length:
                current_motif_set.add(kmer)
            return

        # Pruning check
        if current_len < length - errors_remaining or current_len > length + errors_remaining:
             return

        # 1. Deletion (del)
        if current_len > 0 and current_len >= length - errors_remaining + 1:
            for j in range(current_len):
                new_kmer = kmer[:j] + kmer[j+1:]
                self._generate_neighborhood(new_kmer, errors_remaining - 1, length, current_motif_set)
        
        # 2. Substitution (sub)
        for j in range(current_len):
            original_char = kmer[j]
            for k in range(self.domain_size):
                new_char_code = str(k)
                if new_char_code != original_char:
                    new_kmer = kmer[:j] + new_char_code + kmer[j+1:]
                    self._generate_neighborhood(new_kmer, errors_remaining - 1, length, current_motif_set)
                    
        # 3. Insertion (ins)
        if current_len <= length + errors_remaining - 1:
            for j in range(current_len + 1):
                for k in range(self.domain_size):
                    new_char_code = str(k)
                    new_kmer = kmer[:j] + new_char_code + kmer[j:]
                    self._generate_neighborhood(new_kmer, errors_remaining - 1, length, current_motif_set)


    def search(self):
        n_reads = len(self.reads)
        
        for i, seq in enumerate(self.reads):
            current_seq_id = i + 1
            print(f"Processing sequence {i+1}...")
            
            sequence_candidates: Set[str] = set()

            for q in range(-self.d, self.d + 1):
                k = self.l + q
                
                if k <= 0 or k > len(seq): continue

                for start_index in range(len(seq) - k + 1):
                    kmer = seq[start_index : start_index + k]
                    
                    self._generate_neighborhood(kmer, self.d, self.l, sequence_candidates)

            for candidate_motif in sequence_candidates:
                if self.last_seq_match.get(candidate_motif) != current_seq_id:
                    self.last_seq_match[candidate_motif] = current_seq_id
                    self.motif_counts[candidate_motif] = self.motif_counts.get(candidate_motif, 0) + 1
                    
            print(f"Done. Current candidate pool size: {len(self.motif_counts)}")

        self.motifs = [motif_encoded for motif_encoded, count in self.motif_counts.items() if count == n_reads]
        self.motifs.sort()