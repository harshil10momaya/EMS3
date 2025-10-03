from motif_finder import MotifFinderBase
from utils import Params, WILDCARD_CODE, decode_motif
from motif_tree import MotifTreeBase, MotifTreeFast, MotifTreeSimple
from typing import Dict, List, Optional
import math

class Ems2(MotifFinderBase):
    """
    Implements EMS Version 2 (ems2 and ems2m): Iterative intersection of motif 
    candidates using a Motif Tree (replaces C++ Ems2).
    """
    
    TREE_MAP = {
        "ems2": MotifTreeFast,   
        "ems2m": MotifTreeSimple 
    }

    def __init__(self, input_path: str, l: int, d: int, params: Params, version: str):
        super().__init__(version, input_path, l, d, params)
        
        TreeClass = self.TREE_MAP.get(version)
        if not TreeClass: raise ValueError(f"Unknown version for Ems2: {version}")
             
        self.main_tree: Optional[MotifTreeBase] = None
        self.TreeClass = TreeClass
        self.kmer: str = ""
        self.l_target = str(WILDCARD_CODE) 
        self.leftmost = False
        self.rightmost = False

    # The recursive neighborhood generation functions (delta/sigma/alpha decomposition)

    def _gen_nbrhood_3(self, alpha: int, tree: MotifTreeBase):
        """Handles insertions ('*')."""
        current_kmer = self.kmer 
        current_len = len(current_kmer)
        
        if alpha == 0:
            # Base Case: Motif is complete. Filter out deletion markers ('-')
            t = "".join(c for c in current_kmer if c != '-')
            t = t.replace('*', self.l_target)
            
            if len(t) == self.l:
                tree.insert(t)
            return

        # Recursive Step: Insertion
        for j in range(current_len + 1):
            if j < current_len and current_kmer[j] == '*': continue
            
            # Temporarily insert '*'
            self.kmer = current_kmer[:j] + '*' + current_kmer[j:]
            self._gen_nbrhood_3(alpha - 1, tree) 
            self.kmer = current_kmer # Backtrack


    def _gen_nbrhood_2(self, sigma: int, alpha: int, tree: MotifTreeBase):
        """Handles substitutions ('*')."""
        current_kmer = self.kmer
        current_len = len(current_kmer)

        if sigma == 0:
            self._gen_nbrhood_3(alpha, tree)
            return

        # Recursive Step: Substitution
        for j in range(current_len):
            if current_kmer[j] == '-': continue
            
            t = current_kmer[j]
            # Substitute with '*'
            self.kmer = current_kmer[:j] + '*' + current_kmer[j+1:]
            self._gen_nbrhood_2(sigma - 1, alpha, tree)
            self.kmer = current_kmer[:j] + t + current_kmer[j+1:] # Backtrack


    def _gen_nbrhood(self, delta: int, sigma: int, alpha: int, tree: MotifTreeBase):
        """Handles deletions ('-')."""
        current_kmer = self.kmer
        current_len = len(current_kmer)
        
        if delta == 0:
            self._gen_nbrhood_2(sigma, alpha, tree)
            return

        # Recursive Step: Deletion
        for j in range(current_len):
            original_char = current_kmer[j]
            
            # Substitute with deletion marker '-'
            self.kmer = current_kmer[:j] + '-' + current_kmer[j+1:]
            self._gen_nbrhood(delta - 1, sigma, alpha, tree)
            self.kmer = current_kmer[:j] + original_char + self.kmer[j+1:] # Backtrack, fixing length


    def _gen_all(self, seq: str, tree: MotifTreeBase):
        """Iterates through all possible k-mer lengths and error partitions."""
        m = len(seq)
        
        for q in range(-self.d, self.d + 1):
            k = self.l + q 
            
            if k <= 0 or k > m: continue

            for delta in range(max(0, q), math.floor((self.d + q) / 2) + 1):
                alpha = delta - q 
                sigma = self.d - alpha - delta
                
                for i in range(m - k + 1):
                    self.kmer = seq[i : i + k]
                    
                    self._gen_nbrhood(delta, sigma, alpha, tree)

    def search(self):
        if not self.reads: return
            
        print(f"Processing sequence 0...")
        self.main_tree = self.TreeClass(self.l, self.motifs, "main")
        self._gen_all(self.reads[0], self.main_tree)
        
        for i in range(1, len(self.reads)):
            print(f"Processing sequence {i}...")
            
            tmp_motifs: List[str] = []
            tmp_tree = self.TreeClass(self.l, tmp_motifs, "tmp") 
            
            self._gen_all(self.reads[i], tmp_tree)
            self.main_tree.intersect(tmp_tree)
            
            if not self.main_tree.root['children']: 
                self.motifs.clear()
                print("No common motifs found after intersection. Stopping.")
                return

        self.main_tree.traverse()
        self.motifs.sort()