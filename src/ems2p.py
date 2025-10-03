from motif_finder import MotifFinderBase
from utils import Params, WILDCARD_CODE, uint64, uint32, decode_motif
from motif_data import Motif, Auxif
from typing import Dict, List, Optional, Tuple, Set, Any
import math
import multiprocessing as mp
import random

# ====================================================================
# NbdGenerator (Neighbor Generation) - Corresponds to C++ NbdGenerator
# ====================================================================

class NbdGenerator:
    """Generates the (l, d) neighborhood for a kmer, saving results as Motif/Auxif pairs."""
    def __init__(self, domain_size: int, curr_array: List[Motif], curr_aux_array: List[Auxif], 
                 x: str, l: int, d: int, leftmost: bool, rightmost: bool):
        self.domain_size = domain_size
        self.curr_array = curr_array
        self.curr_aux_array = curr_aux_array
        self.x = x
        self.l = l
        self.d = d
        self.k = len(x)
        self.leftmost = leftmost
        self.rightmost = rightmost
        self.expanded_count = 0

    def _gen_nbrhood_3(self, alpha: int, count: int):
        current_len = len(self.x)

        if alpha == 0:
            if current_len == self.l:
                motif_encoded_str = "".join(str(c) if c != '-' else '' for c in self.x)
                if len(motif_encoded_str) == self.l:
                    self.curr_array.append(Motif(motif_encoded_str))
                    self.curr_aux_array.append(Auxif(motif_encoded_str, WILDCARD_CODE))
                    self.expanded_count += count
            return

        # Insertion
        for j in range(current_len + 1):
            original_char = self.x[j] if j < current_len else ''

            temp_x = self.x[:j] + str(WILDCARD_CODE) + self.x[j:]
            
            old_x = self.x
            self.x = temp_x
            self._gen_nbrhood_3(alpha - 1, count * self.domain_size) 
            self.x = old_x # Backtrack


    def _gen_nbrhood_2(self, sigma: int, alpha: int, count: int):
        """Handles substitutions."""
        current_len = len(self.x)

        if sigma == 0:
            self._gen_nbrhood_3(alpha, count)
            return

        # Substitution
        for j in range(current_len):
            if self.x[j] == '-': continue
            
            original_char = self.x[j]
            
            old_x = self.x
            self.x = self.x[:j] + str(WILDCARD_CODE) + self.x[j+1:]
            self._gen_nbrhood_2(sigma - 1, alpha, count * self.domain_size)
            self.x = old_x # Backtrack


    def _gen_nbrhood(self, delta: int, sigma: int, alpha: int, count: int):
        """Handles deletions."""
        current_len = len(self.x)
        
        if delta == 0:
            self._gen_nbrhood_2(sigma, alpha, count)
            return

        # Deletion
        for j in range(current_len):
            original_char = self.x[j]
            
            old_x = self.x
            self.x = self.x[:j] + self.x[j+1:]
            self._gen_nbrhood(delta - 1, sigma, alpha, count)
            self.x = old_x[:j] + original_char + old_x[j:] # Backtrack


    def generate(self) -> int:
        """Entry point for neighborhood generation."""
        self.expanded_count = 0
        q = self.k - self.l
        
        for delta in range(max(0, q), math.floor((self.d + q) / 2) + 1):
            alpha = delta - q 
            sigma = self.d - alpha - delta
            
            self._gen_nbrhood(delta, sigma, alpha, 1)
            
        return self.expanded_count

# ====================================================================
# Parallel Workload Functions
# ====================================================================

def _radix_sort_and_intersect(main_array_data: List[Motif], 
                              curr_array: List[Motif], curr_aux_array: List[Auxif], 
                              compact_count: int) -> List[Motif]:
    """Sorts the generated motifs and intersects with the current running result."""
    
    combined_data: List[Motif] = curr_array[:compact_count]
    combined_data.sort(key=lambda m: m.data)

    # Remove duplicates
    sorted_motifs: List[Motif] = []
    last_motif_data = -1

    for motif_obj in combined_data:
        if motif_obj.data != last_motif_data:
            sorted_motifs.append(motif_obj)
            last_motif_data = motif_obj.data
            
    # Perform intersection
    if not main_array_data:
        return sorted_motifs
    else:
        result_motifs: List[Motif] = []
        i = j = 0
        
        while i < len(main_array_data) and j < len(sorted_motifs):
            main_motif = main_array_data[i]
            sorted_motif = sorted_motifs[j]
            
            if main_motif < sorted_motif:
                i += 1
            elif main_motif > sorted_motif:
                j += 1
            else:
                result_motifs.append(main_motif)
                i += 1
                j += 1
                
        return result_motifs

class Worker:
    """Manages the workload partitioning for multiprocessing."""
    def __init__(self, domain_size: int, seq: str, main_array: List[Motif], l: int, d: int):
        self.domain_size = domain_size
        self.seq = seq
        self.main_array = main_array
        self.l = l
        self.d = d
        self.m = len(seq)
        self.works: List[Tuple[int, int]] = []
        
        for k in range(self.l - self.d, self.l + self.d + 1):
            if k > 0 and k <= self.m:
                for i in range(self.m - k + 1):
                    self.works.append((i, k))

        random.seed(42) 
        random.shuffle(self.works)

    def get_load(self) -> int: return len(self.works)

    def process_segment(self, start: int, end: int) -> List[Motif]:
        """Entry point for multiprocessing pool: executes a segment of work."""
        curr_array: List[Motif] = []
        curr_aux_array: List[Auxif] = []
        
        for j in range(start, end):
            i, k = self.works[j]
            x = self.seq[i : i + k]
            
            generator = NbdGenerator(self.domain_size, curr_array, curr_aux_array, 
                                     x, self.l, self.d, (i<=0), (i+k>=self.m))
            generator.generate()
        
        compact_count = len(curr_array)
        
        return _radix_sort_and_intersect(self.main_array, 
                                          curr_array, curr_aux_array, compact_count)

class Ems2p(MotifFinderBase):
    """
    Implements EMS Version 2 Parallel: Uses multiprocessing to distribute 
    neighbor generation, sorting, and intersection (replaces C++ Ems2p).
    """
    def __init__(self, input_path: str, l: int, d: int, params: Params):
        super().__init__("ems2p", input_path, l, d, params)
        self.main_array: List[Motif] = [] 

    @staticmethod
    def _merge_motifs(a: List[Motif], b: List[Motif]) -> List[Motif]:
        """Merges two sorted, unique Motif lists (used in the binary merging tree)."""
        c: List[Motif] = []
        i, j = 0, 0
        
        while i < len(a) or j < len(b):
            if j == len(b) or (i < len(a) and a[i] < b[j]):
                motif_to_add = a[i]
                i += 1
            elif i == len(a) or (j < len(b) and b[j] < a[i]):
                motif_to_add = b[j]
                j += 1
            else: # Equal case
                motif_to_add = a[i]
                i += 1
                j += 1
            
            if not c or motif_to_add != c[-1]:
                c.append(motif_to_add)
                
        return c

    def _gen_all(self, seq: str):
        worker = Worker(self.domain_size, seq, self.main_array, self.l, self.d)
        total_load = worker.get_load()
        num_threads = self.params.num_threads
        ind_load = math.ceil(total_load / num_threads)
        
        tasks: List[Tuple[int, int]] = []
        start = 0
        for _ in range(num_threads):
            pos = min(start + ind_load, total_load)
            if pos > start: tasks.append((start, pos))
            start = pos
            
        with mp.Pool(processes=num_threads) as pool:
            results = pool.starmap(worker.process_segment, tasks)

        current_arrays = [r for r in results if r]
        
        # Binary merging tree
        while len(current_arrays) > 1:
            next_arrays = []
            for i in range(0, len(current_arrays), 2):
                if i + 1 < len(current_arrays):
                    merged = self._merge_motifs(current_arrays[i], current_arrays[i+1])
                    next_arrays.append(merged)
                else:
                    next_arrays.append(current_arrays[i])
            current_arrays = next_arrays

        if current_arrays:
            self.main_array = current_arrays[0] 
        elif not self.main_array:
            self.main_array = []


    def search(self):
        for i, seq in enumerate(self.reads):
            print(f"Processing sequence {i} (Parallel) ...")
            
            self._gen_all(seq)

            if not self.main_array:
                self.motifs.clear()
                print("No common motifs found after intersection. Stopping.")
                return

        # Decode the final list of Motif objects
        self.motifs = [m.get_kmer(self.l) for m in self.main_array]
        self.motifs.sort()
        
        print(f"Num threads = {self.params.num_threads}")