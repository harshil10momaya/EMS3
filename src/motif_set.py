import heapq
from typing import List, Tuple, Any
from motif_data import Motif

# ====================================================================
# MOTIFSET (K-way Merge using Heap)
# Replaces custom C++ MotifSet min-heap implementation.
# ====================================================================

class MotifSet:
    """
    Implements a min-heap (priority queue) for K-way merging of sorted Motif lists.
    """
    def __init__(self):
        # Heap stores tuples: (Motif_data, List_Index, Element_Index_in_List)
        self.heap: List[Tuple[int, int, int]] = [] 
        # Stores the original lists: list of [Motif_array, current_pos, end_pos]
        self.data_desc: List[List[Any]] = []
        self.desc_pos: int = 0
        self.last_popped_motif: Motif = Motif()

    def init_add(self, buffer: List[Motif], start_pos: int, end_pos: int):
        """Adds a new sorted list (buffer) segment."""
        if start_pos < end_pos:
            desc_id = self.desc_pos
            self.data_desc.append([buffer, start_pos, end_pos])
            
            first_motif = buffer[start_pos]
            heapq.heappush(self.heap, (first_motif.data, desc_id, start_pos))
            
            self.desc_pos += 1

    def clear(self):
        """Clears the heap and descriptors."""
        self.heap.clear()
        self.data_desc.clear()
        self.desc_pos = 0
        self.last_popped_motif.clear()

    def get_min(self, motif: Motif) -> bool:
        """Extracts the unique minimum motif, advances the source list."""
        if not self.heap:
            return False

        while self.heap:
            min_data, desc_id, element_idx = heapq.heappop(self.heap)

            buffer, current_pos_ref, end_pos = self.data_desc[desc_id]
            extracted_motif = buffer[element_idx]

            # Advance the position in the source list
            self.data_desc[desc_id][1] += 1
            current_pos = self.data_desc[desc_id][1]

            # Push the next element if list is not exhausted
            if current_pos < end_pos:
                next_motif = buffer[current_pos]
                heapq.heappush(self.heap, (next_motif.data, desc_id, current_pos))
            
            # Check for uniqueness against the last returned motif
            if extracted_motif != self.last_popped_motif:
                 motif.set(extracted_motif)
                 self.last_popped_motif.set(extracted_motif)
                 return True
                 
        return False