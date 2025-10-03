from abc import ABC, abstractmethod
import time
from typing import List, Tuple
from utils import Params, read_file, get_out_file, decode_motif, getrusage_maxrss, diffclock, edist
import sys
import os
import re

class MotifFinderBase(ABC):
    """
    Base class for motif finding algorithms.
    """
    def __init__(self, name: str, input_path: str, l: int, d: int, params: Params):
        self.name: str = name
        self.input: str = input_path
        self.l: int = l
        self.d: int = d
        self.params: Params = params
        
        self.reads: List[str] = []
        self.domain: str = ""
        self.motifs: List[str] = [] # Stores decoded motifs (results)
        self.consensus_motif: str = self._extract_consensus_motif(input_path) 

        # Load file (reads are encoded here)
        self.domain, _ = read_file(input_path, self.reads)
        self.domain_size: int = len(self.domain)
        
    def _extract_consensus_motif(self, input_path: str) -> str:
        """
        Attempts to read the consensus motif from the first line of the input file.
        """
        try:
            with open(input_path, 'r') as f:
                first_line = f.readline().strip()
                # Regex to find "Motif [SEQUENCE]"
                match = re.search(r'Motif\s+([ACGTacgt]+)\s+planted', first_line)
                if match:
                    return match.group(1).upper()
                return ""
        except Exception:
            return "" 
        
    @abstractmethod
    def search(self):
        """Must be implemented by derived classes to perform the core search."""
        pass

    def search_write_motifs(self):
        """Performs search and writes results, handling timing and logging."""
        
        output_path = get_out_file(self.input, self.l, self.d, self.name)
        
        # Determine log file path (in the same folder as the output file)
        log_dir = os.path.dirname(output_path)
        log_path = os.path.join(log_dir, "emsTimeMemory.log")
        
        print(f"l      = {self.l}, d = {self.d}")
        print(f"input  = {self.input}")
        print(f"output = {output_path}")
        if self.consensus_motif:
            print(f"Consensus: {self.consensus_motif} (Length {len(self.consensus_motif)})")

        begin_time = time.time()
        self.search()
        elapsed = diffclock(begin_time)
        max_rss = getrusage_maxrss()
        
        # --- Prepare Output with Distance ---
        final_output_lines = []

        for motif in self.motifs:
            decoded_motif = decode_motif(motif, self.domain) if any(c.isdigit() for c in motif) else motif
            
            # Calculate distance only if lengths match (required for edit distance comparison)
            distance = "N/A"
            if self.consensus_motif and len(decoded_motif) == len(self.consensus_motif):
                distance = edist(decoded_motif, self.consensus_motif)
                
            final_output_lines.append(f"{decoded_motif}\tDistance: {distance}")
        
        # Write performance to emsTimeMemory.log
        with open(log_path, 'a') as f:
            f.write(f"{self.name}: ({self.l},{self.d}) Edited Motifs found using {self.params.num_threads} threads:(in {elapsed:.4f} sec, using {max_rss} KB): {len(final_output_lines)}\n")
            
        # Write motifs and distances to the output file
        with open(output_path, 'w') as out:
            for line in final_output_lines:
                out.write(f"{line}\n")
        
        # Console output
        print(f"\r{self.name}: ({self.l},{self.d}) Edited Motifs found (in {elapsed:.4f} sec, using {max_rss} KB): {len(final_output_lines)}        ")
