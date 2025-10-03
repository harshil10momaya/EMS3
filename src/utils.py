import argparse
import sys
import time
import os
import math
from typing import List, Tuple, Dict, Any, Union

# --- Windows Compatibility Fix for Memory Tracking ---
try:
    import resource
    def getrusage_maxrss():
        return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
except ImportError:
    try:
        import psutil
        def getrusage_maxrss():
            # Convert bytes to KB
            return psutil.Process(os.getpid()).memory_info().rss // 1024
    except ImportError:
        # Placeholder if neither resource nor psutil is available
        def getrusage_maxrss():
            return 0
# --- End of Compatibility Fix ---

# ====================================================================
# TYPE DEFINITIONS and CONSTANTS
# ====================================================================

Reads = List[str] 
Motifs = List[str]

DNA_DOMAIN = "ACGT"
ENCODING_MAP = {char: i for i, char in enumerate(DNA_DOMAIN)}
DECODING_MAP = {i: char for i, char in enumerate(DNA_DOMAIN)}
WILDCARD_CODE = len(DNA_DOMAIN)

uchar = int
uint32 = int
uint64 = int

# ====================================================================
# FILE I/O AND SEQUENCE PROCESSING
# ====================================================================

def get_alphabet(raw_strings: List[str]) -> str:
    return DNA_DOMAIN

def encode_strings(raw_reads: List[str], domain: str) -> List[str]:
    encoding_map = {char: str(i) for i, char in enumerate(domain)}
    encoded_reads = []
    for read_seq in raw_reads:
        encoded_seq = "".join(encoding_map.get(char, '0') for char in read_seq.upper())
        encoded_reads.append(encoded_seq)
    return encoded_reads

def read_file(filepath: str, reads: Reads) -> Tuple[str, List[str]]:
    raw_reads = []
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip().upper().replace('U', 'T')
                if line and not line.startswith('>'): 
                    raw_reads.append(line)
        
        domain = get_alphabet(raw_reads) 
        encoded_sequences = encode_strings(raw_reads, domain)
        
        reads.clear()
        reads.extend(encoded_sequences)
        print(f"Read and encoded {len(reads)} sequences.")
        return domain, encoded_sequences
        
    except FileNotFoundError:
        print(f"ERROR: could not open file {filepath}", file=sys.stderr)
        sys.exit(-1)

def decode_motif(encoded_motif: str, domain: str) -> str:
    decoding_map = {i: char for i, char in enumerate(domain)}
    decoded_seq = ""
    for char_code_str in encoded_motif:
        code = int(char_code_str) if char_code_str.isdigit() else WILDCARD_CODE
        if code == WILDCARD_CODE:
             decoded_seq += '*'
        else:
             decoded_seq += decoding_map.get(code, '?')
    return decoded_seq

def remove_extension(filename: str) -> str:
    return os.path.splitext(filename)[0]

def get_out_file(input_path: str, l: int, d: int, name: str) -> str:
    """Generates the output file name and ensures the 'output' directory exists."""
    # 1. Determine the path of the input file
    input_dir = os.path.dirname(input_path)
    
    # 2. Define the new output directory (e.g., 'test/output')
    output_dir = os.path.join(input_dir, 'output')
    
    # 3. Create the directory if it does not exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 4. Generate the base file name
    base_name = os.path.basename(input_path)
    base_no_ext = remove_extension(base_name)
    filename = f"out_{base_no_ext}_{name}_l{l}_d{d}.txt"
    
    # 5. Return the full path to the new file inside the output directory
    return os.path.join(output_dir, filename)

def diffclock(start_time: float) -> float:
    return time.time() - start_time

def edist(s1: str, s2: str) -> int:
    """Calculates the Edit Distance (Levenshtein distance)."""
    len1, len2 = len(s1), len(s2)
    d = [[0] * (len2 + 1) for _ in range(len1 + 1)]

    for i in range(len1 + 1): d[i][0] = i
    for j in range(len2 + 1): d[0][j] = j

    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            cost = 0 if s1[i - 1] == s2[j - 1] else 1
            d[i][j] = min(d[i - 1][j] + 1, d[i][j - 1] + 1, d[i - 1][j - 1] + cost)
            
    return d[len1][len2]

def found_in_seq(candidate: str, seq: str, l: int, d: int) -> bool:
    m = len(seq)
    len_candi = len(candidate)
    
    for q in range(-d, d + 1):
        k = len_candi + q 
        
        if k <= 0 or k > m: continue
            
        for i in range(m - k + 1):
            x = seq[i : i + k]
            if edist(x, candidate) <= d:
                return True
                
    return False

# ====================================================================
# PARAMS structure
# ====================================================================

class Params:
    def __init__(self, l: int = 1, d: int = 1, num_threads: int = 1):
        self.l = l
        self.d = d
        self.num_threads = num_threads