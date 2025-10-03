from typing import Any, List
from utils import WILDCARD_CODE

class Motif:
    """
    Python class mimicking the C++ bit-packed Motif structure (2 bits per char).
    Uses a single large integer (`self.data`) for bit-packed storage.
    Corresponds to C++ `motif.hpp`.
    """
    MAX_L = 32 

    def __init__(self, x: str = None):
        self.data: int = 0
        if x is not None:
            self.set_from_kmer(x)

    def set(self, other: 'Motif'):
        self.data = other.data
        
    def set_from_kmer(self, x: str):
        self.data = 0
        length = len(x)
        for i in range(length):
            char_code = int(x[i]) & 0b11
            j = length - 1 - i
            p = 2 * j
            self.data |= char_code << p

    def get_kmer(self, k: int) -> str:
        """Unpacks the data back to an ENCODED motif string of length k (e.g., "0123")."""
        if k > self.MAX_L:
            raise ValueError(f"Requested length {k} exceeds max Motif length {self.MAX_L}")
        
        unpacked_codes: List[str] = [''] * k
        for i in range(k):
            p = 2 * (k - 1 - i)
            code = (self.data >> p) & 0b11
            unpacked_codes[i] = str(code)
            
        return "".join(unpacked_codes)

    def get_2bits(self, p: int) -> int:
        """Gets the 2-bit code at bit position p."""
        return (self.data >> p) & 0b11

    def clear(self):
        self.data = 0
        
    # Comparison operators (essential for sorting and intersection)
    def __lt__(self, other: 'Motif') -> bool: return self.data < other.data
    def __eq__(self, other: Any) -> bool: return isinstance(other, Motif) and self.data == other.data
    def __hash__(self) -> int: return hash(self.data)

class Auxif:
    """
    Python class mimicking the C++ Auxif structure (1 bit per char, 32-bit).
    Used in Ems2p for tracking wildcards/differences during sorting.
    """
    def __init__(self, x: str = None, wildcard_code: int = WILDCARD_CODE):
        self.data: int = 0 
        if x is not None:
            self.set_from_kmer(x, wildcard_code)

    def set(self, other: 'Auxif'):
        self.data = other.data
        
    def set_from_kmer(self, x: str, wildcard_code: int):
        self.data = 0
        length = len(x)
        for i in range(length):
            char_code = int(x[i]) if x[i].isdigit() else int(x[i])
            j = length - 1 - i
            p = j
            
            if char_code == wildcard_code:
                 self.data |= 1 << p

    def get_bit(self, p: int) -> int:
        """Gets the 1-bit code at bit position p."""
        return (self.data >> p) & 1
    
    def __lt__(self, other: 'Auxif') -> bool: return self.data < other.data