import random
import time
import sys
import os

N = 20  # Number of sequences
M = 600 # Length of each sequence
L = 11  # Motif length
D = 3   # Max edit distance
DNA_ALPHABET = "ACGT"

# --- Global Storage ---
_seqs = []
_motif = ['a'] * L
_edited_motif = [''] * N
_pos = [0] * N

def generate_seqs():
    """Generates N random sequences of length M."""
    global _seqs
    _seqs.clear()
    for _ in range(N):
        temp = [random.choice(DNA_ALPHABET) for _ in range(M)]
        _seqs.append("".join(temp))

def generate_motif():
    """Generates a random consensus motif of length L."""
    global _motif
    _motif = [random.choice(DNA_ALPHABET) for _ in range(L)]
    _motif = "".join(_motif)

def edit_motif(motif: str) -> str:
    """
    Applies deletions (delta), insertions (alpha), and substitutions (beta) 
    to the motif up to a total distance of D.
    """
    temp_motif = list(motif)
    
    # 1. Deletions (delta)
    delta = random.randint(0, D)
    for _ in range(delta):
        if not temp_motif: break
        pos = random.randrange(len(temp_motif))
        temp_motif.pop(pos)

    # 2. Insertions (alpha)
    alpha = random.randint(0, D - delta)
    for _ in range(alpha):
        # Insertions can happen at position 0 up to len(temp_motif)
        pos = random.randrange(len(temp_motif) + 1)
        x = random.choice(DNA_ALPHABET)
        temp_motif.insert(pos, x)

    # 3. Substitutions (beta)
    beta = D - delta - alpha
    for _ in range(beta):
        if not temp_motif: continue
        pos = random.randrange(len(temp_motif))
        original_char = temp_motif[pos]
        
        # Ensure the substitution is a change
        x = original_char
        while x == original_char:
            x = random.choice(DNA_ALPHABET)
            
        temp_motif[pos] = x
        
    return "".join(temp_motif)

def planted_seqs():
    """Edits the motif and plants the edited version into each sequence."""
    global _seqs, _edited_motif, _pos
    
    for i in range(N):
        # 1. Edit the motif (ensures l-d <= length <= l+d)
        edited = edit_motif(_motif)
        _edited_motif[i] = edited
        
        # 2. Determine planting position
        end_pos = M - len(edited)
        if end_pos < 0:
            # Should not happen with M=600, L=11, D=3, but safe check
            _pos[i] = 0
            continue
            
        _pos[i] = random.randrange(end_pos + 1)
        
        # 3. Plant the edited motif
        seq_list = list(_seqs[i])
        seq_list[_pos[i]:_pos[i] + len(edited)] = list(edited)
        _seqs[i] = "".join(seq_list)


def main():
    # Use current time as seed (C++ equivalent: srand(time(NULL)))
    random.seed(int(time.time()))
    
    generate_seqs()
    generate_motif()
    planted_seqs()

    file_name = f"planted_l{L}_d{D}.txt"
    
    # Ensure the test folder exists (in case the user runs this outside the intended structure)
    os.makedirs(os.path.dirname(file_name) or '.', exist_ok=True)
    
    try:
        with open(file_name, 'w') as myFile:
            for i in range(N):
                myFile.write(f"{i} Motif {_motif} planted as {_edited_motif[i]} at position {_pos[i]}\n")
                myFile.write(f"{_seqs[i]}\n")
        print(f"Successfully generated test case: {file_name}")
        print(f"Planted consensus motif was: {_motif}")
        
    except IOError:
        print(f"Unable to open file {file_name}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()