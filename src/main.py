import argparse
import sys
import os
import multiprocessing
from utils import Params
from motif_finder import MotifFinderBase
from ems1 import Ems1
from ems2 import Ems2
from ems2p import Ems2p
from typing import Optional

def usage(argv: str):
    print(f"Usage: {os.path.basename(argv)} [OPTIONS] <input-sequence-file>")
    print("\t-s <version>    Possible versions: 1, 2, 2m, 2p")
    print("\t-l <l>          Length (l) of (l,d) motif")
    print("\t-d <d>          Maximum edit distance (d) of (l,d) motif")
    print("\t-t <int>        Number of processes (threads) for Ems2p (default is CPU count)")
    sys.exit(-1)

def main():
    # Set the start method for multiprocessing compatibility
    if sys.platform.startswith('win'):
         multiprocessing.freeze_support()
    elif sys.platform.startswith('darwin') or sys.platform.startswith('linux'):
         # Fork is default but set_start_method might be needed on macOS/Linux if issues arise
         try:
              multiprocessing.set_start_method('fork', force=True)
         except Exception:
              pass
              
    parser = argparse.ArgumentParser(add_help=False)
    
    parser.add_argument('input', nargs='?', help="Input sequence file")
    parser.add_argument('-s', dest='version', default='2', choices=['1', '2', '2m', '2p'], 
                        help="Algorithm version: 1 (Brute Force), 2 (Fast Trie), 2m (Simple Trie), 2p (Parallel)")
    parser.add_argument('-l', dest='l', type=int, default=1, help="Motif length (l)")
    parser.add_argument('-d', dest='d', type=int, default=1, help="Max edit distance (d)")
    parser.add_argument('-t', dest='num_threads', type=int, 
                        default=multiprocessing.cpu_count(), 
                        help="Number of threads/processes for Ems2p")
    parser.add_argument('-h', '--help', action='store_true', help="Show help message and exit")

    args = parser.parse_args()

    if args.help or args.input is None:
        usage(sys.argv[0])
        return

    params = Params(l=args.l, d=args.d, num_threads=args.num_threads)
    
    if params.l <= 0 or params.d < 0:
        print("Error: l must be > 0 and d must be >= 0.", file=sys.stderr)
        sys.exit(-1)

    ems_instance: Optional[MotifFinderBase] = None
    version = args.version
    input_path = args.input

    try:
        if version == "1":
            ems_instance = Ems1(input_path, params.l, params.d, params)
        elif version == "2":
            ems_instance = Ems2(input_path, params.l, params.d, params, "ems2")
        elif version == "2m":
            ems_instance = Ems2(input_path, params.l, params.d, params, "ems2m")
        elif version == "2p":
            ems_instance = Ems2p(input_path, params.l, params.d, params)
        
        if ems_instance:
            ems_instance.search_write_motifs()

    except Exception as e:
        print(f"An error occurred during execution: {e}", file=sys.stderr)
        sys.exit(-1)


if __name__ == "__main__":
    main()