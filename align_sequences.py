##### !/usr/bin/env python (TODO: change back to #!...)

"""
    usage:
        align_sequences [options] seq1.fa seq2.fa
    where the options are:
        -h,--help : print usage and quit
        -m,--match: score of a match in the alignment [2]
        -x,--mismatch: penalty for a mismatch in the alignment [1]
        -g,--gapopen: penalty for opening a new gap [4]
        -e,--gapextend: penalty for extending a gap [1]
"""

from sys import argv, stderr
from getopt import getopt, GetoptError
import numpy as np

def print_matrix(seq1: str, seq2: str, M: np.ndarray) -> None:
    seq1 = "-" + seq1
    seq2 = "-" + seq2
    M_new: np.ndarray = np.zeros((len(seq1) + 2, len(seq2) + 2))
    M_new = M_new.astype(str)
    M_new[0, :] = [" ", " "] + [" " + str(char) + " " for char in range(len(seq2))]
    M_new[:, 0] = [" ", " "] + [str(char) for char in range(len(seq1))]
    M_new[1, 1:] = np.array([" "] + [" " + char + " " for char in seq2]) # spaces align letters with floats (only works for 3 char float (*.*))
    M_new[1:, 1] = np.array([" "] + [char for char in seq1])
    M_new[2:, 2:] = M
    return M_new

def read_single_contig_fasta(filename):
    # a simple function to read the name and sequence from a file
    # The file is expected to have just one contig/sequence. This function
    # checks the assumption and complains if it is not the case.
    names = []
    sequences = []
    with open(filename, 'r') as f:
        line = f.readline()
        assert line.startswith(">")
        names.append(line.strip().split("\t"))
        sequence = ""
        for line in f:
            if line.startswith(">"):
                sequences.append(sequence)
                names.append(line.strip().split("\t"))
                sequence = ""
            else:
                for x in line.strip():
                    if x not in ["A", "C", "G", "T", "a", "c", "g", "t"]: # added lowercase bases
                        print("Unknown nucleotide {}".format(x), file=stderr)
                        exit(3)
                sequence += line.strip()

    sequences.append(sequence)
    assert len(names) == 1
    assert len(sequences) == 1
    return names[0], sequences[0]

def smith_waterman(seq1: str, seq2: str, match: float, mismatch: float, gapopen: float, gapextend: float, print_matrices: bool = False) -> tuple:
    """
    Smith-Waterman Algorithm Implementation. 
    Given strings seq1 and seq2 using the four score parameters, returns an optimal local alignment between seq1 and seq2.
    match should be provided as a REWARD, and thus should be provided as a POSITIVE int/float to be ADDED to the score.
    mismatch, gapopen, and gapextend should be provided as PENALTIES, and thus should be POSITIVE int/float to be SUBTRACTED from the score.  
    """
    # Variables to return
    max_score: int = 0
    alnseq1: str = ""
    alnseq2: str = ""

    # Define lengths of sequences as n and m
    n: int = len(seq1)
    m: int = len(seq2)

    # TODO: Check args
    if m <= 0 or n <= 0:
        raise ValueError("seq1 and seq2 must not be empty strings!")

    # Initialize matrices (gap matrix is redundant with history matrix, but makes some steps easier to implement and more intuitive)
    score_matrix: np.ndarray = np.zeros((n + 1, m + 1)) # num rows, num cols 
    gap_matrix: np.ndarray = np.zeros((n + 1, m + 1)) # elt i, j = 1 if a gap has been initiated or extended at that position. 0 otherwise
    history_matrix: np.ndarray = np.zeros((n + 1, m + 1))
    history_key: dict = {"match" : 1, "substitution" : 2, "insertion_in_seq1" : 3, "insertion_in_seq2" : 4} # this isn't necessary, but is more readable in my opinion

    # Track the index of score_matrix with the maximum score so far
    max_score_coord: np.ndarray = np.array([0, 0])

    # Iterate through matrix row by row, from left to right. 
    # To fill in element (i, j), use row i-1 and element j-1.
    # Keep track of each origin cell at each iteration
    # Don't start in already initialized row 0 and column 0 
    for i in range(1, n + 1): # iterate through rows  
        for j in range(1, m + 1): # iterate through column 
            # Compute the scores for the 4 candidate operations
            match_score: float = (score_matrix[i-1, j-1] + match) * (seq1[i-1] == seq2[j-1]) # ... * 0 when i and j do not match, ... * 1 when they do
            # If the two current chars don't match, then matching should be impossible. Making the match score is one way to solve this issue. I write an expression equal to -1e100 when True, but this is easier
            if seq1[i-1] != seq2[j-1]:
                match_score = -1e100 # really small int that will never be picked
            sub_score: float = score_matrix[i-1, j-1] - mismatch 
            # If there is a gap at (i-1, j) or (i, j-1), then extend it. If there is no gap, then initate it and extend it.
            ins_in_seq1_score: float = score_matrix[i, j-1] - ((bool(gap_matrix[i, j-1]) * gapextend) + ((not bool(gap_matrix[i, j-1])) * (gapextend+ gapopen))) # side arrow
            ins_in_seq2_score: float = score_matrix[i-1, j] - ((bool(gap_matrix[i-1, j]) * gapextend) + ((not bool(gap_matrix[i-1, j])) * (gapextend+ gapopen))) # down arrow
            
            # Get max score based off of only the 0th value in tuple. Tuple enables us to keep track of which option was chosen
            # If there is a tie, the element earlier in the list is chosen. I have ordered these in the list in the order I want the aligner to handle ties.
            score, history = max([(match_score, history_key["match"]), (sub_score, history_key["substitution"]), (ins_in_seq1_score, history_key["insertion_in_seq1"]), (ins_in_seq2_score, history_key["insertion_in_seq2"]), (0, 0)], key = lambda x:x[0]) 

            # Update alignment and history matrix. 
            score_matrix[i, j] = score
            history_matrix[i, j] = history

            # Update gap matrix if a gap was chosen
            if history == history_key["insertion_in_seq1"] or history == history_key["insertion_in_seq2"]:
                gap_matrix[i, j] = 1
        
            # Update matrix maximum
            if score_matrix[i, j] >= max_score: # >= instead of > so that the maximum furthest towards the 'bottom right' of the matrix wins out if there are multiple maxima
                max_score = score_matrix[i, j]
                max_score_coord = (i, j)

    # Backtrack
    alnseq1, alnseq2 = alignment_from_backtrack(history_matrix, history_key, seq1, seq2, max_score_coord)

    if print_matrices:
        print(f"Max score: {max_score} at {max_score_coord}")
        print(f"Score Matrix: \n{print_matrix(seq1, seq2, score_matrix)}\n")
        print(f"History Matrix: \n{print_matrix(seq1, seq2, history_matrix)}\n")
        print(f"Gap Matrix: \n{print_matrix(seq1, seq2, gap_matrix)}")

    return max_score, alnseq1, alnseq2
    

def alignment_from_backtrack(history_matrix: np.ndarray, history_key: dict, seq1: str, seq2: str, max_score_coord: tuple) -> tuple:
    """
    Given: The the history matrix, history key, seq1, seq2, and coordinate in score matrix with maximum value
    Returns: The aligned sequences of seq1 and seq2
    """
    # Retrieve the maximum value in the matrix and TODO: backtrack
    reached_0: bool = False
    curr_backtrack_coord: np.ndarray = max_score_coord
    visited_coords: list = []
    seq1_gap_positons: list = []
    seq2_gap_positons: list = []

    while not reached_0:
        curr_history: float = history_matrix[tuple(curr_backtrack_coord)]
        visited_coords.append(curr_backtrack_coord) # visit coord

        if curr_history == history_key['match'] or curr_history == history_key['substitution']:
            curr_backtrack_coord = curr_backtrack_coord - np.array([1, 1]) # go up and left to i-1, j-1. update current coordinate
        elif curr_history == history_key['insertion_in_seq1']:
            # Gap in seq1
            seq1_gap_positons.append(curr_backtrack_coord[0]) # log gap
            curr_backtrack_coord = curr_backtrack_coord - np.array([0, 1]) # go left one to i, j-1. update current coordinate
        elif curr_history == history_key['insertion_in_seq2']:
            # Gap in seq2
            seq2_gap_positons.append(curr_backtrack_coord[1]) # log gap
            curr_backtrack_coord = curr_backtrack_coord - np.array([1, 0]) # go up one to i-1, j. update current coordinate
        
        # Check if the updated coord's history is 0. If so, stop. 
        if history_matrix[tuple(curr_backtrack_coord)] == 0:
            reached_0 = True
    
    # Construct Alignment
    # Get bounds of sequences that were aligned 
    alnseq1: str = seq1[visited_coords[-1][0] - 1 : visited_coords[0][0]] # -1 to go from matrix indices to string indices. this is because we add a col and a row of 0s to the seqs when initializing matrices. no -1 in the second index because python indices work like [start, stop)
    alnseq2: str = seq2[visited_coords[-1][1] - 1 : visited_coords[0][1]] # -1 to go from matrix indices to string indices. this is because we add a col and a row of 0s to the seqs when initializing matrices
    
    # Insert gaps
    for pos in seq1_gap_positons: # gap coords will be in decreasing order, so we can insert gaps without chaning coords of subsequent gaps
        # Insert after pos - 1(that is the string coord, not matrix coord)
        alnseq1 = alnseq1[: pos] + "-" + alnseq1[pos:] # convert to string coords by subtracting 1
    
    for pos in seq2_gap_positons: # gap coords will be in decreasing order, so we can insert gaps without chaning coords of subsequent gaps
        alnseq2 = alnseq2[: pos] + "-" + alnseq2[pos:] # convert to string coords by subtracting 1

    return alnseq1, alnseq2


def main(filename1, filename2, match, mismatch, gapopen, gapextend):
    # read the name and sequence from the file
    name1, seq1 = read_single_contig_fasta(filename1)
    name2, seq2 = read_single_contig_fasta(filename2)

    # this function takes as input two nucleotide sequences along with
    # scores for an alignment match, mismatch, opening a new gap, and 
    # extending an existing gap. This should return the maximum alignment
    # score as well as the alignment. For examples see the testdriver script
    max_score, alnseq1, alnseq2 = smith_waterman(seq1, seq2, 
                                  match, mismatch, gapopen, gapextend)
    
    print("Maximum alignment score: {}".format(max_score))
    print("Sequence1 : {}".format(alnseq1))
    print("Sequence2 : {}".format(alnseq2))

if __name__ == "__main__":
    try:
        opts, args = getopt(argv[1:],
                     "hm:x:g:e:",
                     ["help", "match=", "mismatch=", "gapopen=", "gapextend="])
    except GetoptError as err:
        print(err)
        print(__doc__, file=stderr)
        exit(1) 

    match = 2
    mismatch = 1
    gapopen = 4
    gapextend = 1

    for o, a in opts:
        if o in ("-h", "--help"):
            print(__doc__, file=stderr)
            exit()
        elif o in ("-m", "--match"):
            match = float(a)
        elif o in ("-x", "--mismatch"):
            mismatch = float(a)
        elif o in ("-g", "--gapopen"):
            gapopen = float(a)
        elif o in ("-e", "--gapextend"):
            gapextend = float(a)
        else:
            assert False, "unhandled option"

    if len(args) != 2:
        print(__doc__, file=stderr)
        exit(2)

    main(args[0], args[1], match, mismatch, gapopen, gapextend)