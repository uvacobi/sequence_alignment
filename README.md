# Aligning sequences

## Implement the Smith-Waterman algorithm to align two sequences.

For this assignment, implement the `smith_waterman` function in align_sequences.py. The function takes in two nucleotide sequences (seq1 and seq2), and scores/penalty for a match, mismatch, opening of a new gap, and extension of a gap. You should write the `smith_waterman` function to calculate the maximum score of the local alignment and the resulting aligned subsequences from the two input sequences. For example if seq1 is `ACGTTA` and seq2 is `ACTTA`, match=1, mismatch=1, gapopen=1, gapextend=0, then the expected output of the function will be `4, 'ACGTTA','AC-TTA'`. As in the example, all input scores for match, mismatch, gapopen and gapextend are expected to be $\geq 0$ and your code is expected to decide when to add or subtract those values.

Please implement the following gap-penalty $S = O + k*E$, where $S$ is the score for the gap, $O$ is the gap-open penalty, $k$ is the size of the gap and $E$ is the gap-extend penalty. 

## Testing

You can test your work by running:

```
python3 testdriver
```
