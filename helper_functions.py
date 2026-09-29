from Bio.Align import substitution_matrices
sub_matrix = substitution_matrices.load("BLOSUM62")


def global_alignment(seq1, seq2, scoring_function):
    """Global sequence alignment using the Needleman–Wunsch algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> global_alignment("abracadabra", "dabarakadara", lambda x, y: [-1, 1][x == y])
    ('-ab-racadabra', 'dabarakada-ra', 5.0)

    Other alignments are not possible.

    """

    n = len(seq1)
    m = len(seq2)

    alignments = [[0] * (m + 1) for _ in range(n + 1)]
    pred = [[None] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        alignments[i][0] = alignments[i-1][0] + scoring_function(seq1[i-1], "-")
        pred[i][0] = (i-1, 0)

    for j in range(1, m + 1):
        alignments[0][j] = alignments[0][j-1] + scoring_function("-", seq2[j-1])
        pred[0][j] = (0, j-1)

    for i in range(1, n + 1):
        for j in range(1, m + 1):

            match = alignments[i-1][j-1] + scoring_function(seq1[i-1], seq2[j-1])
            seq1_ins = alignments[i-1][j] + scoring_function("-", seq2[j-1])
            seq2_ins = alignments[i][j-1] + scoring_function(seq1[i-1], "-")

            best_score = max(match, seq1_ins, seq2_ins)
            alignments[i][j] = best_score

            if best_score == match:
                pred[i][j] = (i-1, j-1)
            elif best_score == seq1_ins:
                pred[i][j] = (i-1, j)
            else:
                pred[i][j] = (i, j-1)

    i = n
    j = m
    same = 0
    total = 0

    seq1_final = []
    seq2_final = []

    while i > 0 or j > 0:
        pred_i, pred_j = pred[i][j]

        if pred_i == i:
            seq1_final.append("-")
            seq2_final.append(seq2[j-1])

        elif pred_j == j:
            seq1_final.append(seq1[i-1])
            seq2_final.append("-")

        else:
            seq1_final.append(seq1[i-1])
            seq2_final.append(seq2[j-1])
        total += 1

        if seq1[i-1] == seq2[j-1]:
            same += 1
        i, j = pred_i, pred_j

    # Since traversal was from the last index, we have to reverse the strings
    seq1_final.reverse()
    seq2_final.reverse()

    return (
        "".join(seq1_final),
        "".join(seq2_final),
        float(alignments[n][m]),
        float(100 * (same / total))
    )


def local_alignment(seq1, seq2, scoring_function):
    """Local sequence alignment using the Smith-Waterman algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> local_alignment("pending itch", "unending glitch", lambda x, y: [-1, 1][x == y])
    ('ending --itch', 'ending glitch', 9.0)

    Other alignments are not possible.

    """

    n = len(seq1)
    m = len(seq2)

    alignments = [[0] * (m + 1) for _ in range(n + 1)]
    pred = [[None] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        alignments[i][0] = 0

    for j in range(1, m + 1):
        alignments[0][j] = 0

    for i in range(1, n + 1):
        for j in range(1, m + 1):

            match = alignments[i-1][j-1] + scoring_function(seq1[i-1], seq2[j-1])
            seq1_ins = alignments[i-1][j] + scoring_function("-", seq2[j-1])
            seq2_ins = alignments[i][j-1] + scoring_function(seq1[i-1], "-")

            best_score = max(0, match, seq1_ins, seq2_ins)
            alignments[i][j] = best_score

            if best_score < 0:
                continue
            elif best_score == match:
                pred[i][j] = (i-1, j-1)
            elif best_score == seq1_ins:
                pred[i][j] = (i-1, j)
            else:
                pred[i][j] = (i, j-1)
    
    maximum = max(max(row) for row in alignments)

    i = None
    j = None

    for x, row in enumerate(alignments):
        for y, value in enumerate(row):
            if value == maximum:
                i = x
                j = y
                break

    max_score = float(alignments[i][j])

    seq1_final = []
    seq2_final = []

    while (i > 0 or j > 0) and alignments[i][j] > 0:
        if pred[i][j] == None:
            break
            
        pred_i, pred_j = pred[i][j]

        if pred_i == i:
            seq1_final.append("-")
            seq2_final.append(seq2[j-1])

        elif pred_j == j:
            seq1_final.append(seq1[i-1])
            seq2_final.append("-")

        else:
            seq1_final.append(seq1[i-1])
            seq2_final.append(seq2[j-1])

        i, j = pred_i, pred_j

    # Since traversal was from the last index, we have to reverse the strings
    seq1_final.reverse()
    seq2_final.reverse()

    return (
        "".join(seq1_final),
        "".join(seq2_final),
        max_score
    )


## This is an example scoring function, you should implement a version which uses a scoring matrix 
def scoring_function(x, y):
    if x == "-" or y == "-":
        # TODO check scoring gap
        return -5

    return sub_matrix[x, y]
