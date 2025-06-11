# =========================================================================
# This software is licensed Apache License, Version 2.0.
# Python implementation of the Smith-Waterman Algorithm.

# This version is implemented by lopozz 2025-03-01 <https://github.com/lopozz/smith-waterman>
# Implementation reference https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm#Algorithm
# =========================================================================


class SWParams:
    def __init__(self, similarity_score=10, gap_penalty=0, mismatch_penalty=1):
        self.similarity_score = similarity_score
        self.gap_penalty = gap_penalty
        self.mismatch_penalty = mismatch_penalty


def zeros(shape):
    m = []
    for _ in range(shape[0]):
        m.append([])
        for _ in range(shape[1]):
            m[-1].append(0)
    return m


def match_score(a, b, params=None):
    if params is None:
        params = SWParams()

    return (
        params.similarity_score
        if a == b
        else params.gap_penalty
        if a == "-" or b == "-"
        else params.mismatch_penalty
    )


def identity_score(align_s1, align_s2):
    mathches = sum([1 for i, _ in enumerate(align_s1) if align_s1[i] == align_s2[i]])
    identity = float(mathches) / len(align_s1) * 100
    return identity


def water(s1, s2, params=None):
    if params is None:
        params = SWParams()

    # 1. Initialize a scoring matrix H
    m, n = len(s1), len(s2)
    H = zeros((m + 1, n + 1))

    # Initialize traceback path pointer matrix
    P = zeros((m + 1, n + 1))

    # 2. Fill the scoring matrix
    max_score = 0
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = H[i - 1][j - 1] + match_score(s1[i - 1], s2[j - 1])
            insert = H[i][j - 1] + params.gap_penalty
            delete = H[i - 1][j] + params.gap_penalty
            zero = 0
            H[i][j] = max(match, delete, insert, zero)

            # Determine pointer direction
            if H[i][j] == zero:
                P[i][j] = 0  # End of the path
            if H[i][j] == delete:
                P[i][j] = 1  # Trace up
            if H[i][j] == insert:
                P[i][j] = 2  # Trace left
            if H[i][j] == match:
                P[i][j] = 3  # Trace diagonal

            # Update max score position
            if H[i][j] >= max_score:
                max_i, max_j, max_score = i, j, H[i][j]

    # 3. Traceback
    align_s1, align_s2 = [], []
    i, j = max_i, max_j

    while True:
        pointer = P[i][j]
        if pointer == 0:  # End of the path
            break
        if pointer == 3:
            align_s1.append(s1[i - 1])
            align_s2.append(s2[j - 1])
            i -= 1
            j -= 1
        elif pointer == 2:
            align_s1.append("-")
            align_s2.append(s2[j - 1])
            j -= 1
        elif pointer == 1:
            align_s1.append(s1[i - 1])
            align_s2.append("-")
            i -= 1

    align_s1.reverse()
    align_s2.reverse()

    return align_s1, align_s2
