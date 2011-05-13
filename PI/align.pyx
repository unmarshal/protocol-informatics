cdef extern from "Numeric/arrayobject.h":

    struct PyArray_Descr:
        int type_num, elsize
        char type

    ctypedef class Numeric.ArrayType [object PyArrayObject]:
        cdef char *data
        cdef int nd
        cdef int *dimensions, *strides
        cdef object base
        cdef PyArray_Descr *descr
        cdef int flags

import Numeric

def NeedlemanWunsch(seq1, seq2, S, int g, int e):

    cdef int M, N, i, j, t_max, i_max, j_max, dir
    cdef int v1, v2, v3, m, new_len, data, ni, nj

    cdef ArrayType a
    cdef int nrows, ncols
    cdef int *matrix, x

    edits1 = []
    edits2 = []

    M = len(seq1) + 1
    N = len(seq2) + 1

    table = Numeric.zeros([M, N], 'i')

    a = table

    nrows = a.dimensions[0]
    ncols = a.dimensions[1]

    matrix = <int *>a.data

    # Iterate through matrix and score similarities
    for i from 1 <= i < M:
        for j from 1 <= j < N:
            if seq1[i - 1] == seq2[j - 1]:
                matrix[i * ncols + j] = 1
            else:
                matrix[i * ncols + j] = 0

    # Sum the matrix
    i_max = 0
    j_max = 0
    t_max = 0

    for i from 1 <= i < M:
        for j from 1 <= j < N:

            dir = 0

            v1 = matrix[(i - 1) * ncols + (j - 1)]
            v2 = matrix[i * ncols + (j - 1)]
            v3 = matrix[(i - 1) * ncols + j]

            if v1 > 255:
                v1 = v1 >> 8
            if v2 > 255:
                v2 = v2 >> 8
            if v3 > 255:
                v3 = v3 >> 8

            v1 = v1 + matrix[i * ncols + j]
            v2 = v2 - e
            v3 = v3 - e

            if v1 > 0:
                m = v1
            else:
                m = 0

            if v2 > m:
                m = v2

            if v3 > m:
                m = v3

            if m == v1:
                dir = dir | (1 << 0)

            if m == v2:
                dir = dir | (1 << 1)

            if m == v3:
                dir = dir | (1 << 2)

            if m >= t_max:
                t_max = m
                i_max = i
                j_max = j

            m = m << 8
            m = m | dir

            matrix[i * ncols + j] = m

    # Do backtrace through matrix
    i = i_max
    j = j_max

    new_len = 0

    while i and j:
        data = matrix[i * ncols + j]

        if data & (1 << 2):
            ni = i - 1
            nj = j

        if data & (1 << 1):
            ni = i
            nj = j - 1

        if data & (1 << 0):
            ni = i - 1
            nj = j - 1

        new_len = new_len + 1

        i = ni
        j = nj

    new_seq1 = Numeric.zeros((new_len), 'i')
    new_seq2 = Numeric.zeros((new_len), 'i')

    cdef int s1, s2, gaps
    s1 = s2 = new_len
    gaps = 0

    i = i_max
    j = j_max

    while i and j:
        data = matrix[i * ncols + j]

        if data & (1 << 2):
            ni = i - 1
            nj = j
            dir = (1 << 2)

        if data & (1 << 1):
            ni = i
            nj = j - 1
            dir = (1 << 1)

        if data & (1 << 0):
            ni = i - 1
            nj = j - 1
            dir = (1 << 0)

        if dir == (1 << 0):
            s1 = s1 - 1
            s2 = s2 - 1
            new_seq1[s1] = seq1[i - 1]
            new_seq2[s2] = seq2[j - 1]

        if dir == (1 << 1):
            s1 = s1 - 1
            s2 = s2 - 1

            edits1.append(s1)
            new_seq1[s1] = 256 # '_'
            new_seq2[s2] = seq2[j - 1]
            gaps = gaps + 1

        if dir == (1 << 2):
            s1 = s1 - 1
            s2 = s2 - 1

            new_seq1[s1] = seq1[i - 1]
            new_seq2[s2] = 256 # '_'
            edits2.append(s2)
            gaps = gaps + 1

        i = ni
        j = nj

    return (new_seq1, new_seq2, edits1, edits2, t_max, gaps)

def SmithWaterman(seq1, seq2, S, int g, int e):

    cdef int M, N, i, j, t_max, i_max, j_max, dir
    cdef int v1, v2, v3, m, new_len, data, ni, nj

    cdef ArrayType a
    cdef int nrows, ncols
    cdef int *matrix, x

    edits1 = []
    edits2 = []

    M = len(seq1) + 1
    N = len(seq2) + 1

    table = Numeric.zeros([M, N], 'i')

    for i in range(M):
        table[i][0] = 0 - i

    for i in range(N):
        table[0][i] = 0 - i

    a = table

    nrows = a.dimensions[0]
    ncols = a.dimensions[1]

    matrix = <int *>a.data

    # Iterate through matrix and score similarities
    for i from 1 <= i < M:
        for j from 1 <= j < N:
            if seq1[i - 1] == seq2[j - 1]:
                matrix[i * ncols + j] = 2
            else:
                matrix[i * ncols + j] = -1

    # Sum the matrix
    i_max = 0
    j_max = 0
    t_max = 0

    for i from 1 <= i < M:
        for j from 1 <= j < N:

            dir = 0

            v1 = matrix[(i - 1) * ncols + (j - 1)]
            v2 = matrix[i * ncols + (j - 1)]
            v3 = matrix[(i - 1) * ncols + j]

            if v1 > 255 or (v1 & 0xffffff00) == False:
                v1 = v1 >> 8
            if v2 > 255 or (v1 & 0xffffff00) == False:
                v2 = v2 >> 8
            if v3 > 255 or (v1 & 0xffffff00) == False:
                v3 = v3 >> 8

            v1 = v1 + matrix[i * ncols + j]
            v2 = v2 - 2
            v3 = v3 - 2

            if v1 > 0:
                m = v1
            else:
                m = 0

            if v2 > m:
                m = v2

            if v3 > m:
                m = v3

            if m == v1:
                dir = dir | (1 << 0)

            if m == v2:
                dir = dir | (1 << 1)

            if m == v3:
                dir = dir | (1 << 2)

            if m >= t_max:
                t_max = m
                i_max = i
                j_max = j

            m = m << 8
            m = m | dir

            matrix[i * ncols + j] = m

    # Do backtrace through matrix
    i = i_max
    j = j_max

    new_len = 0

    while matrix[i * ncols + j] > 0:
        data = matrix[i * ncols + j]

        if data & (1 << 2):
            ni = i - 1
            nj = j

        if data & (1 << 1):
            ni = i
            nj = j - 1

        if data & (1 << 0):
            ni = i - 1
            nj = j - 1

        new_len = new_len + 1

        i = ni
        j = nj

    new_seq1 = Numeric.zeros((new_len), 'i')
    new_seq2 = Numeric.zeros((new_len), 'i')

    cdef int s1, s2, gaps
    s1 = s2 = new_len
    gaps = 0

    i = i_max
    j = j_max

    while matrix[i * ncols + j] > 0:
        data = matrix[i * ncols + j]

        if data & (1 << 2):
            ni = i - 1
            nj = j
            dir = (1 << 2)

        if data & (1 << 1):
            ni = i
            nj = j - 1
            dir = (1 << 1)

        if data & (1 << 0):
            ni = i - 1
            nj = j - 1
            dir = (1 << 0)

        if dir == (1 << 0):
            s1 = s1 - 1
            s2 = s2 - 1
            new_seq1[s1] = seq1[i - 1]
            new_seq2[s2] = seq2[j - 1]

        if dir == (1 << 1):
            s1 = s1 - 1
            s2 = s2 - 1

            edits1.append(s1)
            new_seq1[s1] = 256 # '_'
            new_seq2[s2] = seq2[j - 1]
            gaps = gaps + 1

        if dir == (1 << 2):
            s1 = s1 - 1
            s2 = s2 - 1

            new_seq1[s1] = seq1[i - 1]
            new_seq2[s2] = 256 # '_'
            edits2.append(s2)
            gaps = gaps + 1

        i = ni
        j = nj

    return (new_seq1, new_seq2, edits1, edits2, t_max, gaps)
