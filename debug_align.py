import numpy as np

MARRAY = 1
IARRAY = 2  
JARRAY = 3

def make_matrix(match_score=5, mismatch_score=-4, n_mismatch_score=-2, n_match_score=-1):
    letters = ['A','T','C','G','N']
    headers = [ord(x) for x in letters]
    mat_size = max(headers) + 1
    nuc_ords = [ord(x) for x in ['A','T','C','G']]
    a = np.zeros((mat_size, mat_size), dtype=np.int64)
    for nuc in nuc_ords:
        for nuc2 in nuc_ords:
            a[nuc,nuc2] = match_score if nuc == nuc2 else mismatch_score
    for nuc in nuc_ords:
        a[nuc,ord('N')] = n_mismatch_score
        a[ord('N'),nuc] = n_mismatch_score
    a[ord('N'),ord('N')] = n_match_score
    return a

def global_align_debug(pystr_seqj, pystr_seqi, matrix, gap_incentive_list, gap_open=-20, gap_extend=-2):
    seqj = pystr_seqj.encode('UTF-8')
    seqi = pystr_seqi.encode('UTF-8')
    gap_incentive = np.array(gap_incentive_list, dtype=np.int64)
    max_j = len(pystr_seqj)
    max_i = len(pystr_seqi)

    mScore = np.zeros((max_i + 1, max_j + 1), dtype=np.int32)
    iScore = np.zeros((max_i + 1, max_j + 1), dtype=np.int32)
    jScore = np.zeros((max_i + 1, max_j + 1), dtype=np.int32)
    mPointer = np.zeros((max_i + 1, max_j + 1), dtype=np.int32)
    iPointer = np.zeros((max_i + 1, max_j + 1), dtype=np.int32)
    jPointer = np.zeros((max_i + 1, max_j + 1), dtype=np.int32)

    min_score = gap_open * max_j * max_i

    mScore[0,1:] = min_score; mScore[1:,0] = min_score; mScore[0,0] = 0
    mPointer[0,1:] = IARRAY; mPointer[1:,0] = JARRAY; mPointer[0,0] = 0

    for i in range(1, max_j+1):
        iScore[0,i] = gap_extend * i + gap_incentive[0]
    iScore[0:,0] = min_score
    iPointer[0,1:] = IARRAY

    for i in range(1, max_i+1):
        jScore[i,0] = gap_extend * i + gap_incentive[0]
    jScore[0,0:] = min_score
    jPointer[1:,0] = JARRAY

    # inner cells
    for i in range(1, max_i):
        ci = seqi[i - 1]
        for j in range(1, max_j):
            cj = seqj[j - 1]
            iFromMVal = gap_open + mScore[i, j - 1] + gap_incentive[i]
            iExtendVal = gap_extend + iScore[i, j - 1] + gap_incentive[i]
            if iFromMVal > iExtendVal: iScore[i,j]=iFromMVal; iPointer[i,j]=MARRAY
            else: iScore[i,j]=iExtendVal; iPointer[i,j]=IARRAY

            jFromMVal = gap_open + mScore[i - 1, j] + gap_incentive[i-1]
            jExtendVal = gap_extend + jScore[i - 1, j]
            if jFromMVal > jExtendVal: jScore[i,j]=jFromMVal; jPointer[i,j]=MARRAY
            else: jScore[i,j]=jExtendVal; jPointer[i,j]=JARRAY

            mVal = mScore[i - 1, j - 1] + matrix[ci,cj]
            iVal = iScore[i - 1, j - 1] + matrix[ci,cj]
            jVal = jScore[i - 1, j - 1] + matrix[ci,cj]
            if mVal > jVal:
                if mVal > iVal: mScore[i,j]=mVal; mPointer[i,j]=MARRAY
                else: mScore[i,j]=iVal; mPointer[i,j]=IARRAY
            else:
                if jVal > iVal: mScore[i,j]=jVal; mPointer[i,j]=JARRAY
                else: mScore[i,j]=iVal; mPointer[i,j]=IARRAY

    # last column (j=max_j, i < max_i) - use gap_extend for new gap
    j = max_j; cj = seqj[j-1]
    for i in range(1, max_i):
        ci = seqi[i-1]
        iFromMVal = gap_extend + mScore[i, j - 1] + gap_incentive[i]
        iExtendVal = gap_extend + iScore[i, j - 1] + gap_incentive[i]
        if iFromMVal > iExtendVal: iScore[i,j]=iFromMVal; iPointer[i,j]=MARRAY
        else: iScore[i,j]=iExtendVal; iPointer[i,j]=IARRAY
        jFromMVal = gap_extend + mScore[i - 1, j] + gap_incentive[i-1]
        jExtendVal = gap_extend + jScore[i - 1, j]
        if jFromMVal > jExtendVal: jScore[i,j]=jFromMVal; jPointer[i,j]=MARRAY
        else: jScore[i,j]=jExtendVal; jPointer[i,j]=JARRAY
        mVal = mScore[i - 1, j - 1] + matrix[ci,cj]
        iVal = iScore[i - 1, j - 1] + matrix[ci,cj]
        jVal = jScore[i - 1, j - 1] + matrix[ci,cj]
        if mVal > jVal:
            if mVal > iVal: mScore[i,j]=mVal; mPointer[i,j]=MARRAY
            else: mScore[i,j]=iVal; mPointer[i,j]=IARRAY
        else:
            if jVal > iVal: mScore[i,j]=jVal; mPointer[i,j]=JARRAY
            else: mScore[i,j]=iVal; mPointer[i,j]=IARRAY

    # last row (i=max_i)
    i = max_i; ci = seqi[i - 1]
    for j in range(1, max_j+1):
        cj = seqj[j - 1]
        iFromMVal = gap_extend + mScore[i, j - 1] + gap_incentive[i]
        iExtendVal = gap_extend + iScore[i, j - 1] + gap_incentive[i]
        if iFromMVal > iExtendVal: iScore[i,j]=iFromMVal; iPointer[i,j]=MARRAY
        else: iScore[i,j]=iExtendVal; iPointer[i,j]=IARRAY
        jFromMVal = gap_extend + mScore[i - 1, j] + gap_incentive[i-1]
        jExtendVal = gap_extend + jScore[i - 1, j]
        if jFromMVal > jExtendVal: jScore[i,j]=jFromMVal; jPointer[i,j]=MARRAY
        else: jScore[i,j]=jExtendVal; jPointer[i,j]=JARRAY
        mVal = mScore[i - 1, j - 1] + matrix[ci,cj]
        iVal = iScore[i - 1, j - 1] + matrix[ci,cj]
        jVal = jScore[i - 1, j - 1] + matrix[ci,cj]
        if mVal > jVal:
            if mVal > iVal: mScore[i,j]=mVal; mPointer[i,j]=MARRAY
            else: mScore[i,j]=iVal; mPointer[i,j]=IARRAY
        else:
            if jVal > iVal: mScore[i,j]=jVal; mPointer[i,j]=JARRAY
            else: mScore[i,j]=iVal; mPointer[i,j]=IARRAY

    print(f"\nmScore final [{max_i},{max_j}] = {mScore[max_i,max_j]}")
    print(f"iScore final = {iScore[max_i,max_j]}")
    print(f"jScore final = {jScore[max_i,max_j]}")

    align_j_list = []; align_i_list = []; matchCount = 0
    i = max_i; j = max_j
    ci = seqi[i - 1]; cj = seqj[j - 1]
    currMatrix = MARRAY
    if mScore[i,j] > jScore[i,j]:
        if mScore[i,j] > iScore[i,j]: currMatrix = MARRAY
        else: currMatrix = IARRAY
    else:
        if jScore[i,j] > iScore[i,j]: currMatrix = JARRAY
        else: currMatrix = IARRAY

    print(f"Starting traceback from ({i},{j}), currMatrix={currMatrix}")

    while i > 0 or j > 0:
        if currMatrix == MARRAY:
            nxt = mPointer[i,j]
            print(f"  M({i},{j}) -> {nxt}, align {chr(cj)}/{chr(ci)}")
            currMatrix = nxt
            align_j_list.append(chr(cj)); align_i_list.append(chr(ci))
            if cj == ci: matchCount += 1
            if i > 1: i -= 1; ci = seqi[i - 1]
            else: i = 0; ci = seqi[i]
            if j > 1: j -= 1; cj = seqj[j - 1]
            else: j = 0; cj = seqj[j]
        elif currMatrix == JARRAY:
            nxt = jPointer[i,j]
            print(f"  J({i},{j}) -> {nxt}, gap_in_read/{chr(ci)}")
            currMatrix = nxt
            align_j_list.append('-'); align_i_list.append(chr(ci))
            if i > 1: i -= 1; ci = seqi[i - 1]
            else: i = 0; ci = seqi[i]
        elif currMatrix == IARRAY:
            nxt = iPointer[i,j]
            print(f"  I({i},{j}) -> {nxt}, {chr(cj)}/gap_in_ref")
            currMatrix = nxt
            align_j_list.append(chr(cj)); align_i_list.append('-')
            if j > 1: j -= 1; cj = seqj[j - 1]
            else: j = 0; cj = seqj[j]

    align_j = ''.join(reversed(align_j_list))
    align_i = ''.join(reversed(align_i_list))
    align_counter = len(align_j_list)
    final_score = 100 * matchCount / float(align_counter)
    return align_j, align_i, round(final_score, 3)

m = make_matrix(5, -4, -2, -1)
print("="*60)
print("Test: ANNG (read/seqj) vs ATCG (ref/seqi), gap_open=-20, gap_extend=-2")
print("Expected by Python tests: read=A-NNG, ref=ATC-G, score=40.0")
r = global_align_debug('ANNG', 'ATCG', m, [0,0,0,0,0], -20, -2)
print(f"\nResult: read={r[0]}, ref={r[1]}, score={r[2]}")
print()
print("="*60)
print("Test: ATTTA (read) vs ATTA (ref), gap_incentive=[0,0,1,0,0]")
print("Expected by Python tests: read=ATTTA, ref=AT-TA, score=80.0")
r2 = global_align_debug('ATTTA', 'ATTA', m, [0,0,1,0,0], -20, -2)
print(f"\nResult: read={r2[0]}, ref={r2[1]}, score={r2[2]}")
