import sys

def empty_matrix(rows, columns): #making a matrix of zeros
    zero_matrix = [[0 for x in range(columns)] for y in range(rows)]
    return zero_matrix

def score_function(char1, char2): #determining the score between any two bases in alignment using the arbitrary score gathered from the user
    if char1 == char2:
        return match_score
    else:
        return mismatch_score

def needleman_wunsch(seq1, seq2):
    n = len(seq1) #columns 
    m = len(seq2) #rows
    score = empty_matrix(m+1, n+1)
    for i in range(0, m + 1): # initialization: first row with progressive penalties
        score[i][0] = gap * i
    for j in range(0, n + 1): # initialization: first column with progressive penalties
        score[0][j] = gap * j
    for i in range(1, m + 1):  # calculating the scores for each of the other cells
        for j in range(1, n + 1): # to fill the cell the function computes the score from diagonal, left or top cell and the maximum among them
            move_diagonally = score[i - 1][j - 1] + score_function(seq1[j-1], seq2[i-1])  
            move_leftward = score[i - 1][j] + gap 
            move_upward = score[i][j - 1] + gap  
            score[i][j] = max(move_diagonally, move_leftward, move_upward)

#### TRACEBACK AND ALIGNMENT #### 
    
    # initialization of the future strings to output the alignment
    align1 = ""
    align2 = ""
    #we start the alignment from the bottom left cell hence we set the dimensions of the matrix as coordinates
    i = m 
    j = n 
    print('The best global alignment score is ',score[m][n])
    while i > 0 and j > 0: # end touching the top or the left edge
        current_position = score[i][j] #initializing a cursor over the matrix located on the last cell
        diagonal = score[i-1][j-1] 
        up = score[i][j-1] #add gap to sequence on rows
        left = score[i-1][j] #add gap to sequence on columns
        
        # backtracking to which cell the current score was calculated from, then i set the coordinates of that cell to continue alignment
       
        if current_position == diagonal + score_function(seq1[j-1], seq2[i-1]): #if the actual score is given by the diagonal, move diagonally = add residue from both sequences
            align1 += seq1[j-1]
            align2 += seq2[i-1]
            i -= 1
            j -= 1
        elif current_position == up + gap: #we continue the alignment from the above cell 
            align1 += seq1[j-1]
            align2 += '-'
            j -= 1
        elif current_position == left + gap: #we continue the alignment from the cell on the left 
            align1 += '-'
            align2 += seq2[i-1]
            i -= 1

    # whenever we arrive at column or row 0 we continue traversing that line knowing that there will be all gaps
    while j > 0: #we are on the first column so we can only move upward
        align1 += seq1[j-1]
        align2 += '-'
        j -= 1
    while i > 0: #we are on the first row so we can only move leftward
        align1 += '-'
        align2 += seq2[i-1]
        i -= 1  
    align1 = align1[::-1] #because we had backtracked 
    align2 = align2[::-1]

    for e in score: #to print the matrix
        print(e)
    return(align1, align2)

if __name__ == "__main__":
    if len(sys.argv) !=6:
        print("Usage: %s sequence1 - sequence2 - match score - mismatch score - gap" % sys.argv[0]) 
        sys.exit()
    seq1 = sys.argv[1].upper()
    seq2 = sys.argv[2].upper()
    match_score = int(sys.argv[3])
    mismatch_score = int(sys.argv[4])
    gap = int(sys.argv[5])
    command_sequence1, command_sequence2 = needleman_wunsch(seq2, seq1) #the swap of the two sequences has been done on purpose for debugging
    
    #layout for the output
    
    print(command_sequence2)
    c = 0
    v = 0
    chars_alignment = '' 
    while c < len(command_sequence1):
        while v < len(command_sequence2):
            if command_sequence1[c] == command_sequence2[v]:
                chars_alignment+='*' #*=match
                c+=1
                v+=1
            elif command_sequence1[c] == '-' or command_sequence2[v] == '-':#-=gap on either sequence
                chars_alignment+=' '
                c+=1
                v+=1
            else:
                chars_alignment+='|' #|=mismatch
                c+=1
                v+=1
    print(chars_alignment)
    print(command_sequence1)
