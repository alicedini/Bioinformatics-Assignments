import sys

def empty_matrix(rows, columns): #making a matrix of zeros
    zero_matrix = [[0 for x in range(columns)] for y in range(rows)] # adds a zero to each column in each row at the end
    return zero_matrix
#notice that with respect to the N&W algorithm we don't initialize progressive gap penalties

def score_function(char1, char2): #determining the score between any two bases in alignment
    if char1 == char2:
        return match_score
    else:
        return mismatch_score

def needleman_wunsch(seq1, seq2):
    n = len(seq1) #on the columns
    m = len(seq2) #on thw rows
    score = empty_matrix(m+1, n+1)
    for i in range(1, m + 1):  # calculating the scores for each of the other cells
        for j in range(1, n + 1): # checks top, left, and diagonal cells 
            move_diagonally = score[i - 1][j - 1] + score_function(seq1[j-1], seq2[i-1])  
            move_leftward = score[i - 1][j] + gap 
            move_upward = score[i][j - 1] + gap  
            score[i][j] = max(move_diagonally, move_leftward, move_upward)
            if score[i][j] < 0: #additional condition of S&W 
                score[i][j] = 0
    
#### TRACEBACK AND ALIGNMENT #### 
    
    # store alignment
    align1 = ""
    align2 = ""
    initial_position_j = 0 #finding the maximum score of the matrix and then storing its coordinates from which traceback will start
    initial_position_i = 0
    best_score = 0
    for j in range(n+1):
        for i in range(m+1):
            if score[i][j] > best_score:
                best_score = score[i][j]
                initial_position_i = i
                initial_position_j = j
    print('The best local alignment score is ',best_score)            
    i = initial_position_i #setting the coordinates of the maximum as starting position of my traceback cursor
    j = initial_position_j
    while i > 0 and j > 0: 
        current_position = score[i][j]
        diagonal = score[i-1][j-1] 
        up = score[i][j-1] 
        left = score[i-1][j] 
        
        # backtracking to which cell the current score was calculated from, then i set the coordinates of that cell to continue alignment
        #whenever it finds a gap, it puts a -, but this time, whenever it arrives at a position with score 0 it ends the execution (theoretically, this occurs when the cells around mine are 0)
        if current_position == 0:
            align1 = align1[::-1] #because we had backtracked 
            align2 = align2[::-1]
            
            for e in score:
                print(e)
            return(align1, align2)
            
        if current_position == diagonal + score_function(seq1[j-1], seq2[i-1]):
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
            
    # top left cell :whenever we arrive at column or row 0 we continue traversing that line knowing that there will be all gaps
    while j > 0:
        align1 += seq1[j-1]
        align2 += '-'
        j -= 1
    while i > 0:
        align1 += '-'
        align2 += seq2[i-1]
        i -= 1  
    align1 = align1[::-1] #because we had backtracked 
    align2 = align2[::-1]

    for e in score:
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
    command_sequence1, command_sequence2 = needleman_wunsch(seq2, seq1) #swap done on purpose

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
            elif command_sequence1[c] == '-' or command_sequence2[v] == '-': #-=gap on either sequence
                chars_alignment+=' '
                c+=1
                v+=1
            else:
                chars_alignment+='|' #|=mismatch
                c+=1
                v+=1
    print(chars_alignment)
    print(command_sequence1)
