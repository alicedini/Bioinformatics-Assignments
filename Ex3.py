''' Dini Alice | Assignment 3 | 830931 '''
import sys
def BWT_Rev_Match(t,p):
    t.upper()
    #Returns list of rotations of input string t
    rotations = [t[i:]+'$'+t[:i] for i in range(len(t)+1)]
    #####Computation of the Burrows-Wheeler matrix#####
    bwt_matrix = sorted(rotations)
    #####Computation of the Burrows-Wheeler transform#####
    bwt = ''.join([bwt_matrix[i][-1] for i in range(len(bwt_matrix))])
    print('The Burrows-Wheeler transform of your genome is: ',bwt)
    #####Computation of the Suffix Array#####
    SA=[rotations.index(i) for i in bwt_matrix]
    #####Computation of the reverse Burrows-Wheeler transform#####
    last = [k[-1] for k in bwt_matrix]
    list_ordered = list(zip(SA,last))
    list_ordered.sort()
    reverse = ''.join([list_ordered[i][1] for i in range(len(list_ordered))])
    if reverse[1:] == t:
        print('The reverse Burrows-Wheeler transform has been correctly computed, and it is ',reverse)
    else:
        print('Something went wrong with the computation of the reverse Burrows-Wheeler transform.')
    #####Computation of the string matching#####
    p = p.upper()
    CS = rank_cumulative(bwt)
    F = ordered_occurrences_ranges(bwt_matrix)
    if p[-1] not in F:
        return 'No occurrences are present.'
    start, stop = F[p[-1]]
    i = len(p)-2
    while i >= 0 and stop > start:
        c = p[i]
        start = F[c][0] + CS[c][start-1]
        stop = F[c][0] + CS[c][stop-1]
        i -= 1
    del F
    del CS
    occurrences = []
    diff = stop - start
    z = start
    while z < start + diff:
        occurrences = occurrences + [(SA[z])]
        z+=1   
    if len(occurrences) == 0:
        return 'No occurrences are present.'
    return str(diff)+' occurrence(s) found starting at position(s): '+str(sorted(occurrences))

def rank_cumulative(bwt):
    #Takes the Burrows-Wheeler transform as argument
    #Stores as keys a character in the BWT, and as values a list of the cumulative number of times
    #it has found the character until that row, comprised that row.
    all_occurrences = {c:0 for c in bwt}
    cumulative_sum = {c:[] for c in bwt}
    for c in bwt:
        all_occurrences[c] += 1
        for c in all_occurrences.keys():
            cumulative_sum[c].append(all_occurrences[c])       
    del all_occurrences
    return cumulative_sum

def ordered_occurrences_ranges(matrix):
    #Takes the Burrows-Wheeler matrix as input
    #Creates a dictionary storing the range of rows where each character starts and ends
    #Keys are characters in the first column; values are lists storing the range.
    F = [x[0] for x in matrix]
    rank = {}
    for i in range(len(F)):
        if F[i] not in rank:
            rank[F[i]] = [i,i+1]
        else:
            rank[F[i]][1] += 1
    return rank

#Output 
if __name__ == "__main__":
    if len(sys.argv) !=3:
        print("Usage: %s sequence - pattern" % sys.argv[0]) 
        sys.exit()
    genome = sys.argv[1].upper()
    pattern = sys.argv[2].upper()
    print(BWT_Rev_Match(genome, pattern))
    print('Done.')
                
        
