''' Dini Alice | Assignment 3 | 830931 '''
import sys
def rotations(t):
    #Takes a string t as argument
    #Returns list of rotations of input string t
    t = t.upper()
    return [t[i:]+'$'+t[:i] for i in range(len(t)+1)]

def BWT_matrix(t):
    #Takes a string t as argument
    #Returns the Burrows-Wheeler matrix.
    #Computes the list of t's rotations, making up a matrix
    matrix = rotations(t)
    #.sort() method works directly on the list without creating a new one, hence faster.
    matrix.sort()
    return matrix

def BWT(matrix):
    #Takes the Burrows-Wheeler matrix as argument
    #Returns the Burrows-Wheeler transform
    return ''.join([matrix[i][-1] for i in range(len(matrix))])

def suffix_array(no_sorted,yes_sorted):
    #Takes the matrix made with the rotations not lexicographically sorted as first argument
    #Takes the Burrows-Wheeler matrix as second argument
    #Returns the suffix array.
    #By iterating in the BW matrix, each element is fetched into the rotations' matrix
    #The index at which such element is found is stored in the suffix array
    #Corresponding to the position of each suffix lexicographically sorted in the matrix.
    return [no_sorted.index(i) for i in yes_sorted]

def reverse(matrix):
    #Takes the Burrows-Wheeler matrix as argument
    #Storing the elements in the last row
    last = [k[-1] for k in matrix]
    #Uses the global variable SA (suffix array)
    #Putting together first column & offset, then storing
    #By default, sorted() sorts according to the numerical values it finds in the tuples
    list_ordered = list(zip(SA,last))
    list_ordered.sort()
    #Then the characters are extracted from the tuples to be arranged in the reversed BWT.
    return ''.join([list_ordered[i][1] for i in range(len(list_ordered))])

def rank_cumulative(bwt):
    #Takes the Burrows-Wheeler transform as argument
    #Stores as keys a character in the BWT, and as values a list of the cumulative number of times
    #it has found the character until that row, comprised that row.
    #In this way we always know how many characters up to one point have been found in BWT of that kind 
    #An auxiliary dictionary containing the number of times a character is present is used.
    #The auxiliary dictionary for the chars' occurrences is initialized
    #At the beginning, the number of occurrences is 0
    all_occurrences = {c:0 for c in bwt}
    #The dictionary for cumulative sums is initialized
    cumulative_sum = {c:[] for c in bwt}
    #By updating the auxiliary dictionary, each time we fill a position of the cumulative_sum dictionary's lists
    #The lists will have all the same length of the BWT standing for the different positions where the computations
    #arrive during the matching.
    for c in bwt:
        all_occurrences[c] += 1
        #Whenever a character is encountered in BWT, its number of occurrences is updated in the auxiliary dictionary
        #We fill each position of the corresponding value in cumulative_sum with the number of occurrences found so far
        #At each iteration, while in the auxiliary dictionary only the value of the considered character is incremented
        #In cumulative_sum we fill the correspondent position (of BWT) in each list with the occurrence of each char
        for c in all_occurrences.keys():
            cumulative_sum[c].append(all_occurrences[c])
    #Frees up memory for already used elements        
    del all_occurrences
    return cumulative_sum

def ordered_occurrences_ranges(matrix):
    #Takes the Burrows-Wheeler matrix as input
    #Creates a dictionary storing the range of rows where each character starts and ends
    #Keys are characters in the first column; values are lists storing the range. 
    F = [x[0] for x in matrix]
    rank = {}
    for i in range(len(F)):
        rank[F[i]] = [i,i+F.count(F[i])]
##        if F[i] not in rank:
##            rank[F[i]] = [i,i+1]
##        else:
##            rank[F[i]][1] += 1
    print(rank)
    return rank

def matching(p):
    p = p.upper()
    #Storing a dictionary with characters as keys and cumulative occurrences as values in lists
    CS = rank_cumulative(bwt)
    #Storing the ranges at which each letter in the input sequence occurs, in lexicographic order, in the first column
    F = ordered_occurrences_ranges(bwt_mat)
    #Checking if the last character of the query is present in the first column
    if p[-1] not in F:
        #If not so, no occurrences are present
        return 'No occurrences are present.'
    #Setting two cursors, one at the beginning and one at the end of the rows having as char the last char in p
    #We consider the range (value in the dictionary of the first column) that has the last character of the query as key
    start, stop = F[p[-1]]
    #Setting i as the second-to-last character since we altready checked the last one
    i = len(p)-2
    while i >= 0 and stop > start:
        #Two exit conditions: if i reaches a value smaller than 0, we traversed all the pattern, and an occurence has been found hopefully
        #If the starting position of the range of rows considered at the moment is greater or equal the final one, no occurrence has been found.
        #c assumes the character at index i, performing a backward search
        c = p[i]
        #start will assume the starting point of the new range of rows that will be considered, by fetching the starting position
        #Of such character occurrences in the first row, since they are in lexicographic order (ie all the characters prior to c); 
        #but according to the rank in CS dictionary, the precise row at which we need to put the starting point will be obtained,
        #Considering the characters "c" that have been encountered so far. The -1 is required to get how many "c" characters actually preceed
        #The character considered in the previous iteration to have a proper rank starting from 0.
        #for example, after the initial check where start and stop are set at F[p[-1]] and given the knowledge of which p[-2] characters are present in BWT,
        #we can calculate the new range, and directly move at that rows.
        #The same is performed for the end point of the range. the two cursors will indeed start counting from the beginning of the occurrences.
        #It updates every time the range in the first column that contains the range start-stop where the suffix p[i:] is found.
        start = F[c][0] + CS[c][start-1]
        stop = F[c][0] + CS[c][stop-1]
        #Going backward in the query
        i -= 1
    #Frees up memory for already used elements
    del F
    del CS
    #The list that will contain the occurrences' ranges is initialized
    occurrences = []
    #Given the lexicographic order, the occurrences will be found all one near another
    #Given the number of occurrences, we can compute the number of rows to be considered 
    #Hence we can map each row from the start to the end of the range to the suffix array
    #The reasoning is based on the correspondence between suffix array and first column indexes
    diff = stop - start
    z = start
    #z is the point from which we start mapping index, so where we arrived from the computation
    while z < start + diff:
        #occurrences is updated with the starting position of each occurrence
        #start is fetched from the corresponding index in the suffix array
        occurrences = occurrences + [(SA[z])]
        z+=1   
    if len(occurrences) == 0:
        return 'No occurrences are present.'
    #Return size of final range which corresponds to the number of occurrences
    #Returns the starting position at which the query is found.
    return str(diff)+' occurrence(s) found starting at position(s): '+str(sorted(occurrences))
               
#Output & global variables
if __name__ == "__main__":
    if len(sys.argv) !=3:
        print("Usage: %s sequence - pattern" % sys.argv[0]) 
        sys.exit()
    genome = sys.argv[1].upper()
    pattern = sys.argv[2].upper()
    #Storing a series of global variables that otherwise would be computed more than once
    #The Burrows-Wheeler transform is printed out.
    rot_mat = rotations(genome)
    bwt_mat = BWT_matrix(genome)
    bwt = BWT(bwt_mat)
    SA = suffix_array(rot_mat,bwt_mat)
    #Frees up memory for already used elements
    del rot_mat
    rev = reverse(bwt_mat)
    print('The Burrows-Wheeler transform of your genome is: ',bwt)
    #Checks whether the computation of the reversed Burrows-Wheeler transform works.
    if rev[1:] == genome.upper():
        print('The reverse Burrows-Wheeler transform has been correctly computed, and it is ',rev)
    else:
        print('Something went wrong with the computation of the reverse Burrows-Wheeler transform.')
    #Frees up memory for already used elements
    del genome
    del rev
    #Executes the matching without taking further parameters
    print(matching(pattern))
    #Frees up memory as a final clean up
    del bwt_mat
    del bwt
    del SA
    print('Done.')
                
        
