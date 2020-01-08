import itertools, random

def Overlap(seq1,seq2):
    #The list containing overlapping prefixes and suffixes is built
    ov = [seq1[-i:] for i in range(min(len(seq1),len(seq2))) if seq1[-i:]==seq2[:i]]

    #Since we may not have found overlaps, we put a condition
    #In case there are any, the function finds the longest possible
    if len(ov)>0:
        max_overlap = ov[0] 
        for k in ov[1:]:
            if len(k)>len(max_overlap):
                max_overlap = k
        return len(max_overlap)
    else:
        return 0

def SCS(l, threshold=3):
    #SCS takes as argument a list of strings l, and threshold of overlap, fixed at 3, but can be modified 
    #If the list contains only one string, either we have reached the shortest common superstring
    #Or no other strings are available for merging
    if len(l)==1:
        return l[0]

    #Generates all pairs to find the best match
    pairs=[pair for pair in itertools.permutations(l,2)]
        
    #Creates a dictionary where keys are pairs of sequences
    #Values are the overlapping score associated to such pair    
    matches_pair = {i:Overlap(i[0],i[1]) for i in pairs}

    #Finds the best score that can be obtained by all the prefix-suffix overlap scores computed     
    maximum_score = max(list(matches_pair.values()))
    possibilities = [w for w in matches_pair if matches_pair[w] == maximum_score]

    #If the same maximum score has been obtained more than once,
    #The pair of sequences to be joined is randomly chosen
    maximum_score_tuple = (maximum_score, random.choice(possibilities))
        
    #Checks that the overlap is greater than the chosen threshold
    #If not so, the two strings get concatenated
    if maximum_score < threshold:
        new_string = maximum_score_tuple[1][0]+ maximum_score_tuple[1][1]
    else:
        new_string = maximum_score_tuple[1][0]+ maximum_score_tuple[1][1][maximum_score:]

    #Removes the pair of sequences that have given rise to the new merged string
    #The new string is then added at the beginning of the list
    l.remove(maximum_score_tuple[1][0]) 
    l.remove(maximum_score_tuple[1][1])
    l = [new_string] + l

    #Recursively, the function is called until we obtain only the shortest common superstring 
    return SCS(l)

'''To test the module'''
##def getGenome(length=1000):
##    genome=''.join(random.choice('AGCT') for i in range(length))
##    return genome
##
##def getSubstrings(seq, length=100):
##    L=[]
##    for i in range(len(seq)-length+1):
##        L.append(seq[i:i+length])
##    return L
##
###User inputs the length of the genome 
##test_genome=getGenome(int(input('Length of the genome: ')))
##
###User inputs the length of the substrings with which SCS is computed
##stringset=getSubstrings(test_genome,int(input('Length of the read: ')))
##
###The randomly generated genome is reported 
##print(test_genome)
##
###The shortest common superstring is generated and reported
##output=SCS(stringset)
##print(output)
##
###For consistency, it checks if the shortest common superstring corresponds.
###True = genome and SCS correspond, False = they don't.
##print(output==test_genome)




