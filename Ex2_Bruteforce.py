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
    
def singular_SCS(l):
    #singular_SCS takes as argument a list of strings l
    #Computes the algorithm in one possible set, iterating from left to right
    #If the list contains only one string, either we have reached the shortest common superstring
    #Or no other strings are available for merging
    if len(l)>1:
        match = Overlap(l[0],l[1])

        #Merges the first two strings according to their prefix-suffix overlap
        new_string = l[0] + l[1][match:]

        #Remove the first merged string 
        l.remove(l[0])

        #Removes the second merged string
        #What was at index 1 now becomes index 0
        l.remove(l[0]) 
        l = [new_string] + l

        #Recursively calls the function 
        return singular_SCS(l)    
    else:
        return l

def SCS(l):
    #Creates all the possible permutations
    #List type elements are modifiable, allowing .remove() method
    total_permut = [list(sequence) for sequence in itertools.permutations(l)]

    #Calls the function singular_SCS for each of the strings' permutations
    #Obtains all the possible shortest common supestrings that can be obtained 
    #Computes the best, shortest common supestring 
    #Initializing as best shortest superstring the first one
    best = singular_SCS(total_permut[0])
    for i in total_permut[1:]:
        comparate = singular_SCS(i)
        if len(comparate) < len(best):
            best =  comparate
    return ''.join(best)

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

    
    
        
