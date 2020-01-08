Dini Alice - Assignment no. 3

The file Ex3.py provides the output of:
-the Burrows-Wheeler transform
-a check on the consistency of the function that computes the reverse BWT
-the number of occurrences of the pattern string in the main sequence
-the positions where the eventual occurrences have been found
-if no occurrences exist, a message will report that
-if the function that computes the reverse BWT doesn't provide the expected result, 
a message will report that as well.

The file can be run by the command line from the directory of storage, calling:

Ex3.py [sequence] [query] 

sequence is the string on which all the computations will be performed;
query is the pattern string which will be checked against the first one for the matching.
No quotes are needed. The program is case insensitive.

Below, a detailed description of what each piece of code does can be read, if something results unclear from the code.

def BWT_Rev_Match(t,p):
    #Returns list of rotations of input string t
    
    #####Computation of the Burrows-Wheeler matrix#####
    #Returns the Burrows-Wheeler matrix, which is just the sorted list of the string t's rotations

    
    #####Computation of the Burrows-Wheeler transform#####
    #Returns the Burrows-Wheeler transform selecting the last column's elements of BWT matrix

    
    #####Computation of the Suffix Array#####
    #By iterating in the BW matrix, each element is fetched into the rotations' matrix
    #The index at which such element is found in this latter is stored in the suffix array
    #Corresponding to the position of each suffix lexicographically sorted in the matrix.
    
    #####Computation of the reverse Burrows-Wheeler transform#####
    #Storing the elements in the last row of the BW matrix
    #Uses the suffix array, putting together first column & offset, then storing
    #By default, sorted() sorts according to the numerical values it finds in the tuples

    #####Computation of the string matching#####

    #CS : Storing a dictionary with characters as keys and cumulative occurrences as values in lists

    #F : Storing the ranges at which each letter in the input sequence occurs, in lexicographic order, in the first column

    #Checking if the last character of the query is present in the first column
        #If not so, no occurrences are present

    #Setting two cursors, one at the beginning and one at the end of the rows having as char the last char in p

    #We consider the range (ie value) in F dictionary that has the last character of the query as key
    #Setting i as the second-to-last character since we altready picked the last one

        #Two exit conditions: if i reaches a value smaller than 0, we traversed all the pattern, and an occurence has been found
        #If the starting position of the range of rows considered at the moment is greater or equal the final one, no occurrence has been found.
        #c assumes the character at index i, performing a backward search
	
	#start is everytime updated considering both the starting position of the row that start with the character p[i] and the number of occurrences
	#of such character found up to that moment in BWT
	#so by knowing how many times that character was found in the BWT, for the LF mapping property we can traduce occurrences in a range
        #The -1 is required to get where "c" row occurrences actually start, to have a proper rank to be used as a range to start, and not a cumulative sum.
        #The same is performed for the end point of the range, which will reach the end-point of the rows we have to actually consider
	#both start and stop start counting from the beginning of the occurrences' rows in the first row of BW matrix, given the lexicographic sorting
	#we can use ranges since occurrences of a pattern will always be found one close to the other in the BW matrix
        #In simple words, it updates every time the range in the first column that contains the range start-stop where the suffix p[i:] is found

        #Going backward in the query

    #del frees up memory for already used elements
 
    #occurrences: the list that will contain the occurrences' offsets is initialized

    #Given the lexicographic order, the occurrences will be found all one near another
    #Given the number of occurrences, we can compute the number of rows to be considered 
    #Hence we can map each row from the start to the end of the range to the suffix array
    #The reasoning is based on the correspondence between suffix array and first column indexes
 
    #z is the point from which we start mapping index, so where the occurrences start in the BW matrix
 
        #occurrences is updated with the starting position of each occurrence
        #start is fetched from the corresponding index in the suffix array
        
    #Return size of final range which corresponds to the number of occurrences
    #Returns the starting position at which the query is found.

def rank_cumulative(bwt):
    #Takes the Burrows-Wheeler transform as argument
    #Stores as keys a character in the BWT, and as values a list of the cumulative number of times
    #it has found the character until that row, comprised that row.
    #In this way we always know how many characters up to one point have been found in BWT of that kind 
    #An auxiliary dictionary containing the number of times a character is present is used.
    #The auxiliary dictionary for the chars' occurrences is initialized
    #At the beginning, the number of occurrences is 0
    
    #The dictionary for cumulative sums is initialized

    #By updating the auxiliary dictionary, each time we fill a position of the cumulative_sum dictionary's lists
    #The lists will have all the same length of the BWT; each index corresponds to the different positions where the computations
    #arrived during the matching, in the start and stop values characterizing the range of rows considered.
 
        #Whenever a character is encountered in BWT, its number of occurrences is updated in the auxiliary dictionary
        #We fill each position of the corresponding value in cumulative_sum with the number of occurrences found so far
        #At each iteration, while in the auxiliary dictionary only the value of the considered character is incremented
        #In cumulative_sum we fill the correspondent position (of BWT) in each list with the occurrence of each char
       
    #Frees up memory for already used elements        
   

def ordered_occurrences_ranges(matrix):
    #Takes the Burrows-Wheeler matrix as input
    #Creates a dictionary storing the range of rows where each character starts and ends
    #Keys are characters in the first column; values are lists storing the range. 

 
