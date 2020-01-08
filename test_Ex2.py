import random, itertools, os, Ex2, Ex2_Bruteforce
'''This module allows to test both the greedy SCS algorithm and the bruteforce solution to the problem.
A "sequences.fasta" file can be taken as input: it should contain some genome sequences.
This is done for the mere purpose of testing the algorithms, so even if FASTA files generally contain reads,
in this case it is required an hypothetical FASTA file containing, line after line, the sequences that will
be first splitted into substrings; from these latters, the shortest common superstring is obtained.
If the file exists, the program allows to choose which one of the algorithms to use.
If such file doesn't exist, it creates a new fasta file according to the desires of the user,
randomly generating sequences requiring as input the length of the sequences, along with
the number of sequences of that length the user wants to generate.
In both the cases, the user can provide as input the length of the substrings that will be generated to
compute the strongest common supestring.
An evaluation of the consistency of the output sequence is also provided:
True = SCS is the same as the input; False = SCS is different from the input.

PS: To get reasonably fast results, when choosing to run the bruteforce
algorithm with this test, choose higher length of reads/substrings than
usual (e.g.for a sequence of 15 nucleotides, reads with length = 9 should be fine,
while with Ex2_Bruteforce also reads with length = 6 could be used).'''

def getGenome(length=1000):
    genome=''.join(random.choice('AGCT') for i in range(length))
    return genome

def getSubstrings(seq, length=100):
    L=[]
    for i in range(len(seq)-length+1):
        L.append(seq[i:i+length])
    return L

#Opens the file with the genomes
#Generates substrings for each of the genome
#Runs the SCS greedy algorithm on each of the sequences
#It finally checks if the reconstructed sequences are the same as the input ones
def readFASTA1():
    with open('sequences.fasta') as rd:
        contents = rd.readlines()
        for i in contents:
            result = Ex2.SCS(getSubstrings(i,readlen))
            print(result,result==i)
    return 'Done'

#Opens the file with the genomes
#Generates substrings for each of the genome
#Runs the SCS greedy algorithm - Bruteforce approach on each of the sequences
#It finally checks if the reconstructed sequence is the same
def readFASTA2():
    with open('sequences.fasta') as rd:
        contents = rd.readlines()
        for i in contents:
            result = Ex2_Bruteforce.SCS(getSubstrings(i,readlen))
            print(result,result==i) #checks if the reconstructed sequence is the same  
    return 'Done'

#Storing the name for the file to be fetched 
file_name = "sequences.fasta"

#Gets current directory
cur_dir = os.getcwd()

#Lists the files contained in the current directory
file_list = os.listdir(cur_dir)

#Checks if it exists a file named "sequences.fasta" to be used
#If the file does not exist, it creates a new one
#The user indicates:
#-The length of the genomes
#-The number of genomes with that length to be analysed
if file_name not in file_list:
    genseq=int(input('Length of the genome: '))
    nseq=int(input('How many random sequences with this length do you want? '))
    file=open('sequences.fasta','w+') 
    for i in range(nseq):
        file.write(getGenome(genseq+1)+'\n')
    file.close()

#The user can choose which approach he wants to use    
choice=input('Select 1 to run the greedy algorithm. Select 2 to run the bruteforce greedy algorithm. ')
readlen=int(input('Length of the reads: '))
if choice=='1':
    print(readFASTA1())
else:
    print(readFASTA2())
    
