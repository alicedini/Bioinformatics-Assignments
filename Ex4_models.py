''' Dini Alice | Assignment 4 - Inside & Outside models | 830931 '''
'''Warning: different modules are required to properly execute the routines.'''

import numpy, random, statistics, re, gzip
from urllib.request import urlopen

def models():
    #####Inside model#####
    '''For every line of the textfile containing all the coordinates
    we first extract only those referring to chromosome 22.
    Splitting them allows to obtain clear separation of data:
    since the text file is tab separated, we select data
    in column 1 & 2 referring to start and stop.
    Since the files are directly taken from the web, they have to be decoded
    using UTF-8; in particular, the FASTA file storing chromosome 22 requires
    also to be unzipped.'''

    #Fetches chromosome 22 FASTA zipped file.
    chr22 = urlopen('ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr22.fa.gz')
    
    #Extracts the file and stores it as lines of bytes, in a list; the header is excluded.
    with gzip.open(chr22, 'rb') as file:
        chr22_unzipped = file.readlines()[1:]
        
    #All the lines are now decoded in string format, set uppercase,
    #and cleaned from newline characters.
    stringlist = [x.decode('utf-8') for x in chr22_unzipped]
    chr22 = ''.join(stringlist).upper().replace('\n','')

    #Fetches CpG islands coordinates txt file (tab separated type) and directly stores all the lines of bytes in a list.
    coordinates = urlopen('http://www.haowulab.org/software/makeCGI/model-based-cpg-islands-hg19.txt').readlines()
    #Again all the lines are decoded as before.
    stringlist = [x.decode('utf-8') for x in coordinates]
    #All the coordinates which pertain to CpG islands in chromosome 22 are stored.
    #The columns at which coordinates are found in the txt file correspond to precise positions (indexes) in the elements of the list.
    all_coordinates = [(i[1],i[2]) for i in [i.split() for i in stringlist if i[0:5]=='chr22']]
    
    sequences_inside_model = []
    for i in all_coordinates:
        start = int(i[0])
        stop = int(i[1])+1
        sequences_inside_model = sequences_inside_model + [chr22[start:stop]] 
    
    #####Outside model#####
    len_chr = len(chr22)
    length_CpGs = [len(i) for i in sequences_inside_model]

    #####Computation of the average CpG length#####
    '''These lines of code have been used to compute the average, which
    is used by the above function "scan_genome". The user is free to change
    such parameter.
    
    average_CpG = int(statistics.mean(length_CpGs))
    '''
    
    sequences_outside_model = []
    #Random sequences are gathered from chromosome 22 to be used for the construction
    #of the outside model; as many random sequences as many CpG islands considered.
    for i in length_CpGs:
        start = random.randint(0,len_chr)
        end = start + i
        #Ns are removed by convention.
        sequences_outside_model = sequences_outside_model + [chr22[start:end].replace('N','')]

    #Calls the function found below.
    inside = frequencies(sequences_inside_model)
    outside = frequencies(sequences_outside_model)

    #Once the models have been computed, they are saved in the same directory of this file.
    numpy.save('inside_file', inside)
    numpy.save('outside_file', outside)

    return inside, outside

def frequencies(sequences):
    '''This function takes a list of sequences as argument, and for each of them
    it computes the frequency of each possible dinucleotide; then, it computes
    the probability of finding a certain nucleotide given the previous one, which is
    known, as relative frequency. The list of sequences undergoing the computations should make up
    a larger genomic sequence of which the overall absolute frequencies are computed,
    to finally yield the matrix containing in each cell the probability Y nucleotide
    given X (hence the matrix has on rows the nucleotide X which is given and came before Y on columns:
    P(base|previous base), to compute the probability for a sequence to be inside or outside a CpG island.'''

    #Absolute frequencies are computed as all the sequences taken as argument made up a single chromosome or genome.
    
    totAA = sum([len(re.findall('(?=AA)',i)) for i in sequences])
    totAT = sum([len(re.findall('(?=AT)',i)) for i in sequences])
    totAG = sum([len(re.findall('(?=AG)',i)) for i in sequences])
    totAC = sum([len(re.findall('(?=AC)',i)) for i in sequences])

    totAX = totAA + totAT + totAG + totAC

    totTA = sum([len(re.findall('(?=TA)',i)) for i in sequences])
    totTT = sum([len(re.findall('(?=TT)',i)) for i in sequences])
    totTG = sum([len(re.findall('(?=TG)',i)) for i in sequences])
    totTC = sum([len(re.findall('(?=TC)',i)) for i in sequences])

    totTX = totTA + totTT + totTG + totTC
    
    totGA = sum([len(re.findall('(?=GA)',i)) for i in sequences])
    totGT = sum([len(re.findall('(?=GT)',i)) for i in sequences])
    totGG = sum([len(re.findall('(?=GG)',i)) for i in sequences])
    totGC = sum([len(re.findall('(?=GC)',i)) for i in sequences])

    totGX = totGA + totGT + totGG + totGC
    
    totCA = sum([len(re.findall('(?=CA)',i)) for i in sequences])
    totCT = sum([len(re.findall('(?=CT)',i)) for i in sequences])
    totCG = sum([len(re.findall('(?=CG)',i)) for i in sequences])
    totCC = sum([len(re.findall('(?=CC)',i)) for i in sequences])

    totCX = totCA + totCT + totCG + totCC

    
    #The relative sequences are computed by dividing the absolute frequency of each
    #Dinucleotide by the absolute frequencies of all dinucleotides
    #With a particular one.
    
    #I row - The first nucleotide is on the columns; the second nucleotide on rows
    P_AA = totAA / totAX
    P_CA = totAC / totAX
    P_GA = totAG / totAX
    P_TA = totAT / totAX

    #II row
    P_AC = totCA / totCX
    P_CC = totCC / totCX
    P_TC = totCT / totCX
    P_GC = totCG / totCX

    #III row
    P_AT = totTA / totTX
    P_CT = totTC / totTX            #Example: C is on column and T on row
    P_TT = totTT / totTX
    P_GT = totTG / totTX

    #IV row
    P_AG = totGA / totGX
    P_CG = totGC / totGX
    P_TG = totGT / totGX
    P_GG = totGG / totGX

    probabilities = [[P_AA, P_CA, P_GA, P_TA],[P_AC, P_CC, P_GC, P_TC],[P_AG, P_CG, P_GG, P_TG],[P_AT, P_CT, P_GT, P_TT]]

    #The values are rounded off by convention, keeping just 2 decimals.
    matrix = numpy.array(probabilities)
    return numpy.around(matrix, decimals=2)
