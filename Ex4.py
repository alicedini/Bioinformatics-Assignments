''' Dini Alice | Assignment 4 | 830931 '''

import numpy, math, os, Ex4_models
from tabulate import tabulate

def query(s):
    '''Given a sequence s, query(s) computes its log-ratio score
    given the inside and the outside models to predict if
    it pertains to a CpG island or not.'''
    
    #Checks if the models are already present in the directory where this file
    #is stored: if not so, it computes them importing the function "models()".
    cur_dir = os.getcwd()
    file_list = os.listdir(cur_dir)
    file_name1, file_name2 = 'inside_file.npy','outside_file.npy'
    if file_name1 not in file_list or file_name2 not in file_list: Ex4_models.models()

    s = s.upper()

    #The models are loaded for the computations.
    inside, outside = numpy.load('inside_file.npy'), numpy.load('outside_file.npy')
    inside, outside = numpy.array(inside), numpy.array(outside)

    #Conversion of numpy arrays to Python's lists.
    inside, outside = inside.tolist(), outside.tolist()

    if len(s) == 1: return 0.25
    else:
        prob_in, prob_out  = math.log(0.25), math.log(0.25)                   
        #P(C|A) = prob of AC, coordinates in dictionary are fetched row by column;
        #In the models, each cell contains the probability of having a certain
        #nucleotide given the previous one, hence probability of dinucleotides.
        coordinates_matrix = {'AA':[0,0], 'AC':[0,1], 'AG':[0,2], 'AT':[0,3], 'CA':[1,0], 'CC':[1,1], 'CG':[1,2], 'CT':[1,3], 'GA':[2,0], 'GC':[2,1], 'GG':[2,2], 'GT':[2,3], 'TA':[3,0], 'TC':[3,1], 'TG':[3,2], 'TT':[3,3]}
        i = 1
        #Scanning the sequence for the dimers, starting at the second nucleotide.
        while i <= len(s)-2:
            coord = coordinates_matrix[s[i:i+2]]
            #Probabilities are gathered from the dictionary thanks to row by column mapping.
            prob_in+=math.log(inside[coord[0]][coord[1]])
            prob_out+=math.log(outside[coord[0]][coord[1]])
            i+=1

    #Properties of logarithms.        
    ratio = prob_in - prob_out
    return ratio

def scan_genome(genome, average_CpG = 567):
    '''The computation of the default average CpG length has been computed using
    the lines of codes exposed in Ex_models.py, but the user is able to
    change this value in order to use the window he wants to scan the input sequence.'''
    
    if len(genome) < average_CpG:
        average_CpG = int(input('Since your genome is smaller than the default CpG island average length (567), insert the windows length you would like to use. '))
        while average_CpG <= 0:
            average_CpG = int(input('A sliding window must have a size greater than zero. Insert a new value. '))
    windows = []
    i = 0
    #number_of_windows counts how many windows are actually checked inside the query.
    number_of_windows = 0
    #The while loop proceeds until there's not enough space to fit a window.
    while i + average_CpG <= len(genome):
        number_of_windows += 1
        window = genome[i:i + average_CpG]
        score = query(window)
        index = i
        windows = windows + [(index,score,window)]
        i+=1
    #positives contains all the windows of the genome which have positive score.
    positives = [[j[0], j[1]] for j in windows if j[1]>0]
    return positives, number_of_windows

#User interaction
genome=input('Insert the genome: ')
choice=input('Press 1 to obtain the prediction of your sequence being inside or outside a CpG island, or press 2 to visualize the distribution of scores according to a fixed window. In this latter case, a sequence longer than 1000 nucleotides is suggested, in order to obtain an appropriate result.')

if not genome: print('Empty sequence inserted.')
elif choice == '1':
    result = query(genome)
    print('The log-ratio score of your sequence given the inside and outside models previously computed is: ',result)
    if result > 0:
        print('Your sequence is more likely to be inside a CpG island.')
    elif result < 0:
        print('Your sequence is more likely to be outside a CpG island.')
    else:
        print('The probability of your sequence of being inside a CpG island equals the probability of being outside one of them.')
else:
    tot_result = scan_genome(genome)
    result = tot_result[0]
    print('Out of ',tot_result[1],' scanned windows, ',len(result),' positive scores have been obtained.')
    #tabulate is the function that allows to output the data with the desired layout.
    print(tabulate(result, headers = ['Offset', 'Score'], tablefmt = 'orgtbl'))



