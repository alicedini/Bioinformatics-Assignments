import numpy, random, math, os, statistics
#global variable required for the windows.
average_CpG = 0

def query(s):
    cur_dir = os.getcwd()
    file_list = os.listdir(cur_dir)
    file_name1 = 'inside_file.npy'
    file_name2 = 'outside_file.npy'
    if file_name1 not in file_list or file_name2 not in file_list:
        models()
        
    s = s.upper()
    inside = numpy.load('inside_file.npy')
    outside = numpy.load('outside_file.npy')
    inside = numpy.array(inside)*3
    outside = numpy.array(outside)*3
    inside = inside.tolist()
    outside = outside.tolist()

    if len(s) == 0:
        return 'No sequence'
    elif len(s) == 0:
        return 0.25
    else:
        prob_in = 0.25*3 #'REVERSE' order because p(C|A) = prob of AC, dict row by column coord
        prob_out = 0.25*3
        coordinates_matrix = {'AA':[0,0], 'AC':[0,1], 'AG':[0,2], 'AT':[0,3], 'CA':[1,0], 'CC':[1,1], 'CG':[1,2], 'CT':[1,3], 'GA':[2,0], 'GC':[2,1], 'GG':[2,2], 'GT':[2,3], 'TA':[3,0], 'TC':[3,1], 'TG':[3,2], 'TT':[3,3]}
        i = 1
        while i <= len(s)-2:
            coord = coordinates_matrix[s[i:i+2]]
            prob_in*=inside[coord[0]][coord[1]]
            prob_out*=outside[coord[0]][coord[1]]
            i+=1
    prob_in = prob_in/3
    prob_out = prob_out/3
    ratio = prob_in/prob_out
    prediction = math.log(ratio)
    return prediction

def scan_genome(genome):
    models()
    windows = []
    i = 0
    while i+average_CpG<=len(genome):
        window = genome[i:i+average_CpG]
        score = query(window)
        index = i
        windows = windows + [(index,score,window)]
        i+=1

    positives = []
    for j in windows:
        if j[1]>0:
            positives = positives + ['Offset: ',j[0],' with score: ',round(j[1],2)]
    return positives
        
def models():
    #####Inside model#####
    #For every line of the textfile containing all the coordinates
    #we first extract only those referring to chromosome 22
    #splitting them allows to obtain clear separation of data
    #since the text file is tab separated, we select data
    #in column 1 & 2 referring to start and stop
   
    file = open('chr22.fa','r')
    complete_chr22 = file.readlines()[1:]
    chr22 = ''.join(complete_chr22).upper().replace('\n','')

    coordinates = open("model-based-cpg-islands-hg19.txt", "r")
    all_coordinates = [(i[1],i[2]) for i in [i.split() for i in coordinates.readlines() if i[0:5]=='chr22']]

    sequences_inside_model = []
    for i in all_coordinates:
        start = int(i[0])
        stop = int(i[1])+1
        sequences_inside_model = sequences_inside_model + [chr22[start:stop]]
    print(sequences_inside_model[0])    

    #####Outside model#####
    global average_CpG
    length_CpGs = [len(i) for i in sequences_inside_model]
    average_CpG = int(statistics.mean(length_CpGs)) 
    
    len_chr = len(chr22)
    
    sequences_outside_model = []
    
    for i in length_CpGs:
        start = random.randint(0,len_chr)
        end = start + i
        sequences_outside_model = sequences_outside_model + [chr22[start:end]]

    inside = frequencies(sequences_inside_model)
    outside = frequencies(sequences_outside_model)
    numpy.save('inside_file', inside)
    numpy.save('outside_file', outside)

    return inside, outside

def frequencies(sequences):

    totAA = sum([i.count('AA') for i in sequences])
    totAT = sum([i.count('AT') for i in sequences])
    totAG = sum([i.count('AG') for i in sequences])
    totAC = sum([i.count('AC') for i in sequences])

    totAX = totAA + totAT + totAG + totAC

    totTA = sum([i.count('TA') for i in sequences])
    totTT = sum([i.count('TT') for i in sequences])
    totTG = sum([i.count('TG') for i in sequences])
    totTC = sum([i.count('TC') for i in sequences])

    totTX = totTA + totTT + totTG + totTC
    
    totGA = sum([i.count('GA') for i in sequences])
    totGT = sum([i.count('GT') for i in sequences])
    totGG = sum([i.count('GG') for i in sequences])
    totGC = sum([i.count('GC') for i in sequences])

    totGX = totGA + totGT + totGG + totGC
    
    totCA = sum([i.count('CA') for i in sequences])
    totCT = sum([i.count('CT') for i in sequences])
    totCG = sum([i.count('CG') for i in sequences])
    totCC = sum([i.count('CC') for i in sequences])

    totCX = totCA + totCT + totCG + totCC

    #I row
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
    P_CT = totTC / totTX #C is on column and T on rows
    P_TT = totTT / totTX
    P_GT = totTG / totTX

    #IV row
    P_AG = totGA / totGX
    P_CG = totGC / totGX
    P_TG = totGT / totGX
    P_GG = totGG / totGX

    probabilities = [[P_AA, P_CA, P_GA, P_TA],[P_AC, P_CC, P_GC, P_TC],[P_GA, P_GC, P_GG, P_TG],[P_AT, P_CT, P_GT, P_TT]]
    for i in range(len(probabilities)):
        for j in range(len(probabilities[i])):
            probabilities[i][j]=round(probabilities[i][j],2) 
    matrix = numpy.array(probabilities)
    return matrix

genome=input('Insert the genome: ')
choice=input('Press 1 to obtain the prediction of your sequence being inside or outside a CpG island, or press 2 to visualize the distribution of scores according to a fixed window.')
if choice=='1':
    result = query(genome)
    if result>0:
        print('Your sequence is more likely to be inside a CpG island.')
    elif result<0:
        print('Your sequence is more likely to be outside a CpG island.')
    else:
        print('The probability of your sequence of being inside a CpG island equals the probability of being outside one of them.')
else:
    result = scan_genome(genome)
    for i in result:
        print(i)


