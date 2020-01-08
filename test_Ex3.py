import random
import Ex3
def getGenome():
    genome=''.join(random.choice('AGCT') for i in range(20))
    p1 = ''.join(random.choice('AGCT') for i in range(5))
    position1 = random.randint(0,20)
    position2 = random.randint(0,20)
    while position2 == position1:
        position2 = random.randint(0,20)
    position3 = random.randint(0,20)
    while position2 == position3 or position1 == position3:
        position3 = random.randint(0,len(20))
    p2 = ''.join(random.choice('AGCT') for i in range(5))
    while p2==p1:
        p2 = ''.join(random.choice('AGCT') for i in range(5))
    genome = genome[:position1]+p1+genome[position1:]
    genome = genome[:position2]+p2+genome[position2:]
    genome = genome[:position3]+p2+genome[position3:]
    p3 = ''.join(random.choice('AGCT') for i in range(5))
    while p3 in genome:
        p3 = ''.join(random.choice('AGCT') for i in range(5))
    return genome, p1, p2, p3

genome, p1, p2, p3 = getGenome()
print('The BWT is: ',Ex3.BWT(genome))
print('The reverse is: ',Ex3.reverse(genome))
print('Does the query string ',p1,' match?',Ex3.matching(genome, p1))
print('Does the query string ',p2,' match?',Ex3.matching(genome, p2))      
print('Does the query string ',p3,' match?',Ex3.matching(genome, p3))
