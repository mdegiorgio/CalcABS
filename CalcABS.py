#This script takes the parsed allele count table as input, and output the ABS scores
from math import log, floor
import sys, optparse

#Calculating pair-wise variation between 2 populations for ABS. 
#a, b are the indeces of populations. n is simulation number, g is neutral or selected
#This function will return/write theta values for all the positions betwenn the 2 pop
def variation(l, a,b, C_x1): #0rep 1pos 2x1 3n1
    i = [a,b]
    l = l.strip().split('\t')
    total = [l[2*x + C_x1 + 1] for x in i]
    drv = [l[2*x + C_x1] for x in i]
    try:
        freq = [float(drv[x])/float(total[x]) for x in [0,1]]
        n = [float(x)/2 for x in total] # number of diploid individuals in the sample
    except(ValueError):
        print('ValueError: Please make sure your count data are all numbers.')
    het = [2*x*(1-x) for x in freq]
    #between-pop var.
    a_l = (freq[0]-freq[1])**2 - (n[0]*het[0]+n[1]*het[1])*sum(n)/(4*n[0]*n[1]*(sum(n)-1))
    #total var.
    al_bl = (freq[0]-freq[1])**2 + (n[0]*het[0]+n[1]*het[1])*(4*n[0]*n[1]-sum(n))/(4*n[0]*n[1]*(sum(n)-1))
    return(a_l, al_bl)


#Improved function to calculate all 6 pairs of pop in a line, and return a dictionary (dashed-table)
def varPairs(l, C_x1): #l: 0chr 1pos 2x1 n1 x2 n2 x3 n3 x4 n4
    Pairs = {}
    pairs = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
    for (a,b) in pairs:
        (al, albl) = variation(l, a,b, C_x1)
        Pairs[(a,b)] = (al, albl)
    return Pairs

#Calculate FST on a sliding window basis
#windict is the dictionary
def sumAB(windict):
    #Initialization
    A = {}; AB = {}
    pairs = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
    for p in pairs:
        A[p] = 0
        AB[p] = 0
    for pos in list(windict.keys()):
        for p in pairs:
            (al,albl) = windict[pos][p]
            if al != 'NA' and albl != 'NA':
                A[p] += al
                AB[p] += albl
    return (A, AB)

def FstToABS(A, AB): #inputs are dictionaries
    Dxy = {}
    pairs = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
    XYZ = ['12|34','13|24','14|23']
    for p in pairs:
        Fst = 0
        if A[p] > 0 and AB[p] > 0:
            Fst = A[p]/AB[p]
        if Fst >= 1:
                print(p, Fst)
                Fst = 0.9
        Dxy[p] = -1*log(1-Fst)
    X = Dxy[(0,1)]+Dxy[(2,3)]
    Y = Dxy[(0,2)]+Dxy[(1,3)]
    Z = Dxy[(0,3)]+Dxy[(1,2)]
    ABS = 0.5*(max(X,Y,Z)-min(X,Y,Z))
    MinPair = XYZ[[X,Y,Z].index(min(X,Y,Z))]
    return (ABS, MinPair)

#Topology is fixed to be one of 12|34, 13|24, and 14|23
def Fst2ABS_topoFix(A, AB, topo):
    Dxy = {}
    pairs = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
    XYZ = ['12|34','13|24','14|23']
    for p in pairs:
        Fst = 0
        if A[p] > 0 and AB[p] > 0:
            Fst = A[p]/AB[p]
        if Fst >= 1:
                print(p, Fst)
                Fst = 0.9
        Dxy[p] = -1*log(1-Fst)
    #Summing up distinct pairs
    X = Dxy[(0,1)]+Dxy[(2,3)]
    Y = Dxy[(0,2)]+Dxy[(1,3)]
    Z = Dxy[(0,3)]+Dxy[(1,2)]
    #Fixing the given topology
    Dlist = [X,Y,Z]
    M = Dlist[XYZ.index(topo)] #the minimum pair
    Dlist.remove(M) #remove the minimum pair from the three pairs
    ABS = 0.5*(max(Dlist) - M)
    return ABS

#Given a string of the following format, return the position of this line
#l: 0rep 1pos 2x1 n1 x2 n2 x3 n3 x4 n4
def getPos(line, C_pos):
    l = line.strip().split('\t')
    pos = float(l[C_pos]) #output in MS is non-integar
    return pos

def getChr(line, C_chr):
    l = line.strip().split('\t')
    ch = l[C_chr]
    return ch

#'infile' and 'outfile' are file directories of input and output
#construct a list of library: {pos:{(pair):(A,AB)}}
def slidingWin(infile, outfile, win, step, C_ch, C_pos, C_x1, topo = ''):
    score = open(outfile,'w')
    window = {}
    count = open(infile, 'r')
    l = count.readline() #skip the header
    l = count.readline() #first line
    window[getPos(l, C_pos)] = varPairs(l,C_x1) #Storing al & albl
    #print(getPos(l, C_pos))
    l = l.strip().split('\t')
    if C_ch != '':
        score.write('chr\tmidPos\tscore\tnumSites\ttopology\n')
        ch = l[C_ch] #keep it as a string to allow for characters.
    else:
        score.write('midPos\tscore\tnumSites\ttopology\n')
    #Settle the window parameter
    start = float(l[C_pos])
    start = int(floor(start/(win/2))*win/2) #give it more zero's...
    end = int(start + win)
    mid = int(start + win/2)
    n_step = 0
    n_line = 0
    n_stride = 0
    NewWin = {}
    LastWin = {}
    #Initializing numerator and denominator
    A = {}; AB = {}
    pairs = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
    for p in pairs:
        A[p] = 0
        AB[p] = 0
    #After all initial window parameters are set
    while l != '':
        l = count.readline()
        n_line += 1
        if l =='': 
            #print 'EOF. step:', n_step, mid
            continue
        while getPos(l, C_pos) <= end:
            window[getPos(l, C_pos)] = varPairs(l, C_x1)
            l = count.readline()
            n_line += 1
            if l == '':
                break #break out from current while loop
        else:
            NewWin[getPos(l, C_pos)] = varPairs(l, C_x1) #grab this line into the next window
        n_step += 1
        #in case the window still contain sites it shouldn't grab
        currentList = [pos for pos in list(window.keys()) if (pos <= end)&(pos > start)]
        if list(window.keys()) != currentList:
            Extra = dict((k,v) for (k,v) in list(window.items()) if k > end)
            #Syntax in python2.7+: Extra = {k:v for (k,v) in window if k > end}
            NewWin.update(Extra)
            window = dict((k,v) for (k,v) in list(window.items()) if k in currentList)
            #window = {k:v for (k,v) in window if k in currentList}
        if len(window) == 0:
            s = 'NA'
        else:
            (A, AB) = sumAB(window)
            if topo != '':
                MP = topo
                s = Fst2ABS_topoFix(A,AB, MP)
            else:
                (s, MP) = FstToABS(A, AB)
        if C_ch != '':
            ol = '%s\t%s\t%s\t%s\t%s\n' % (ch, mid,s,len(window), MP)
        else:
            ol = '%s\t%s\t%s\t%s\n' % (mid,s,len(window), MP)
        score.write(ol) #write the output line
        #Take a new step, reset window parameters and delete old lines
        #in case this extra line (the current line) is many steps away:
        if l == '':
            break
        if getPos(l, C_pos) > end+step:
            stride = floor((getPos(l, C_pos)-end)/step) #the number of steps to skip
            n_step += int(stride)
            n_stride += 1
            start = int(start + step * stride)
            end = int(end + step*stride)
            mid = int(start + win/2)
        else:
            start = start + step
            end = int(end + step)
            mid = int(start + win/2)
        #print start, end
        LastWin = window
        NewWin.update(window)
        window = NewWin
        NewWin = {}

def str_to_float(l):
    for i in range(len(l)):
        if l[i] != 'NA':
            l[i] = float(l[i])
    return l

#Parse the index
def parse_index(indices):
    indices = indices.split(',')
    try:
        C_ch = int(indices[0]) -1
    except(ValueError):
        C_ch = indices[0]
    C_pos = int(indices[1]) -1
    C_x1 = int(indices[2]) -1
    return(C_ch, C_pos, C_x1)

#Things to check: 1. if chromosome name stays the same; 2. if the file contain header, 3. if the file is tab delimited, 4. if data is complete, 5. if there're monomorphic sites
def error_check(input_file, indices):
    #Checking c,p,x
    if indices.count(',') != 2:
        print('Error. Please indicate column indices in the correct format: c,p,x')
        sys.exit()
    else:
        (C_ch, C_pos, C_x1) = parse_index(indices)
    #Checking input file
    Nline = 1
    ch = ''
    with open(input_file) as f:
        l = next(f) #header
        colnum = len(l.strip().split('\t'))
        if colnum < 9:
            print('Error: input format incorrect. Please make sure your file contains all necessary info and is tab-delimited.')
            sys.exit()
        for l in f:
            Nline += 1
            l = l.strip().split('\t')
            if len(l) != colnum:
                print('Error: please make sure every row has same number of columns. The file should be tab-delimited.')
            try: 
                float(l[C_pos])
            except(ValueError):
                print('line %s: %s' % (Nline, ' '.join(l)))
                print('Error: make sure positions are represented in numbers.')
                sys.exit()
            if Nline >2:
                if float(l[C_pos]) <= pos:
                    print('Error: make sure the positions are listed in ascending order.')
                    sys.exit()
            pos = float(l[C_pos])
            if ch != '':
                if l[C_ch] != ch:
                    print(ch, l[C_ch])
                    print('Error: Chromosome name inconsistent. Please make sure this file only contains data from one chromosome.')
                    sys.exit()
            if C_ch != '':
                ch = l[C_ch]
            try:
                count = [float(l[2*x+C_x1]) for x in range(4)]
                samsize = [float(l[2*x+C_x1+1]) for x in range(4)]
            except(ValueError):
                print('line %s: %s' % (Nline, ' '.join(l)))
                print('Error: allele counts can only be numbers.')
                sys.exit()
            if sum(count) == 0 or sum(count) == sum(samsize):
                print('line %s: %s' % (Nline, ' '.join(l)))
                print('Error: please check your input for <indice> and do not include monomorphic sites.')
                sys.exit()
            elif sum(count) > sum(samsize):
                print('line %s: %s' % (Nline, ' '.join(l)))
                print('Error: please check your input for <indice> and the allele count cannot exceed the total number of observed alleles.')
                sys.exit()
    print('Error check finished.')

def main():
    #================================================================
    #Parsing arguments
    parser = optparse.OptionParser(usage='%prog -i <input file> -o <output file> -n <c,p,x> -w <window size> -s <step size> [--check] [--fix <pair>]')
    parser.add_option('-i','--input', dest='infile', help = 'Path and name of your input file.')
    parser.add_option('-o','--output', dest='outfile', help = 'Path and name of your output file.')
    parser.add_option('-n','--indices', dest='index', help = 'Index numbers for the columns of chromosome name c, locus positions p, and allele counts of the first population x, respectively. Format of this argument should be \'c,p,x\' or \',p,x\', without space.\n')
    parser.add_option('-w','--window', dest='win', type = 'float', help='Length of the sliding window, in the number of nucleotides (nt), for ABS scan.\n')
    parser.add_option('-s','--step', dest='step', type = 'float', help='Length of each step the sliding window takes, in the number of nucleotides (nt), while scanning the data with ABS.\n')
    parser.add_option('--check', action='store_true', default=False, help='Option to check the file format.')
    parser.add_option('--fix', action = 'store', dest='pair', type = 'string', help='Option to fix the topology formed by population 1, 2, 3 and 4. Argument can be \'12\', \'13\', \'14\', \'23\', \'24\', or \'34\', identifying the pair of target populations.\n')

    #====================================
    opt,v = parser.parse_args(sys.argv[1:])

    (C_ch, C_pos, C_x1) = parse_index(opt.index)

    if opt.check == True:
        if opt.infile and opt.index:
            print('Checking input file format...')
            error_check(opt.infile, opt.index)
            print('Format correct. Continue to compute ABS scores...\n')
        else:
            print('Please make sure to provide input file and column index.')
            sys.exit()
    elif opt.infile and opt.outfile and opt.index and opt.win and opt.step:
        print('Calculating ABS scores based on input file: %s.\nOutput written to %s.' % (opt.infile, opt.outfile))
        print('Sliding window of size %s nt, taking %s nt step.' % (opt.win, opt.step))
        print('column No. %s is positions. Allele count data start at column %s' % (C_pos + 1, C_x1 + 1))
        if opt.pair:
            print('Topology fixed. Target sister populations are', opt.pair)
            if opt.pair in ['12','34']:
                topo = '12|34'
            elif opt.pair in ['13','24']:
                topo = '13|24'
            elif opt.pair in ['14','23']:
                topo = '14|23'
            else:
                print('Error: please provide the correct argument for target populaiton pair. Use \'-h\' or \'--help\' for more information.')
        else:
            topo = ''
        slidingWin(opt.infile, opt.outfile, opt.win, opt.step, C_ch, C_pos, C_x1, topo = topo)
    else:
        print('To initiate the scan, please make sure to include input file (-i), indices of critical columns (-n), output path (-o), scanning window size (-w) and step size (-s).')

    print('Scan finished.')



if __name__ == '__main__':
    main()
