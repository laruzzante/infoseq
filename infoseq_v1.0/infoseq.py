#!/usr/bin/env python3

version = '1.0'

'''
    infoseq
    Author: Livio Ruzzante
    Date: 19th October 2018
    Version: 1.0
'''


import argparse
import sys



'''
    Parsing arguments
'''

parser = argparse.ArgumentParser(description='## Computing Sequence Summary Statistics ##\n')

parser.add_argument('inputFile', metavar='input_file', type=str,
                    help='Input file containing the FASTA sequence.')
args = parser.parse_args()



'''
    Reading input file
'''

try:
    seq = open(args.inputFile,'r')
except OSError as err:
    print("\nOS error: {0}".format(err))
    print("\nERROR: might be missing FASTA file input argument.")
    print("Usage example: < infoseq sequenceX.fa >\n")
    sys.exit()
except:
    print("Unexpected error:", sys.exc_info()[0])
    raise
    sys.exit()



'''
    Creating dictionary of sequences to scaffold IDs
'''

lines = seq.readlines()
seq2scaffDict = {}
seq = ""

for i in range(0,len(lines)):
    if lines[i].startswith('>'):
        scaff = lines[i].strip().strip('>').split(' ')[0]
    elif i < len(lines)-1 and lines[i+1].startswith('>'):
        seq = seq + lines[i].strip()
        seq2scaffDict[scaff] = seq
        seq = ""
    else:
        seq = seq + lines[i].strip()
        seq2scaffDict[scaff] = seq



'''
    Computing basic sequence info: total length [bps], number of scaffolds
'''

seqLength2scaffDict = {}
scaffLengths = []
totLength = 0
scaffoldNumber = len(seq2scaffDict.keys())
scaffoldMore1kbNumber = 0

for scaff,seq in seq2scaffDict.items():
    L = len(seq)
    seqLength2scaffDict[scaff] = L
    scaffLengths.append(L)
    totLength += L
    if L > 1000:
        scaffoldMore1kbNumber += 1



'''
    Computing N50 length [bps]
'''

sortedScaffLengths = sorted(scaffLengths,reverse=True)
halfLength = totLength/2
cumsum = []
i = 0

for l in sortedScaffLengths:
    cumsum.append(l)
    while cumsum[i] < halfLength and i < len(sortedScaffLengths):
        i += 1
        cumsum.append(sum(sortedScaffLengths[0:i]))

    N50 = sortedScaffLengths[i-1]



'''
    Printing output to terminal
'''


print('''            ___''')
print(''' _        /'___)''')
print('''(_)  ___  | (__   _     ___    __     _ _''')
print('''| |/' _ `\| ,__)/'_`\ /',__) /'__`\ /'_` )''')
print('''| || ( ) || |  ( (_) )\__, \(  ___/( (_) |''')
print('''(_)(_) (_)(_)  `\___/'(____/`\____)`\__, |''')
print('''                                       | |''')
print('''                                       (_)''')
print('v'+version)

print('\nInput Sequence File: ', args.inputFile)
print('Total Length [bps]: ', totLength)
print('Number of Scaffolds: ', scaffoldNumber)
print('Number of Scaffolds >1kb: ', scaffoldMore1kbNumber)
print('Longest Scaffold [bps]: ', sortedScaffLengths[0])
print('N50 [bps]: ', N50,'\n')
