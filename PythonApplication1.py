from Bio import SeqIO
from collections import defaultdict
from matplotlib import pyplot as plt 
from tabulate import tabulate
import gzip
import re
import time
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import io

from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

pattern_names_label = []
numOfRecordsLen = 0
table_entries = [['Text number and pattern used', 'Number of skipped alignments', 'Heuristic name']]
skipped_alignments_by_pattern = {"Heuristic1": 0, "Heuristic2" : 0, "Heuristic 1 and 2": 0, "Boyer Moore": 0}


def populate_table_entries(listOfSkippedAlignments, pattern_names_label):
    global table_entries
    number_of_items = len(listOfSkippedAlignments['Heuristic1'])
    for i in range(number_of_items):
        table_entries.append([pattern_names_label[i], listOfSkippedAlignments['Heuristic1'][i], 'Heuristic 1'])
        table_entries.append([pattern_names_label[i], listOfSkippedAlignments['Heuristic2'][i], 'Heuristic 2'])
        table_entries.append([pattern_names_label[i], listOfSkippedAlignments['Heuristic 1 and 2'][i], 'Heuristic 1 and 2'])
        table_entries.append([pattern_names_label[i], listOfSkippedAlignments['Boyer Moore'][i], 'Boyer Moore'])
        
        
def preprocessingForHeuristic2(bpos,pattern):
    charBeforeBorderLookup = defaultdict(list)
    i = 1
    while i < len(bpos)-1:
        if pattern[bpos[i]:]:
            charBeforeBorderLookup[pattern[bpos[i]:]].append({'previousChar':pattern[i-1], 'index':i-1})
        i+=1

    for key in charBeforeBorderLookup:
        charBeforeBorderLookup[key].reverse()
    return charBeforeBorderLookup

def preprocessingForHeuristic1(pattern):
    charPositionTable = defaultdict(list)
    i = len(pattern) - 1
    while(i >= 0):
        charPositionTable[pattern[i]].append(i)
        i-=1
    return charPositionTable


def getSequencesFromFile(file):
    global numOfRecordsLen
    with gzip.open(file, "rt") as handle:
        list_of_seq = list(SeqIO.parse(handle, "fasta"))
        numOfRecordsLen = len(list_of_seq)
        for record in list_of_seq:
            yield str(record.seq).upper()

class UserTests:
    @staticmethod
    def GetSequences():
        return {"AAAAAAAAAAAAAAAA",
                "TFGACGAAACGAGTAGCSFGATAGACGA",
                "CTATCGAAGTAGCCGATTAGC",
                "AAA"
                "TCGATGCG",
                "AACCAACCAC",
                "ACACACACACACA",
                "ACCAACCA",
                "ACACCACCAACA",
                "TCAGCGCGCTAGCGACTCGCTCAAGCATCGATCGACTGATCGGCCAACGCGAGCGACG",
                "TCGCGCTAGCATCGATCGATCGTAGCA",
                "AGCGCGAGCATAGCGCATACGTACG",
                "TCGAGCTGCTAGCACGGCATGACTATCGCA",
                "CAGGCTTAGCTGACTAT", 
                "TGCATATCGATCTGAAAGCGCAGTGCATACGTCAG",
                "GCATGACTGATCGCATGCTGAC"
            }

    @staticmethod
    def GetPatterns():
        return {
                "CGA",
                "GATAGACGA",
                "AACCACCAC",
                "CABADABA",
                "ACTACTAC",
                "GCTG",
                "A",
                "AA",
                "ACA",
                "ACCA",
                "ACACA",
                "CA",
                "CCA",
                "CTACTA",
                "GAC",
                "CAG",
                "TCG",
                "AGCT"
            }

    @staticmethod
    def PerformTests():
        global pdf
        global table_entries
        sequence_number = 1
        for sequence in UserTests.GetSequences():
            print("Searching sequence: " + sequence)
            start = time.time()
            listOfSkippedAlignments = {"Heuristic1": [], "Heuristic2": [], "Heuristic 1 and 2": [], "Boyer Moore": []}
            table_entries = [['Sequence number and pattern used', 'Number of skipped alignments', 'Heuristic name']]
            pattern_names_label.clear()
            with PdfPages("plots_" + str(sequence_number) + ".pdf") as pdf:
                for pattern in UserTests.GetPatterns():
                    pattern_names_label.append("Sequence: " + sequence + '\n' + "Pattern: " +pattern)
                    searchPattern(sequence, pattern, listOfSkippedAlignments)
                showCharts(listOfSkippedAlignments)
                tableText = tabulate(table_entries, headers='firstrow', tablefmt='fancy_grid', showindex = True)
                with io.open(r"table_" + str(sequence_number) + ".txt", 'w', encoding='utf-8') as tableFile:
                    tableFile.write(tableText)
            print("Searching for sequence: "+str(sequence_number)+" finished! Time spent: ", time.time() - start)
            sequence_number+=1

def preprocess_bad_character(pattern):
    bad_char_table = {}
    for i, c in  enumerate(pattern):
        bad_char_table[c] = i
    return bad_char_table

def preprocess_strong_suffix(pattern, m):
  
    # m is the length of pattern
    i = m
    j = m + 1

    bpos = [0] * (m + 1)
  
    # Initialize all occurrence of shift to 0
    shift = [0] * (m + 1)

    bpos[i] = j
  
    while i > 0:
          
        '''if character at position i-1 is 
        not equivalent to character at j-1, 
        then continue searching to right 
        of the pattern for border '''
        while j <= m and pattern[i - 1] != pattern[j - 1]:
              
            ''' the character preceding the occurrence 
            of t in pattern P is different than the 
            mismatching character in P, we stop skipping
            the occurrences and shift the pattern 
            from i to j '''
            if shift[j] == 0:
                shift[j] = j - i
  
            # Update the position of next border
            j = bpos[j]
              
        ''' p[i-1] matched with p[j-1], border is found. 
        store the beginning position of border '''
        i -= 1
        j -= 1
        bpos[i] = j

    return shift, bpos

# Preprocessing for case 2
def preprocess_case2(shift, bpos, pattern, m):
    j = bpos[0]
    for i in range(m + 1):
          
        ''' set the border position of the first character 
        of the pattern to all indices in array shift
        having shift[i] = 0 '''
        if shift[i] == 0:
            shift[i] = j
              
        ''' suffix becomes shorter than bpos[0], 
        use the position of next widest border
        as value of j '''
        if i == j:
            j = bpos[j]

    return shift;

class Heuristic1:
    
    def get_last_occurence(self, character, table):
        if(character in table):
            return table[character]
        else:
            return -1
    
    def getNextShift(self, i, text, pattern, table, foundList, text_length, pattern_length):
        current_text_index = i
        j = pattern_length - 1
        while(j>=0 and pattern[j] == text[i+j]):
            j=j-1
        if(j<0):
            #yield i
            foundList.append(i)
            return 1
            #i = i + 1
        else:
            if(table[text[i+j]]):
                for item in table[text[i+j]]:
                    if item < j:
                        return j - item
                return j + 1
            else:
               return j + 1

    def search(self, text, pattern, listOfSkippedAlignments):
        global skipped_alignments_by_pattern
        table = preprocessingForHeuristic1(pattern)
        text_length = len(text)
        pattern_len = len(pattern)
        i = 0
        skipped_alignments = 0
        foundList = []
        text_length = len(text)
        pattern_length = len(pattern)
        while(i <= text_length - pattern_length):
            shift = self.getNextShift(i, text, pattern, table, foundList, text_length, pattern_length)
            skipped_alignments += 1
            i += shift
        skipped_alignments = len(text) - skipped_alignments
        #listOfSkippedAlignments["Heuristic1"].append(skipped_alignments)
        skipped_alignments_by_pattern["Heuristic1"] += skipped_alignments

        #print("     Total skipped alignments by heuristic 1 is: " + str(skipped_alignments))
        return foundList

class Heuristic2:
  
    def getNextShift(self, s, text, pattern, shift, bpos, charBeforeBorderLookup, foundList):
  
        j = len(pattern) - 1
          
        ''' Keep reducing index j of pattern while characters of 
            pattern and text are matching at this shift s'''
        while j >= 0 and pattern[j] == text[s + j]:
            j -= 1
              
        ''' If the pattern is present at the current shift, 
            then index j will become -1 after the above loop '''
        if j < 0:
            foundList.append(s)
            return shift[0]
        else:
            '''pattern[i] != pattern[s+j] so shift the pattern 
            shift[j+1] times '''
                
            foundCharMatch = False
            foundHashLookup = False
            if pattern[bpos[j]:]:
                listOfPrevChr = charBeforeBorderLookup[pattern[bpos[j]:]]
                for item in listOfPrevChr:
                    foundHashLookup = True
                    # The condition item['index'] < j is very important because if we didn't do that check, we could get a negative shift. That could happen if the
                    # same suffix is found in the matched part of the pattern. In other words, if the mismatched character is to the left of the suffix in the pattern.
                    if (item['previousChar'] == text[s+j-1]) and (item['index'] < j):
                        return j - item['index']
                        foundCharMatch = True
                        break

            if not foundHashLookup:
               return shift[j + 1]
            else:
                if not foundCharMatch:
                   return shift[0]

    def search(self, text, pattern, listOfSkippedAlignments):
        global skipped_alignments_by_pattern
        # s is shift of the pattern with respect to text
        s = 0
        pattern_length = len(pattern)
        text_length = len(text)
        # Do the preprocessing
        shift, bpos = preprocess_strong_suffix(pattern, pattern_length)
        shift = preprocess_case2(shift, bpos, pattern, pattern_length)

        # Heuristic optimization
        
        charBeforeBorderLookup = preprocessingForHeuristic2(bpos,pattern)
        skipped_alignments = 0
        foundList = []
        while s <= text_length - pattern_length:
            shift_amount = self.getNextShift(s, text, pattern, shift, bpos, charBeforeBorderLookup, foundList)
            skipped_alignments += 1
            s += shift_amount
        skipped_alignments = len(text) - skipped_alignments
        #listOfSkippedAlignments["Heuristic2"].append(skipped_alignments)
        skipped_alignments_by_pattern["Heuristic2"] += skipped_alignments
        #print("     Total skipped alignments by heuristic 2 is: " + str(skipped_alignments))
        return foundList


class Heuristic1and2:
    
    def search(self, text, pattern, listOfSkippedAlignments):
        global skipped_alignments_by_pattern
        # s is shift of the pattern with respect to text
        s = 0
        pattern_length = len(pattern)
        text_length = len(text)
        # Do preprocessing for heuristic 1
        table = preprocessingForHeuristic1(pattern)
        # End preprocessing for heuristic 1
        # Do preprocessing for heuristic 2
        shift, bpos = preprocess_strong_suffix(pattern, pattern_length)
        shift = preprocess_case2(shift, bpos, pattern, pattern_length)
        charBeforeBorderLookup = preprocessingForHeuristic2(bpos,pattern)
        # Finished preprocessing for heuristic 2
        heuristic1 = Heuristic1()
        heuristic2 = Heuristic2()
        skipped_alignments = 0
        foundListSuffix = []
        foundListCharacter = []
        while s <= text_length - pattern_length:
            shift_amount = max(heuristic2.getNextShift(s, text, pattern, shift, bpos, charBeforeBorderLookup, foundListSuffix), heuristic1.getNextShift(s, text, pattern, table, foundListCharacter, text_length, pattern_length))
            skipped_alignments += 1
            s += shift_amount
        skipped_alignments = len(text) - skipped_alignments
        #listOfSkippedAlignments["Heuristic 1 and 2"].append(skipped_alignments)
        skipped_alignments_by_pattern["Heuristic 1 and 2"] += skipped_alignments
        #print("     Total skipped alignments by heuristic 1 and 2 is:  "+ str(skipped_alignments))
        return foundListSuffix

class BadCharacterAndGoodSuffixRuleHeuristic:
    def getNextShiftGoodSuffix(self, s, text, pattern, shift, bpos, foundList):
  
        j = len(pattern) - 1
          
        ''' Keep reducing index j of pattern while characters of 
            pattern and text are matching at this shift s'''
        while j >= 0 and pattern[j] == text[s + j]:
            j -= 1
              
        ''' If the pattern is present at the current shift, 
            then index j will become -1 after the above loop '''
        if j < 0:
            foundList.append(s)
            return shift[0]
        else:
            '''pattern[i] != pattern[s+j] so shift the pattern 
            shift[j+1] times '''
            return shift[j+1]
    def getNextShiftBadChar(self, s, text, pattern, table, foundList, pattern_length ):
        j = pattern_length - 1
        while(j>=0 and pattern[j] == text[s+j]):
            j=j-1
        if(j<0):
            foundList.append(s)
            return 1
        else:
            bad_char_pos = -1
            if(text[s+j] in table):
                bad_char_pos = table[text[s+j]]
            return max(1, j-bad_char_pos)

    def search(self, text, pattern, listOfSkippedAlignments):
        global skipped_alignments_by_pattern
        # s is shift of the pattern with respect to text
        s = 0
        pattern_length = len(pattern)
        text_length = len(text)

        # Do preprocessing
        shift, bpos = preprocess_strong_suffix(pattern, pattern_length)
        shift = preprocess_case2(shift, bpos, pattern, pattern_length)
        bad_char_table = preprocess_bad_character(pattern)
        skipped_alignments = 0
        foundListSuffix = []
        foundListCharacter = []
        while s <= text_length - pattern_length:
            shift_amount = max(self.getNextShiftGoodSuffix(s, text, pattern, shift, bpos, foundListSuffix), self.getNextShiftBadChar(s, text, pattern, bad_char_table, foundListCharacter, pattern_length))
            skipped_alignments += 1
            s += shift_amount
        skipped_alignments = len(text) - skipped_alignments
        #listOfSkippedAlignments["Boyer Moore"].append(skipped_alignments)
        skipped_alignments_by_pattern["Boyer Moore"] += skipped_alignments
        #print("     Total skipped alignments by Boyer Moore full algorithm is: " + str(skipped_alignments))
        return foundListSuffix

class NoHeuristic:
    def search(self, text, pattern):
        return [a.start() for a in list(re.finditer("(?=" + pattern + ")", text))]
    
class BoyesMooreAlgorithm:
    Heuristics = {"Heuristic1":Heuristic1(),
                "Heuristic2":Heuristic2(),
                "Heuristic1and2":Heuristic1and2(),

                "BadCharacterAndGoodSuffixRuleHeuristic":BadCharacterAndGoodSuffixRuleHeuristic(),
                 "NoHeuristic":NoHeuristic()}
    
    def search(self, heuristic, text, pattern, listOfSkippedAlignments):
        if (heuristic in self.Heuristics):
            start = time.time()
            #print("     Starting search using: '" + heuristic + "'")
            result = list(self.Heuristics[heuristic].search(text,pattern, listOfSkippedAlignments))
            #print("     Number of matches: " + str(len(result)))
            #print("     Matches at: " + str(result))
            #print("     Searching finished! Time spent: ", time.time() - start)
        

def searchPattern(text, pattern, listOfSkippedAlignments):
    #print("Searching text: \"" + text + "\" for pattern: \"" + pattern + "\"")
    bm = BoyesMooreAlgorithm()
    bm.search("Heuristic1", text, pattern, listOfSkippedAlignments)
    bm.search("Heuristic2", text, pattern, listOfSkippedAlignments)
    bm.search("Heuristic1and2", text, pattern, listOfSkippedAlignments)
    bm.search("BadCharacterAndGoodSuffixRuleHeuristic",  text, pattern, listOfSkippedAlignments)
    #bm.search("NoHeuristic",  text, pattern)
    #TO DO: Add assertion of lists with pattern found indices

def showCharts(listOfSkippedAlignments):
    global pattern_names_label
    
    # The label locations
    x = np.arange(len(pattern_names_label))
    # The width of the bars
    width = 0.2

    fig, ax = plt.subplots()
    rects1 = ax.bar(x, listOfSkippedAlignments["Heuristic1"], width, label='Heuristic1')
    rects2 = ax.bar(x + width, listOfSkippedAlignments["Heuristic2"], width, label='Heuristic2')
    rects3 = ax.bar(x + 2*width, listOfSkippedAlignments["Heuristic 1 and 2"], width, label='Heuristic 1 and 2')
    rects4 = ax.bar(x + 3*width, listOfSkippedAlignments["Boyer Moore"], width, label='Boyer Moore')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Skipped alignments')
    ax.set_title('Skipped alignments by heuristic and pattern')
    ax.set_xticks([x + 1.5*width for x in range(len(listOfSkippedAlignments["Heuristic1"]))])
    ax.set_xticklabels(pattern_names_label)
    ax.legend()
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    ax.bar_label(rects1, padding=3)
    ax.bar_label(rects2, padding=3)
    ax.bar_label(rects3, padding=3)
    ax.bar_label(rects4, padding=3)

    fig.tight_layout()     
    fig.set_size_inches(9.69, 6.27)
    pdf.savefig(fig, dpi = 1200)
    plt.show()

    populate_table_entries(listOfSkippedAlignments, pattern_names_label)
    pattern_names_label.clear()
    listOfSkippedAlignments["Heuristic1"].clear()
    listOfSkippedAlignments["Heuristic2"].clear()
    listOfSkippedAlignments["Heuristic 1 and 2"].clear()
    listOfSkippedAlignments["Boyer Moore"].clear()

def searchGenomesFromFiles():
    global listOfSkippedAlignments
    global table_entries
    global pattern_names_label
    global pdf
    global skipped_alignments_by_pattern
    numFile = 0
    listOfFiles = [r"chr1c.fna.gz", r"chrX.fna.gz", r"GCA_002588565.1_ASM258856v1_genomic.fna.gz" ]
    fileNames = ["Coffea arabica chr 1C", "Mus pahari chr X", "Coelastrella genome"]
    listOfPatterns = [["ATGCATG", "TCTCTCTA", "TTCACTACTCTCA"], ["ATGATG", "CTCTCTA", "TCACTACTCTCA"], ["ACGATGAGCGTGCTCCGCA", "CTCGACTAACGCCTA", "TATCCGGCGACGTCGGT"]]
    for file in listOfFiles:
        start = time.time()
        print("Searching file: " + file)
        listOfSkippedAlignments = {"Heuristic1": [], "Heuristic2": [], "Heuristic 1 and 2": [], "Boyer Moore": []}
        table_entries = [['Text number and pattern used', 'Number of skipped alignments', 'Heuristic name']]
        with PdfPages("plots_" + file + ".pdf") as pdf:
            for pattern in listOfPatterns[numFile]:
                pattern_names_label.append("Genome and chromosome :\n" + fileNames[numFile] + '\n' + "Pattern: " + pattern)
                for seq in getSequencesFromFile(file):
                    searchPattern(seq, pattern, listOfSkippedAlignments)
                listOfSkippedAlignments["Heuristic1"].append(skipped_alignments_by_pattern["Heuristic1"])
                listOfSkippedAlignments["Heuristic2"].append(skipped_alignments_by_pattern["Heuristic2"])
                listOfSkippedAlignments["Heuristic 1 and 2"].append(skipped_alignments_by_pattern["Heuristic 1 and 2"])
                listOfSkippedAlignments["Boyer Moore"].append(skipped_alignments_by_pattern["Boyer Moore"])
                skipped_alignments_by_pattern = {"Heuristic1": 0, "Heuristic2" : 0, "Heuristic 1 and 2": 0, "Boyer Moore": 0}
            showCharts(listOfSkippedAlignments)
            tableText = tabulate(table_entries, headers='firstrow', tablefmt='fancy_grid', showindex = True)
            with io.open(r"table_" + file + ".txt", 'w', encoding='utf-8') as tableFile:
                tableFile.write(tableText)
        print("Searching in file: "+file+" finished! Time spent: ", time.time() - start)
        numFile += 1

searchGenomesFromFiles()
#UserTests.PerformTests()