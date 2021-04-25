from Bio import SeqIO
from collections import defaultdict
from matplotlib import pyplot as plt 
from tabulate import tabulate
import gzip
import re
import time
import numpy as np

pattern_names_label = []
plot_counter = 0
numOfRecordsLen = 0
currentRecord = 0
table_entries = [['Sequence number and pattern used', 'Number of skipped alignments', 'Heuristic name']]
listOfFiles = [r"C:\Users\Aleksandar\source\repos\PythonApplication1\PythonApplication1\GCA_003713225.1_Cara_1.0_genomic.fna.gz", r"C:\Users\Aleksandar\source\repos\PythonApplication1\PythonApplication1\GCA_003957725.1_ASM395772v1_genomic.fna.gz", r"C:\Users\Aleksandar\source\repos\PythonApplication1\PythonApplication1\GCA_900095145.2_PAHARI_EIJ_v1.1_genomic.fna.gz" ]
listOfPatterns = [["ATGCATG", "TCTCTCTA", "TTCACTACTCTCA"], ["ATGATG", "CTCTCTA", "TCACTACTCTCA"], ["ACGATG", "CTCGACTA", "TCACTACTAATTCG"]]

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
    #ako hocemo da imamo prethodno slovo od poslednjeg u nizu ponavljajucih suffixa, mada je nepotrebno, posto ako na njima padne svakako postoji sledeci
    #a ako njih prodje, nema potrebe da ih vise razmatramo u okviru iste iteracije
    #bposset = set(bpos)
    #for item in bposset:
    #    u bpossetu se nalaze jedinstveni indeksi; ovo se radi da bi se dodao karakter koji se nalazi prije poslednjeg u nizu ponavljajucih suffixa, a ovo uspijeva
    #    zahvaljujuci cinjenici da svaki od suffixa u bpos nizu sadrzi indeks posljednjeg suffixa; iteracijom kroz ovaj bposset dodajemo prethodnike posljednjih suffixa
    #    if item <= m and pattern[item:]:
    #        charBeforeBorderLookup[pattern[item:]].append({'previousChar':pattern[item-1], 'index':item-1})
    #        charBeforeBorderLookup[pattern[item:]].reverse()
    #print(bposset)
        
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
                "SFGATFGACGAAACGAGTAGCSFGATAGACGA",
                "CTATCGAAGTAGCCGATTAGC",
                ""
            }

    @staticmethod
    def GetPatterns():
        return {#"A", check_first_n_chars vece od len(pattern)?
                #"AA", isto kao i gore?
                "CGA",
                "SFGATAGACGA"
            }

    @staticmethod
    def PerformTests():
        i = 1
        listOfSkippedAlignments = {"Heuristic1": [], "Heuristic2": [], "Heuristic 1 and 2": [], "Boyer Moore": []}
        for sequence in UserTests.GetSequences():
            for pattern in UserTests.GetPatterns():
                searchPattern(sequence, pattern, i, listOfSkippedAlignments, len(UserTests.GetSequences())*len(UserTests.GetPatterns()) == i)
                i+=1

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
  
    # initialize all occurrence of shift to 0
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
            #preimenovati u broj iteracija. Posto je kod naivnog algoritma broj iteracija jednak duzini teksta, broj preskocenih poredjenja je jednak razlici broja iteracija
            #naivnog algoritma i broja iteracija heuristike
            skipped_alignments += 1
            i += shift
        skipped_alignments = len(text) - skipped_alignments
        listOfSkippedAlignments["Heuristic1"].append(skipped_alignments)
        print("     Total skipped alignments by heuristic 1 is: " + str(skipped_alignments))
        return foundList

class Heuristic2:
  
    '''Search for a pattern in given text using 
    Boyer Moore algorithm with Good suffix rule '''
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
                    #uslov da je index manji od jot je VRLO BITAN. Zato sto bez tog uslova bi u patternu sa ponavljajucim sufixima, mogli da dobijemo negativan pomjeraj
                    #ako se pattern ne poklopi na slovu koje je lijevo od jednog od ponavljajucih suffixa
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
        # s is shift of the pattern with respect to text
        s = 0
        pattern_length = len(pattern)
        text_length = len(text)
        # do preprocessing
        shift, bpos = preprocess_strong_suffix(pattern, pattern_length)
        shift = preprocess_case2(shift, bpos, pattern, pattern_length)

        #heuristic optimization
        
        charBeforeBorderLookup = preprocessingForHeuristic2(bpos,pattern)
        skipped_alignments = 0
        foundList = []
        while s <= text_length - pattern_length:
            shift_amount = self.getNextShift(s, text, pattern, shift, bpos, charBeforeBorderLookup, foundList)
            skipped_alignments += 1
            s += shift_amount
        skipped_alignments = len(text) - skipped_alignments
        listOfSkippedAlignments["Heuristic2"].append(skipped_alignments)
        print("     Total skipped alignments by heuristic 2 is: " + str(skipped_alignments))
        return foundList


class Heuristic1and2:
    '''Search for a pattern in given text using 
    Boyer Moore algorithm with Good suffix rule '''

    def search(self, text, pattern, listOfSkippedAlignments):
        # s is shift of the pattern with respect to text
        s = 0
        pattern_length = len(pattern)
        text_length = len(text)
        # do preprocessing for heuristic 1
        table = preprocessingForHeuristic1(pattern)
        # end preprocessing for heuristic 1
        # do preprocessing for heuristic 2
        shift, bpos = preprocess_strong_suffix(pattern, pattern_length)
        shift = preprocess_case2(shift, bpos, pattern, pattern_length)
        charBeforeBorderLookup = preprocessingForHeuristic2(bpos,pattern)
        #finished preprocessing for heuristic 2
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
        listOfSkippedAlignments["Heuristic 1 and 2"].append(skipped_alignments)
        print("     Total skipped alignments by heuristic 1 and 2 is:  "+ str(skipped_alignments))
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
            #yield i
            foundList.append(s)
            return 1
            #i = i + 1
        else:
            bad_char_pos = -1
            if(text[s+j] in table):
                bad_char_pos = table[text[s+j]]
            return max(1, j-bad_char_pos)

    def search(self, text, pattern, listOfSkippedAlignments):
        # s is shift of the pattern with respect to text
        s = 0
        pattern_length = len(pattern)
        text_length = len(text)

        # do preprocessing
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
        listOfSkippedAlignments["Boyer Moore"].append(skipped_alignments)
        print("     Total skipped alignments by Boyer Moore full algorithm is: " + str(skipped_alignments))
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
            print("     Starting search using: '" + heuristic + "'")
            result = list(self.Heuristics[heuristic].search(text,pattern, listOfSkippedAlignments))
            print("     Number of matches: " + str(len(result)))
            #print("     Matches at: " + str(result))
            print("     Searching finished! Time spent: ", time.time() - start)
        else:
            print("     Heuristic '" + heuristic + "' is not implemented!")


def searchPattern(text, pattern, i, listOfSkippedAlignments, forcePlot):
    #print("Searching text: \"" + text + "\" for pattern: \"" + pattern + "\"")
    bm = BoyesMooreAlgorithm()
    bm.search("Heuristic1", text, pattern, listOfSkippedAlignments)
    bm.search("Heuristic2", text, pattern, listOfSkippedAlignments)
    bm.search("Heuristic1and2", text, pattern, listOfSkippedAlignments)
    bm.search("BadCharacterAndGoodSuffixRuleHeuristic",  text, pattern, listOfSkippedAlignments)
    #bm.search("NoHeuristic",  text, pattern)
    #TO DO: Add assertion of lists with pattern found indices
    pattern_names_label.append("Sequence " + str(i) +"\nPattern: " + pattern)

    showCharts(listOfSkippedAlignments, forcePlot)

numOfCharts = 1;

def showCharts(listOfSkippedAlignments, forcePlot = False):
    global plot_counter
    global numOfCharts

    if plot_counter == 2 or forcePlot:
        plt.clf()
        plt.close()
       
        x = np.arange(len(pattern_names_label))  # the label locations
        width = 0.2  # the width of the bars

        fig, ax = plt.subplots()
        numOfCharts += 1
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

        ax.bar_label(rects1, padding=3)
        ax.bar_label(rects2, padding=3)
        ax.bar_label(rects3, padding=3)
        ax.bar_label(rects4, padding=3)

        fig.tight_layout()
        plt.show()
        populate_table_entries(listOfSkippedAlignments, pattern_names_label)
        pattern_names_label.clear()
        listOfSkippedAlignments["Heuristic1"].clear()
        listOfSkippedAlignments["Heuristic2"].clear()
        listOfSkippedAlignments["Heuristic 1 and 2"].clear()
        listOfSkippedAlignments["Boyer Moore"].clear()
    plot_counter = (plot_counter + 1)%3

#UserTests.PerformTests()

for file in listOfFiles:
    listOfSkippedAlignments = {"Heuristic1": [], "Heuristic2": [], "Heuristic 1 and 2": [], "Boyer Moore": []}
    currentRecord = 0
    table_entries = [['Sequence number and pattern used', 'Number of skipped alignments', 'Heuristic name']]
    for seq in getSequencesFromFile(file):
       currentRecord += 1
       for pattern in listOfPatterns:
           searchPattern(seq, pattern, currentRecord, listOfSkippedAlignments, currentRecord == numOfRecordsLen)
    print(tabulate(table_entries, headers='firstrow', tablefmt='fancy_grid', showindex = True))
#TO DO: Ako se bude imalo vremena, izdvojiti preprocesiranje za heuristiku 2 u zasebnu funkciju