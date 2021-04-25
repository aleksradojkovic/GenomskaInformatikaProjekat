from matplotlib import pyplot as plt 
import time
import re
from collections import defaultdict
import numpy as np

listOfSkippedalignments = defaultdict(list)
listOfSkippedAlignments = {"Heuristic1": [], "Heuristic2": [], "Heuristic 1 and 2": [], "Boyer Moore": []}
pattern_names_label = []


def preprocess_bad_character(pattern):
    bad_char_table = {}
    for i, c in  enumerate(pattern):
        bad_char_table[c] = i
    return bad_char_table

def preprocess_strong_suffix(pat, m):
  
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
        while j <= m and pat[i - 1] != pat[j - 1]:
              
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
def preprocess_case2(shift, bpos, pat, m):
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
    def check_first_n_characters(self, pattern, text, n, j):
        i = 0
        for i in range(n):
            if(text[j+i] == pattern[i]):
                break
        return i+1
    
    def get_last_occurence(self, character, table):
        if(character in table):
            return table[character]
        else:
            return -1
    
        #s (i) u ovom kontekstu dokle smo stigli u tekstu, indeks kursora na tekst (odkale poceti search za oviu iteraciju)
        #ova funkcija vraca pomjeraj za heuristiku 1. Kad  je chars_to_skip jednak jedinici, ova heuristika postaje bad character rule
    def getNextShift(self, i, text, pattern, chars_to_skip, table, foundList):
        current_text_index = i
        j = len(pattern) - 1
        while(j>=0 and pattern[j] == text[i+j]):
            j=j-1
        if(j<0):
            #yield i
            foundList.append(i)
            return 1
            #i = i + 1
        else:
            current_text_index = current_text_index + len(pattern) - j
            #i = i + max(self.check_first_n_characters(pattern, text, chars_to_skip, current_text_index), j-self.get_last_occurence(text[i+j],table))
            return max(self.check_first_n_characters(pattern, text, chars_to_skip, current_text_index), j-self.get_last_occurence(text[i+j],table))

    def search(self, text, pattern, chars_to_skip):
        table = preprocess_bad_character(pattern)
        i = 0
        skipped_alignments = 0
        foundList = []
        while(i<=len(text) - len(pattern)):
            shift = self.getNextShift(i, text, pattern, chars_to_skip, table, foundList)
            #preimenovati u broj iteracija. Posto je kod naivnog algoritma broj iteracija jednak duzini teksta, broj preskocenih poredjenja je jednak razlici broja iteracija
            #naivnog algoritma i broja iteracija heuristike
            skipped_alignments += 1
            i += shift
        skipped_alignments = len(text) - skipped_alignments
        listOfSkippedAlignments["Heuristic1"].append(skipped_alignments)
        print("     Total skipped alignments by heuristic 1 with n = " + str(chars_to_skip) +" is: " + str(skipped_alignments))
        return foundList

class Heuristic2:
  
    '''Search for a pattern in given text using 
    Boyer Moore algorithm with Good suffix rule '''
    def getNextShift(self, s, text, pat, shift, bpos, charBeforeBorderLookup, foundList):
  
        j = len(pat) - 1
          
        ''' Keep reducing index j of pattern while characters of 
            pattern and text are matching at this shift s'''
        while j >= 0 and pat[j] == text[s + j]:
            j -= 1
              
        ''' If the pattern is present at the current shift, 
            then index j will become -1 after the above loop '''
        if j < 0:
            foundList.append(s)
            return shift[0]
        else:
            '''pat[i] != pat[s+j] so shift the pattern 
            shift[j+1] times '''
                
            foundCharMatch = False
            foundHashLookup = False
            if pat[bpos[j]:]:
                listOfPrevChr = charBeforeBorderLookup[pat[bpos[j]:]]
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
                   return len(pat)

    def search(self, text, pat, n):
        # s is shift of the pattern with respect to text
        s = 0
        m = len(pat)
        # do preprocessing
        shift, bpos = preprocess_strong_suffix(pat, m)
        shift = preprocess_case2(shift, bpos, pat, m)

        #heuristic optimization
        
        charBeforeBorderLookup = defaultdict(list)
        
        #for i in range(len(pattern)):
        #   charBeforeBorderLookup = [{pattern[a.start()-1], a.start()-1}  for a in list(re.finditer(pattern[i:], pattern))]

        i = 1
        while i < len(bpos)-1:
            if pat[bpos[i]:]:
                charBeforeBorderLookup[pat[bpos[i]:]].append({'previousChar':pat[i-1], 'index':i-1})
            i+=1
        #ako hocemo da imamo prethodno slovo od poslednjeg u nizu ponavljajucih suffixa, mada je nepotrebno, posto ako na njima padne svakako postoji sledeci
        #a ako njih prodje, nema potrebe da ih vise razmatramo u okviru iste iteracije
        #bposset = set(bpos)
        #for item in bposset:
        #    u bpossetu se nalaze jedinstveni indeksi; ovo se radi da bi se dodao karakter koji se nalazi prije poslednjeg u nizu ponavljajucih suffixa, a ovo uspijeva
        #    zahvaljujuci cinjenici da svaki od suffixa u bpos nizu sadrzi indeks posljednjeg suffixa; iteracijom kroz ovaj bposset dodajemo prethodnike posljednjih suffixa
        #    if item <= m and pat[item:]:
        #        charBeforeBorderLookup[pat[item:]].append({'previousChar':pat[item-1], 'index':item-1})
        #        charBeforeBorderLookup[pat[item:]].reverse()
        #print(bposset)
        
        for key in charBeforeBorderLookup:
            charBeforeBorderLookup[key].reverse()
        skipped_alignments = 0
        foundList = []
        while s <= len(text) - len(pat):
            shift_amount = self.getNextShift(s, text, pat, shift, bpos, charBeforeBorderLookup, foundList)
            skipped_alignments += 1
            s += shift_amount
        skipped_alignments = len(text) - skipped_alignments
        listOfSkippedAlignments["Heuristic2"].append(skipped_alignments)
        print("     Total skipped alignments by heuristic 2 = " + str(skipped_alignments))
        return foundList


class Heuristic1and2:
    '''Search for a pattern in given text using 
    Boyer Moore algorithm with Good suffix rule '''

    def search(self, text, pat, n):
        # s is shift of the pattern with respect to text
        s = 0
        m = len(pat)
        # do preprocessing for heuristic 1
        bad_char_table = preprocess_bad_character(pat)
        # end preprocessing for heuristic 1
        # do preprocessing for heuristic 2
        shift, bpos = preprocess_strong_suffix(pat, m)
        shift = preprocess_case2(shift, bpos, pat, m)
        charBeforeBorderLookup = defaultdict(list)
        i = 1
        while i < len(bpos)-1:
            if pat[bpos[i]:]:
                charBeforeBorderLookup[pat[bpos[i]:]].append({'previousChar':pat[i-1], 'index':i-1})
            i+=1 
        for key in charBeforeBorderLookup:
            charBeforeBorderLookup[key].reverse()
        #finished preprocessing for heuristic 2
        heuristic1 = Heuristic1()
        heuristic2 = Heuristic2()
        skipped_alignments = 0
        foundListSuffix = []
        foundListCharacter = []
        while s <= len(text) - len(pat):
            shift_amount = max(heuristic2.getNextShift(s, text, pat, shift, bpos, charBeforeBorderLookup, foundListSuffix), heuristic1.getNextShift(s, text, pat, n, bad_char_table, foundListCharacter))
            skipped_alignments += 1
            s += shift_amount
        skipped_alignments = len(text) - skipped_alignments
        listOfSkippedAlignments["Heuristic 1 and 2"].append(skipped_alignments)
        print("     Total skipped alignments by heuristic 1 and 2 with n = " + str(n) +" is: " + str(skipped_alignments))
        return foundListSuffix

class BadCharacterAndGoodSuffixRuleHeuristic:
    def getNextShiftGoodSuffix(self, s, text, pat, shift, bpos, foundList):
  
        j = len(pat) - 1
          
        ''' Keep reducing index j of pattern while characters of 
            pattern and text are matching at this shift s'''
        while j >= 0 and pat[j] == text[s + j]:
            j -= 1
              
        ''' If the pattern is present at the current shift, 
            then index j will become -1 after the above loop '''
        if j < 0:
            foundList.append(s)
            return shift[0]
        else:
            '''pat[i] != pat[s+j] so shift the pattern 
            shift[j+1] times '''
            return shift[j+1]

    def search(self, text, pat, n):
        # s is shift of the pattern with respect to text
        s = 0
        m = len(pat)
        n = len(text)

        # do preprocessing
        shift, bpos = preprocess_strong_suffix(pat, m)
        shift = preprocess_case2(shift, bpos, pat, m)
        bad_char = Heuristic1()
        bad_char_table = preprocess_bad_character(pat)
        skipped_alignments = 0
        foundListSuffix = []
        foundListCharacter = []
        while s <= len(text) - len(pat):
            shift_amount = max(self.getNextShiftGoodSuffix(s, text, pat, shift, bpos, foundListSuffix), bad_char.getNextShift(s, text, pat, 1, bad_char_table, foundListCharacter))
            skipped_alignments += 1
            s += shift_amount
        skipped_alignments = len(text) - skipped_alignments
        listOfSkippedAlignments["Boyer Moore"].append(skipped_alignments)
        print("     Total skipped alignments by Boyer Moore full algorithm = " + str(skipped_alignments))
        return foundListSuffix

class NoHeuristic:
    def search(self, text, pattern, n):
        return [a.start() for a in list(re.finditer(pattern, text))]
    
class BoyesMooreAlgorithm:
    Heuristics = {"Heuristic1":Heuristic1(),
                "Heuristic2":Heuristic2(),
                "Heuristic1and2":Heuristic1and2(),

                "BadCharacterAndGoodSuffixRuleHeuristic":BadCharacterAndGoodSuffixRuleHeuristic(),
                 "NoHeuristic":NoHeuristic()}
    
    def search(self, heuristic, text, pattern, n):
        if (heuristic in self.Heuristics):
            start = time.time()
            print("     Starting search using: '" + heuristic + "'")
            result = list(self.Heuristics[heuristic].search(text,pattern, n))
            print("     Matches at: " + str(result))
            print("     Searching finished! Time spent: ", time.time() - start)
        else:
            print("     Heuristic '" + heuristic + "' is not implemented!")


def searchPattern(text, pattern, skip_param_heuristic1):
    print("Searching text: \"" + text + "\" for pattern: \"" + pattern + "\"")
    bm = BoyesMooreAlgorithm()
    bm.search("Heuristic1", text, pattern, skip_param_heuristic1)
    bm.search("Heuristic2", text, pattern, 1)
    bm.search("Heuristic1and2", text, pattern, skip_param_heuristic1)
    bm.search("BadCharacterAndGoodSuffixRuleHeuristic",  text, pattern, 1)
    #TO DO: Add assertion of lists with pattern found indices
    pattern_names_label.append(pattern)

searchPattern( "AAAAAAAAAAAAAAAA", "A", 2)
searchPattern( "SFGATFGACGAAACGAGTAGCSFGATAGACGA", "SFGATAGACGA", 2)
searchPattern( "CTATCGAAGTAGCCGATTAGC", "CGA", 2)

x = np.arange(len(pattern_names_label))  # the label locations
width = 0.8  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(x - 2*width/4, listOfSkippedAlignments["Heuristic1"], width, label='Heuristic1')
rects2 = ax.bar(x - width/4, listOfSkippedAlignments["Heuristic2"], width, label='Heuristic2')
rects3 = ax.bar(x + width/4, listOfSkippedAlignments["Heuristic 1 and 2"], width, label='Heuristic 1 and 2')
rects4 = ax.bar(x + 2*width/4, listOfSkippedAlignments["Boyer Moore"], width, label='Boyer Moore')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Skipped alignments')
ax.set_title('Skipped alignments by heuristic and pattern')
ax.set_xticks(x)
ax.set_xticklabels(pattern_names_label)
ax.legend()

ax.bar_label(rects1, padding=3)
ax.bar_label(rects2, padding=3)
ax.bar_label(rects3, padding=3)
ax.bar_label(rects4, padding=3)

fig.tight_layout()
plt.plot(heuristicNames, listOfSkippedAlignments)
plt.show()
#TO DO: Ako se bude imalo vremena, izdvojiti preprocesiranje za heuristiku 2 u zasebnu funkciju