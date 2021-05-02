import CppOptimization

objHr1 = CppOptimization.Heur1()
objHr2 = CppOptimization.Heur2()
objHr12 = CppOptimization.Heur12()
objHrNtv = CppOptimization.HeurNative()

objHr1.search("ABCACAABABABABBABABBBABBABABBCACABABABABABAB", "CACA")
print(objHr1.getMatches())

objHr2.search("ABCACAABABABABBABABBBABBABABBCACABABABABABAB", "CACA")
print(objHr2.getMatches())

#print(objHr12.search("ABCACAABABABABBABABBBABBABABBCACABABABABABAB", "CACA"))

objHrNtv.search("ABCACAABABABABBABABBBABBABABBCACABABABABABAB", "CACA")
print(objHrNtv.getMatches())
