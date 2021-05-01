import CppOptimization

objHr1 = CppOptimization.Heur1()
objHr2 = CppOptimization.Heur2()
objHr12 = CppOptimization.Heur12()
objHrNtv = CppOptimization.HeurNative()
print(objHr1.search("ABCACAABABABABBABABBBABBABABBCACABABABABABAB", "CACA"))
print(objHr2.search("ABCACAABABABABBABABBBABBABABBCACABABABABABAB", "CACA"))
print(objHr12.search("ABCACAABABABABBABABBBABBABABBCACABABABABABAB", "CACA"))
print(objHrNtv.search("ABCACAABABABABBABABBBABBABABBCACABABABABABAB", "CACA"))
