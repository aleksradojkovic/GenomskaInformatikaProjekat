#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <array>
#include <algorithm>

namespace py = pybind11;

using MatchList = std::vector<unsigned>;


class Algorithm
{
protected:
    MatchList m_FoundList;

    virtual void PreProcess(const std::string& pattern) = 0;
    virtual int GetNextShift(int s, const std::string& text, const std::string& pattern) = 0;
public:
    void search(const std::string& text, const std::string& pattern)
    {
        int s = 0, shift_amount;
        
        PreProcess(pattern);

        while (s <= text.length() - pattern.length())
        {
            shift_amount = GetNextShift(s, text, pattern);
            s += shift_amount;
        }
    }

    const MatchList& GetMatches() const { return m_FoundList; }

};

void preprocess_bad_character(const std::string& pattern, std::array<int, 256>& lastOccurence)
{
    for (unsigned i = 0; i < pattern.length(); ++i)
    {
        lastOccurence[pattern.at(i)] = i;
    }
}

std::vector<std::vector<size_t>> preprocessingForHeuristic1(const std::string& pattern)
{
    std::vector<std::vector<size_t>> table{ 256 };

    for (int i = pattern.length() - 1; i >= 0; --i)
    {
        table[pattern[i]].emplace_back(i);
    }

    return table;
}


struct PrevChar
{
    PrevChar(char prevChar, int ind) : previousChar(prevChar), index(ind) {}
    char previousChar;
    int index;

};

std::map<std::string, std::vector<PrevChar>> preprocessingForHeuristic2(const std::vector<int>& bpos, const std::string& pattern)
{
    std::map<std::string, std::vector<PrevChar>> charBeforeBorderLookup;
    unsigned i = 1;
    while (i < bpos.size() - 1)
    {
        if (!pattern.substr(bpos[i]).empty())
        {
            auto it = charBeforeBorderLookup.find(pattern.substr(bpos[i]));
            if (it == charBeforeBorderLookup.end())
            {
                charBeforeBorderLookup.insert({ pattern.substr(bpos[i]), {PrevChar(pattern[i - 1], i - 1)} });
            }
            else
            {
                it->second.push_back(PrevChar(pattern[i - 1], i - 1));
            }
        }
        ++i;
    }

    for (auto& item : charBeforeBorderLookup)
    {
        std::reverse(item.second.begin(), item.second.end());
    }

    return charBeforeBorderLookup;
}

void preprocess_strong_suffix(const std::string& pattern, std::vector<int>& shift, std::vector<int>& bpos)
{
    bpos.resize(pattern.length() + 1);
    shift.resize(pattern.length() + 1);

    int i = pattern.length();
    int j = i + 1;

    bpos[i] = j;

    while (i > 0)
    {
        while (j <= pattern.length() && pattern[i - 1] != pattern[j - 1])
        {
            if (shift[j] == 0)
            {
                shift[j] = j - i;
            }
           
            j = bpos[j];
        }

        i -= 1;
        j -= 1;
        bpos[i] = j;
    }
}

void preprocess_case2(std::vector<int>& shift, const std::vector<int>& bpos, const std::string& pattern)
{
    int j = bpos[0];
    
    for (int i = 0; i <= pattern.length(); i++)
    {
        if (shift[i] == 0)
        {
            shift[i] = j;
        }

        if (i == j)
        {
            j = bpos[j];
        }
    }
}

class Heuristic1_cpp
{
public:
    unsigned GetNextShift(unsigned i, const std::string& text, const std::string& pattern, const std::vector<std::vector<size_t>>& table, std::vector<size_t>& foundList)
    {
        int j = pattern.length() - 1;

        while (j >= 0 and pattern[j] == text[i + j])
            --j;

        if (j < 0)
        {
            foundList.emplace_back(i);
            return 1;
        }
        else
        {
            if (table[text[i + j]].size())
            {
                for (unsigned int cnt = 0; cnt < table[text[i + j]].size(); ++cnt)
                {
                    const size_t& item = table[text[i + j]].at(cnt);
                    if (item < j)
                    {
                        return j - item;
                    }
                }
            }
            return j + 1;
        }
    }

    std::vector<size_t> search(const std::string& text, const std::string& pattern)
    {
        size_t i = 0, shift = 0;
        std::vector<std::vector<size_t>> table = preprocessingForHeuristic1(pattern);
        std::vector<size_t> foundList;

        while (i <= text.length() - pattern.length())
        {
            shift = GetNextShift(i, text, pattern, table, foundList);
            i += shift;
        }

        return foundList;
    }

};

class Heuristic2_cpp
{
public:
    unsigned GetNextShift(int s, const std::string& text, const std::string& pattern, const std::vector<int>& shift, const std::vector<int>& bpos, const std::map<std::string, std::vector<PrevChar>>& charBeforeBorderLookup, std::vector<size_t>& foundList)
    {
        int j = pattern.length() - 1;
        while (j >= 0 and pattern[j] == text[s + j])
            --j;
        
        if (j < 0)
        {
            foundList.emplace_back(s);
            return shift[0];
        }
        else
        {
            bool foundCharMatch = false, foundHashLookup = false;

            if (pattern.substr(bpos[j]).size())
            {
                auto it = charBeforeBorderLookup.find(pattern.substr(bpos[j]));
                if (it != charBeforeBorderLookup.end())
                {
                    const std::vector<PrevChar>& listOfPrevChar = it->second;
                    
                    for (const PrevChar& item : listOfPrevChar)
                    {
                        if (item.previousChar == text[s + j - 1] && item.index < j)
                        {
                            foundCharMatch = true;
                            return j - item.index;
                        }
                    }
                }
            }

            if (!foundHashLookup)
            {
                return shift[j + 1];
            }
            else if (!foundCharMatch)
            {
                return shift[0];
            }
            else
            {
                //should never happen
                assert(false);
                return 1;
            }
        }
    }

    std::vector<size_t> search(const std::string& text, const std::string& pattern)
    {
        int s = 0, shift_amount;
        
        std::vector<int> shift, bpos;
        
        std::vector<size_t> foundList;
        
        preprocess_strong_suffix(pattern, shift, bpos);
        preprocess_case2(shift, bpos, pattern);

        std::map<std::string, std::vector<PrevChar>> charBeforeBorderLookup = preprocessingForHeuristic2(bpos, pattern);

        while (s <= text.length() - pattern.length())
        {
            shift_amount = GetNextShift(s, text, pattern, shift, bpos, charBeforeBorderLookup, foundList);
            s += shift_amount;
        }

        return foundList;
    }
};

class Heuristic1and2_cpp
{
public:
    std::vector<size_t> search(const std::string& text, const std::string& pattern)
    {
        int s = 0, shift_amount = 0;
        std::vector<std::vector<size_t>> table = preprocessingForHeuristic1(pattern);
        std::vector<size_t> foundListSuffix, foundListCharacter;

        std::vector<int> shift, bpos;

        preprocess_strong_suffix(pattern, shift, bpos);
        preprocess_case2(shift, bpos, pattern);

        std::map<std::string, std::vector<PrevChar>> charBeforeBorderLookup = preprocessingForHeuristic2(bpos, pattern);

        Heuristic1_cpp heuristic1 = Heuristic1_cpp();
        Heuristic2_cpp heuristic2 = Heuristic2_cpp();

        while (s <= text.length() - pattern.length())
        {
            shift_amount = std::max(heuristic2.GetNextShift(s, text, pattern, shift, bpos, charBeforeBorderLookup, foundListSuffix), heuristic1.GetNextShift(s, text, pattern, table, foundListCharacter));
            s += shift_amount;
        }
        return foundListSuffix;
    }
};

class BadCharacterAndGoodSuffixRuleHeuristic
{
public:
    unsigned GetNextShiftGoodSuffix(unsigned int s, const std::string& text, const std::string& pattern, const std::vector<int>& shift, const std::vector<int>& bpos, std::vector<int>& foundList)
    {
        int j = pattern.length() - 1;
        while (j >= 0 && pattern[j] == text[s + j])
        {
            --j;
        }

        if (j < 0)
        {
            foundList.push_back(s);
            return shift[0];
        }
        else
        {
            return shift[j + 1];
        }
    }


    unsigned GetNextShiftBadChar(unsigned int s, const std::string& text, const std::string& pattern, const std::array<int, 256>& lastOccurence, std::vector<int>& foundList)
    {
        int j = pattern.length() - 1;
        while (j >= 0 and pattern[j] == text[s + j])
        {
            --j;
        }

        if (j < 0)
        {
            foundList.push_back(s);
        }
        else
        {
            if (lastOccurence[text[s + j]] != -1)
            {
                return std::max(1, j - lastOccurence[text[s + j]]);
            }
        }   
        return 1;
    }

    std::vector<int> search(const std::string& text, const std::string& pattern)
    {
        unsigned s = 0, shift_amount = 0;
        std::vector<int> foundListSuffix, foundListCharacter;

        std::vector<int> shift, bpos;

        std::array<int, 256> lastOccurence;
        std::fill(lastOccurence.begin(), lastOccurence.end(), -1);

        preprocess_strong_suffix(pattern, shift, bpos);
        preprocess_case2(shift, bpos, pattern);
        preprocess_bad_character(pattern, lastOccurence);

        while (s <= text.length() - pattern.length())
        {
            shift_amount = std::max(GetNextShiftGoodSuffix(s, text, pattern, shift, bpos, foundListSuffix), GetNextShiftBadChar(s, text, pattern, lastOccurence, foundListCharacter));
            s += shift_amount;
        }

        return foundListCharacter;
    }
};


PYBIND11_MODULE(CppOptimization, m) {
    py::class_<Heuristic1_cpp>(m, "Heur1")
        .def(py::init<>())
        .def("search", &Heuristic1_cpp::search);
    py::class_<Heuristic2_cpp>(m, "Heur2")
        .def(py::init<>())
        .def("search", &Heuristic2_cpp::search);
    py::class_<Heuristic1and2_cpp>(m, "Heur12")
        .def(py::init<>())
        .def("search", &Heuristic1and2_cpp::search);
    py::class_<BadCharacterAndGoodSuffixRuleHeuristic>(m, "HeurNative")
        .def(py::init<>())
        .def("search", &BadCharacterAndGoodSuffixRuleHeuristic::search);
}