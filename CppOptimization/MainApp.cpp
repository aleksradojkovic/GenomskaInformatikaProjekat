#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <string>
#include <vector>
#include <iostream>

namespace py = pybind11;

std::vector<std::vector<size_t>> preprocessingForHeuristic1(const std::string& pattern)
{
    std::vector<std::vector<size_t>> table{ 256 };

    for (int i = pattern.length() - 1; i >= 0; --i)
    {
        table[pattern[i]].emplace_back(i);
    }

    return table;
}

class Heuristic1_cpp
{
public:
    unsigned GetNextShift(unsigned i, const std::string& text, const std::string& pattern, const std::vector<std::vector<size_t>>& table, std::vector<size_t>& foundList)
    {
        unsigned j = pattern.length() - 1;

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


PYBIND11_MODULE(CppOptimization, m) {
    py::class_<Heuristic1_cpp>(m, "Heur1")
        .def(py::init<>())
        .def("search", &Heuristic1_cpp::search);
}