#ifndef CSVROW
#define CSVROW

// Necessary cpp libraries
#include <stdio.h>
#include <iostream>
#include <string> // the C++ String Class, easier to deal with than char arrays for this application
#include <fstream>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <cmath> // math utility functions
#include <ctime> // timers for performance testing
#include <iterator>
#include <cstring>

class CSVRow
{
public:
   std::string const& operator[](std::size_t index) const
   {
      return m_data[index];
   }
   std::size_t size() const
   {
      return m_data.size();
   }
   void readNextRow(std::istream& str)
   {
      std::string         line;
      std::getline(str, line);

      std::stringstream   lineStream(line);
      std::string         cell;

      m_data.clear();
      while (std::getline(lineStream, cell, ','))
      {
         m_data.push_back(cell);
      }
      // This checks for a trailing comma with no data after it.
      if (!lineStream && cell.empty())
      {
         // If there was a trailing comma then add an empty element.
         m_data.push_back("");
      }
   }
private:
   std::vector<std::string>    m_data;
};

#endif