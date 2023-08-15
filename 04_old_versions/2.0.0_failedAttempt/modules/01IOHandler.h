#ifndef IOHANDLER
#define IOHANDLER

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

#include "01CSVRow.h"

class IOHandler
{
public:
   IOHandler() {};
   double Cells(long row, long col);
   double fCells(long row, long col);
};

#endif