#include "01IOHandler.h"
#include "00MainProgram.h"

// returns a reference to the "cell" in question for reading or writing
// for IO code equivalency with original VBA version
double IOHandler::Cells(long row, long col)
{
   MainProgram garisom;
   double outDbl = 0.0;
   if (row >= 0 && row < DATAFILE_MAXROWS)
   {
      if (col >= 0 && col < DATAFILE_MAXCOLS)
      {
         return garisom.dataCells[row][col];
      }
   }

   return garisom.dummyDouble;
}

double IOHandler::fCells(long row, long col)
{
   MainProgram garisom;
   double outDbl = 0.0;
   if (row >= 0 && row < MAX_SUMMARY_ROWS)
   {
      if (col >= 0 && col < MAX_SUMMARY_COLS)
      {
         return garisom.finalOutCells[row][col];
      }
   }
   return garisom.dummyDouble;
}
