#include "cell.hpp"


bool operator==(int cellId, Cell cell)
{
    if(cellId == cell.id)
        return true;
    else
        return false;
}

// bool operator==(Cell &cell, int &cellId)
// {
//     if(cellId == cell.id)
//         return true;
//     else
//         return false;
// }