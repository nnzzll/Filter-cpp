#include "FilterLib.h"


void InitArray(ndarray *&arr, int ndim)
{
    int i;
    arr = new ndarray;
    arr->data = nullptr;
    arr->ndim = ndim;
    arr->size = new int[ndim];
    arr->length = 0;
    for (i = 0; i < ndim; i++)
        arr->size[i] = 0;
}


void FreeArray(ndarray *&arr)
{
    if (arr)
    {
        if (arr->data)
            delete[] arr->data;
        if (arr->size)
            delete[] arr->size;
        delete arr;
    }
}