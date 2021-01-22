#pragma once


struct ndarray
{
    double *data;
    int *size;
    int length;
    unsigned char ndim;
};

template <typename T1, typename T2>
void copyArrayProperty(T1 A, T2 B)
{
    B->length = A->length;
    for (int i = 0; i < A->ndim; i++)
    {
        B->size[i] = A->size[i];
    }
}


void InitArray(ndarray *&arr, int ndim);
void FreeArray(ndarray *&arr);