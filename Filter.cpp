#include "FilterLib.h"
#include <cmath>
#include <iostream>

using namespace std;

int getMedian(double *arr, int length)
{
    int i, j;
    double temp;
    for (i = 0; i < length; i++)
        for (j = 0; j < length - i - 1; j++)
            if (arr[j] > arr[j + 1])
            {
                temp = arr[j + 1];
                arr[j + 1] = arr[j];
                arr[j] = temp;
            }

    if (length % 2 == 0)
        temp = (arr[length / 2] + arr[length / 2 + 1]) / 2;
    else
        temp = arr[(length - 1) / 2];
    return temp;
}

ndarray *GaussianFilter_3d(ndarray *arr, int kernel_size, double sigma)
{
    ndarray *arr_out;
    InitArray(arr_out, 3);
    copyArrayProperty(arr, arr_out);
    arr_out->data = new double[arr->length];
    int edge = int(kernel_size / 2);
    double value;
    /* 计算高斯滤波核 */
    int x, y, _x, _y;
    double kernel_sum = 0;
    double *kernel;
    kernel = new double[kernel_size * kernel_size];
    for (x = 0; x < kernel_size; x++)
    {
        for (y = 0; y < kernel_size; y++)
        {
            _x = x - edge;
            _y = y - edge;

            kernel[kernel_size * x + y] = exp(-(_x * _x + _y * _y) / (2 * sigma * sigma)) / (2 * 3.14159265 * sigma * sigma);
            kernel_sum += kernel[kernel_size * x + y];
        }
    }

    for (x = 0; x < kernel_size; x++)
    {
        for (y = 0; y < kernel_size; y++)
        {
            kernel[kernel_size * x + y] /= kernel_sum;
        }
    }

    /* 高斯滤波 */
    for (int k = 0; k < arr->size[0]; k++)
        for (int i = 0; i < arr->size[1]; i++)
            for (int j = 0; j < arr->size[2]; j++)
            {
                value = 0;
                for (int di = -edge; di < edge + 1; di++)
                    for (int dj = -edge; dj < edge + 1; dj++)
                        if (((i + di) >= 0) && ((j + dj) >= 0) && ((i + di) < arr->size[1]) && ((j + dj) < arr->size[2]))
                            value += arr->data[k * arr->size[1] * arr->size[2] + (i + di) * arr->size[2] + j + dj] * kernel[(di + edge) * kernel_size + dj + edge];
                arr_out->data[k * arr->size[1] * arr->size[2] + i * arr->size[2] + j] = value;
            }
    return arr_out;
}

ndarray *MedianFilter_3d(ndarray *arr, int kernel_size)
{
    ndarray *arr_out;
    InitArray(arr_out, 3);
    copyArrayProperty(arr, arr_out);
    arr_out->data = new double[arr->length];
    int edge = int(kernel_size / 2);
    double *kernel;
    kernel = new double[kernel_size * kernel_size];
    int count;
    for (int k = 0; k < arr->size[0]; k++)
        for (int i = 0; i < arr->size[1]; i++)
            for (int j = 0; j < arr->size[2]; j++)
            {
                count = 0;
                for (int di = -edge; di < edge + 1; di++)
                    for (int dj = -edge; dj < edge + 1; dj++)
                        if (((i + di) >= 0) && ((j + dj) >= 0) && ((i + di) < arr->size[1]) && ((j + dj) < arr->size[2]))
                        {
                            kernel[count] = arr->data[k * arr->size[1] * arr->size[2] + (i + di) * arr->size[2] + j + dj];
                            count++;
                        }
                arr_out->data[k * arr->size[1] * arr->size[2] + i * arr->size[2] + j] = getMedian(kernel, count);
            }
    return arr_out;
}

ndarray *AverageFilter_3d(ndarray *arr, int kernel_size)
{
    ndarray *arr_out;
    InitArray(arr_out, 3);
    copyArrayProperty(arr, arr_out);
    arr_out->data = new double[arr->length];
    int edge = int(kernel_size / 2);
    int count;
    double value;
    for (int k = 0; k < arr->size[0]; k++)
        for (int i = 0; i < arr->size[1]; i++)
            for (int j = 0; j < arr->size[2]; j++)
            {
                count = 0;
                value = 0;
                for (int di = -edge; di < edge + 1; di++)
                    for (int dj = -edge; dj < edge + 1; dj++)
                        if (((i + di) >= 0) && ((j + dj) >= 0) && ((i + di) < arr->size[1]) && ((j + dj) < arr->size[2]))
                        {
                            value += arr->data[k * arr->size[1] * arr->size[2] + (i + di) * arr->size[2] + j + dj];
                            count++;
                        }
                arr_out->data[k * arr->size[1] * arr->size[2] + i * arr->size[2] + j] = value/count;
            }
    return arr_out;
}

extern "C" __declspec(dllexport) void GammaCorrection(double *arr_input, double *arr_output, int *size, double gamma)
{
    int length;
    length = size[0] * size[1] * size[2];

    for (int i = 0; i < length; i++)
        arr_output[i] = pow(arr_input[i] / 255., gamma);
}

extern "C" __declspec(dllexport) void GaussianFilter(double *arr_input, double *arr_output, int *size, int kernel_size, double sigma)
{
    ndarray *arr;
    InitArray(arr, 3);
    arr->size = size;
    arr->length = arr->size[0] * arr->size[1] * arr->size[2];
    arr->data = arr_input;

    ndarray *filtered = GaussianFilter_3d(arr, kernel_size, sigma);
    for (int i = 0; i < arr->length; i++)
        arr_output[i] = filtered->data[i];
    arr->data = NULL;
    arr->size = NULL;
    FreeArray(arr);
    FreeArray(filtered);
}

extern "C" __declspec(dllexport) void MedianFilter(double *arr_input, double *arr_output, int *size, int kernel_size)
{
    ndarray *arr;
    InitArray(arr, 3);
    arr->size = size;
    arr->length = arr->size[0] * arr->size[1] * arr->size[2];
    arr->data = arr_input;

    ndarray *filtered = MedianFilter_3d(arr, kernel_size);
    for (int i = 0; i < arr->length; i++)
        arr_output[i] = filtered->data[i];
    arr->data = NULL;
    arr->size = NULL;
    FreeArray(arr);
    FreeArray(filtered);
}

extern "C" __declspec(dllexport) void AverageFilter(double *arr_input, double *arr_output, int *size, int kernel_size)
{
    ndarray *arr;
    InitArray(arr, 3);
    arr->size = size;
    arr->length = arr->size[0] * arr->size[1] * arr->size[2];
    arr->data = arr_input;

    ndarray *filtered = AverageFilter_3d(arr, kernel_size);
    for (int i = 0; i < arr->length; i++)
        arr_output[i] = filtered->data[i];
    arr->data = NULL;
    arr->size = NULL;
    FreeArray(arr);
    FreeArray(filtered);
}