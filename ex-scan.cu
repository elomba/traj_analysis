// Filename: csort.cu
// nvcc -c -arch sm_13 csort.cu
//#include <thrust/device_vector.h>
//#include <thrust/sort.h>
#include <thrust/scan.h>
#include <thrust/device_ptr.h>
#include <thrust/execution_policy.h>
extern "C" {
//Sort for integer arrays
void scan_int_wrapper( int *data, int N, int start)
{
// Wrap raw pointer with a device_ptr
thrust::device_ptr <int> dev_ptr(data);
// Use device_ptr in Thrust sort algorithm
//  int data[6] = {1, 0, 2, 2, 1, 3};
thrust::exclusive_scan(thrust::device,  dev_ptr,  dev_ptr+N,  dev_ptr, start);
}

//Sort for float arrays
void scan_float_wrapper( float *data, int N, float start)
{
thrust::device_ptr <float> dev_ptr(data);
 thrust::exclusive_scan(thrust::device, dev_ptr, dev_ptr+N,  dev_ptr, start);
}

//Sort for double arrays
void scan_double_wrapper( double *data, int N, double start)
{
thrust::device_ptr <double> dev_ptr(data);
 thrust::exclusive_scan(thrust::device, dev_ptr, dev_ptr+N,  dev_ptr, start);
}

}
