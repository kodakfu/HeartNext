#ifndef BASE_H
#define BASE_H
#include<thrust/host_vector.h>
#include<thrust/device_vector.h>
typedef thrust::host_vector<double> host_vector_double;
typedef thrust::device_vector<double> device_vector_double;
typedef thrust::host_vector<int> host_vector_int;
typedef thrust::device_vector<int> device_vector_int;
#endif
