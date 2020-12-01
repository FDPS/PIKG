#ifndef PIKG_CUDA_POINTER
#define PIKG_CUDA_POINTER

#include <cassert>
#include <pikg_vector.hpp>

namespace PIKG{
  template <typename T>
  class CUDAPointer{
  public:
    //T* hst = nullptr;
    T* dev = nullptr;
    U64 _size = 0;
    U64 _limit = 0;

    ~CUDAPointer(){
      free();
    }
    /*    
    void free_host(){
      if(hst != nullptr){
	free(hst);
	hst = nullptr;
      } 
    }
    */
    void free_dev(){
      if(dev != nullptr){
	cudaFree(dev);
	dev = nullptr;
      } 
    }
    void free(){
      //free_host();
      free_dev(); 
    }

    void malloc(const U64 n){
      free();
      //hst = (T*)malloc(n*sizeof(T));
      cudaMalloc((void**)&dev,n*sizeof(T));
      _size = n;
      _limit = n;

      //assert(hst != nullptr);
      assert(dev != nullptr);
    }
    void resize(const U64 n){
      reserve(n);
      _size = n;
    }

    void reserve(const U64 n){
      if(size==0 || _limit < n) malloc(n);
    }

    void h2d(const T* hst,const U64 n = _size){
      cudaMemcpy(dev,hst,n*sizeof(T),cudaMemcpyHostToDevice);
    }
    void d2h(T* hst) const {
      cudaMemcpy(hst,dev,size*sizeof(T),cudaMemcpyDeviceToHost);
    }
  };
};

#endif // PIKG_CUDA_POINTER
