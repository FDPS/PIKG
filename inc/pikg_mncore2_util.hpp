#ifndef H_PIKG_MNCORE2_UTIL
#define H_PIKG_MNCORE2_UTIL

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include<mncl/host/cl/cl.h>

cl_mem clCreateBufferWithAttributes(cl_context context, cl_mem_flags flags, int dram_id, int64_t addr, size_t size, void* host_ptr, cl_int* errcode_ret);

#define PIKG_MN_CHECK(STATUS) if(STATUS != CL_SUCCESS){ std::cout << "error: at line" << __LINE__ << ", PIKG_MN_CHECK_FAILURE"  << std::endl; }

cl_kernel load_kernel_from_file(cl_context& context, std::string filename, std::string kernelname, cl_int& status){

  std::string vsm;
  std::ifstream ifs(filename);
  if(!ifs.is_open()){
    std::cout << "error: failed to open vsm file!" << std::endl;
    exit(-1);
  }
  for(std::string line; std::getline(ifs,line);){
    vsm += line + "\n";
  }

  cl_program program = clCreateProgramWithIL(context, reinterpret_cast<const unsigned char*>(vsm.c_str()), std::string(vsm).size(), &status);
  cl_kernel ret = clCreateKernel(program, kernelname.c_str(), &status); assert(status == CL_SUCCESS);
  assert(status == CL_SUCCESS);

  return ret;
}

#endif
