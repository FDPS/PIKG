#ifndef H_PIKG_MNCORE2_UTIL
#define H_PIKG_MNCORE2_UTIL

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include<mncl/host/cl/cl.h>
//#include<mncl/host/cl/mncl_ext.h>
//#include<mncl/host/cl/util.h>
//#include<codegen/base/log.h>

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

  //std::cout << vsm << std::endl;

  cl_program program = clCreateProgramWithIL(context, reinterpret_cast<const unsigned char*>(vsm.c_str()), std::string(vsm).size(), &status);
  cl_kernel ret = clCreateKernel(program, kernelname.c_str(), &status); assert(status == CL_SUCCESS);
  assert(status == CL_SUCCESS);

  return ret;
}

#if 0
void print_dram(cl_context& context, cl_command_queue& queue, const size_t size, const size_t offset = 0){
  double* dram = new double[size];
  cl_int status;
  cl_mem dramBuffer = clCreateBufferWithAttributes(context, CL_MEM_READ_ONLY, 0, sizeof(PIKG::F64)*offset, sizeof(PIKG::F64)*size, nullptr, &status); PIKG_MN_CHECK(status);

  status = clEnqueueReadBuffer(queue, dramBuffer, CL_TRUE, 0, sizeof(PIKG::F64)*size, dram, 0, nullptr, nullptr); PIKG_MN_CHECK(status);
  const size_t chunk = 8;
  for(size_t i=0;i<size/chunk;i++){
    std::cout << std::setw(9) << chunk*i;
    for(size_t j=0;j<chunk;j++) std::cout << " " << std::setw(12) << std::setprecision(8) << dram[chunk*i+j];
    std::cout << std::endl;
  }

  status = clReleaseMemObject(dramBuffer); PIKG_MN_CHECK(status);
  delete[] dram;
}

void write_r_inst(std::ostream& os, const cl_mem& mem, const int length, const uint64_t* payload){
  cl_uint chip_id;
  clGetMemObjectInfo(mem,CL_MEM_DRAM_ID,  sizeof(cl_uint),&chip_id,NULL);
  cl_uint begin;
  clGetMemObjectInfo(mem,CL_MEM_DRAM_ADDR,sizeof(cl_uint),&begin,NULL);
  begin = begin / 8; // clGetMemObjectInfo returns byte address

  for(int iter=0;iter<(length+2047)/2048;iter++){
    const int length_tmp = std::min(2048, length-iter*2048);
    os << "r " << std::hex << std::setw(1)                      << chip_id
       << " "  << std::hex << std::setw(9) << std::setfill('0') << begin + 2048*iter
       << " "  << std::hex << std::setw(3) << std::setfill('0') << length_tmp
       << " ";
    for(int i=0;i<length_tmp;i++)
       os << std::hex << std::setw(16) << std::setfill('0') << payload[i + 2048*iter];
    os << std::endl;
  }
}
#endif
#endif
