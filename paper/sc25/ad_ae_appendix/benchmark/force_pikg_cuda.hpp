#include "user-defined.hpp"
#include "pikg_cuda_pointer.hpp"
#include "pikg_vector.hpp"
struct EpiGPU{
float3 pos;
};
struct EpjGPU{
float mass;
float3 pos;
};
struct ForceGPU{
float3 acc;
float pot;
};
enum{
  N_THREAD_GPU = 32,
  N_WALK_LIMIT = 1000,
  NI_LIMIT     = N_WALK_LIMIT*1000,
  NJ_LIMIT     = N_WALK_LIMIT*10000,
};
inline __device__ ForceGPU inner_kernel(
				     EpiGPU epi,
				     EpjGPU epj,
				     ForceGPU force)
{
float __fkg_tmp0;
float __fkg_tmp1;
float mr3_inv;
float mr_inv;
float r2;
float r2_inv;
float r_inv;
float3 rij;
rij.x = (epi.pos.x-epj.pos.x);
rij.y = (epi.pos.y-epj.pos.y);
rij.z = (epi.pos.z-epj.pos.z);
__fkg_tmp1 = (rij.x*rij.x+0.0009765625f);
__fkg_tmp0 = (rij.y*rij.y+__fkg_tmp1);
r2 = (rij.z*rij.z+__fkg_tmp0);
r_inv = rsqrt(r2);
r2_inv = (r_inv*r_inv);
mr_inv = (epj.mass*r_inv);
mr3_inv = (r2_inv*mr_inv);
force.acc.x = (force.acc.x - mr3_inv*rij.x);
force.acc.y = (force.acc.y - mr3_inv*rij.y);
force.acc.z = (force.acc.z - mr3_inv*rij.z);
force.pot = (force.pot-mr_inv);
  return force;
}
__device__ ForceGPU ForceKernel_1walk(
				    EpjGPU *jpsh,
				    const EpiGPU ipos,
				    const int     id_walk,
				    const int    *ij_disp,
				    const EpjGPU *epj, 
				    ForceGPU accp){
  const int tid = threadIdx.x;
  const int j_head = ij_disp[id_walk  ];
  const int j_tail = ij_disp[id_walk+1];
  for(int j=j_head; j<j_tail; j+=N_THREAD_GPU){
    jpsh[tid] = epj[j+tid];
    if(j_tail-j < N_THREAD_GPU){
      for(int jj=0; jj<j_tail-j; jj++){
        accp = inner_kernel( ipos, jpsh[jj], accp);
      }
    }else{
#pragma unroll 32
      for(int jj=0; jj<N_THREAD_GPU; jj++){
        accp = inner_kernel( ipos, jpsh[jj], accp);
      }
    }
  }
  return accp;
}
__device__ ForceGPU ForceKernel_2walk(
				    EpjGPU jpsh[2][N_THREAD_GPU],
				    const EpiGPU  ipos,
				    const int     id_walk,
				    const int     iwalk0,
				    const int     iwalk1,
				    const int    *ij_disp,
				    const EpjGPU *epj, 
				    ForceGPU accp){
  const int jbeg0 = ij_disp[iwalk0];
  const int jbeg1 = ij_disp[iwalk1];
  const int jend0 = ij_disp[iwalk0 + 1];
  const int jend1 = ij_disp[iwalk1 + 1];
  const int nj0   = jend0 - jbeg0;
  const int nj1   = jend1 - jbeg1;
  const int nj_longer  = nj0 > nj1 ? nj0 : nj1;
  const int nj_shorter = nj0 > nj1 ? nj1 : nj0;
  const int walk_longer= nj0 > nj1 ? 0 : 1;
  const int jbeg_longer = nj0 > nj1 ? jbeg0 : jbeg1;
  const int mywalk = id_walk==iwalk0 ? 0 : 1;
  const int tid = threadIdx.x;
  for(int j=0; j<nj_shorter; j+=N_THREAD_GPU){
    jpsh[0][tid] = epj[jbeg0 + j + tid];
    jpsh[1][tid] = epj[jbeg1 + j + tid];
    if(nj_shorter-j < N_THREAD_GPU){
      for(int jj=0; jj<nj_shorter-j; jj++){
        accp = inner_kernel( ipos, jpsh[mywalk][jj], accp);
      }
    }else{
#pragma unroll 32
      for(int jj=0; jj<N_THREAD_GPU; jj++){
        accp = inner_kernel( ipos, jpsh[mywalk][jj], accp);
      }
    }
  }
  for(int j=nj_shorter; j<nj_longer; j+=N_THREAD_GPU){
    jpsh[0][tid] = epj[jbeg_longer + j + tid];
    int jrem = nj_longer - j;
    if(jrem < N_THREAD_GPU){
      for(int jj=0; jj<jrem; jj++){
	     if(mywalk == walk_longer)
          accp = inner_kernel( ipos, jpsh[0][jj], accp);
      }
    }else{
#pragma unroll 32
      for(int jj=0; jj<N_THREAD_GPU; jj++){
	     if(mywalk == walk_longer)
          accp = inner_kernel( ipos, jpsh[0][jj], accp);
      }
    }
  }
  return accp;
}
__device__ ForceGPU ForceKernel_multiwalk(
					const EpiGPU ipos,
					const int     id_walk,
					const int    *ij_disp,
					const EpjGPU *epj, 
					ForceGPU accp){
  const int j_head = ij_disp[id_walk  ];
  const int j_tail = ij_disp[id_walk+1];

  for(int j=j_head; j<j_tail; j++){
    EpjGPU jp = epj[j];
    accp = inner_kernel( ipos, jp, accp);
  }
  return accp;
}
__global__ void CalcGravityEpEp_cuda(
                  const int    * ij_disp,
                  const int    * walk,
                  const EpiGPU * epi,
                  const EpjGPU * epj, 
                  ForceGPU     * force){
  int tid = blockDim.x * blockIdx.x + threadIdx.x;
  const EpiGPU ip = epi[tid];
  const int id_walk = walk[tid];
  ForceGPU accp;
    accp.acc.x = 0.0f;
    accp.acc.y = 0.0f;
    accp.acc.z = 0.0f;
    accp.pot = 0.0f;

  int t_head = blockDim.x * blockIdx.x;
  int t_tail = t_head + N_THREAD_GPU - 1;
  int nwalk_in_block = 1 + (walk[t_tail] - walk[t_head]);

  __shared__ EpjGPU jpsh[2][N_THREAD_GPU];

  if(1 == nwalk_in_block){
    accp = ForceKernel_1walk(jpsh[0], ip, id_walk, ij_disp, epj, accp);
  } else if(2 == nwalk_in_block){
    int iwalk0 = walk[t_head];
    int iwalk1 = walk[t_tail];
    accp = ForceKernel_2walk(jpsh, ip, id_walk, iwalk0, iwalk1, ij_disp, epj, accp);
  } else{
    accp = ForceKernel_multiwalk(ip, id_walk, ij_disp, epj, accp);
  }
  force[tid] = accp;
}
static PIKG::CUDAPointer<EpiGPU>   dev_epi;
static PIKG::CUDAPointer<EpjGPU>   dev_epj;
static PIKG::CUDAPointer<ForceGPU> dev_force;
static PIKG::CUDAPointer<int>  ij_disp;
static PIKG::CUDAPointer<int>   walk;
void initialize_CalcGravityEpEp(){
}
PIKG::S32 DispatchCalcGravityEpEp(const PIKG::S32          tag,
                   const PIKG::S32          n_walk,
                   const EPIGrav    *epi[],
                   const PIKG::S32          n_epi[],
                   const EPJGrav    *epj[],
                   const PIKG::S32          n_epj[])
{
  assert(n_walk <= N_WALK_LIMIT);
  static bool init_call = true;
  if(init_call){
    dev_epi  .allocate(NI_LIMIT);
    dev_epj  .allocate(NJ_LIMIT);
    dev_force.allocate(NI_LIMIT);
    ij_disp  .allocate(N_WALK_LIMIT+2);
    walk     .allocate(NI_LIMIT);
    init_call = false;
  }
  ij_disp[0] = 0;
  for(int k=0; k<n_walk; k++){
    ij_disp[k+1]  = ij_disp[k] + n_epj[k];
  }
  ij_disp[n_walk+1] = ij_disp[n_walk];

  assert(ij_disp[n_walk] < NJ_LIMIT);
  ij_disp.h2d(n_walk + 2);


  int ni_tot = 0;
  int nj_tot = 0;
  for(int iw=0; iw<n_walk; iw++){
    for(int i=0; i<n_epi[iw]; i++){
      dev_epi[ni_tot].pos.x = epi[iw][i].pos.x;
      dev_epi[ni_tot].pos.y = epi[iw][i].pos.y;
      dev_epi[ni_tot].pos.z = epi[iw][i].pos.z;
      walk[ni_tot] = iw;
      ni_tot++;
    }
    for(int j=0; j<n_epj[iw]; j++){
      dev_epj[nj_tot].mass = epj[iw][j].mass;
      dev_epj[nj_tot].pos.x = epj[iw][j].pos.x;
      dev_epj[nj_tot].pos.y = epj[iw][j].pos.y;
      dev_epj[nj_tot].pos.z = epj[iw][j].pos.z;
      nj_tot++;
    }
  }
  assert(ni_tot < NI_LIMIT);
  int ni_tot_reg = ni_tot;
  if(ni_tot_reg % N_THREAD_GPU){
    ni_tot_reg /= N_THREAD_GPU;
    ni_tot_reg++;
    ni_tot_reg *= N_THREAD_GPU;
  }
  for(int i=ni_tot; i<ni_tot_reg; i++){
    walk[i] = n_walk;
  }

  walk.h2d(ni_tot_reg);
  dev_epi.h2d(ni_tot_reg);
  dev_epj.h2d(nj_tot);

  int nblocks  = ni_tot_reg / N_THREAD_GPU;
  int nthreads = N_THREAD_GPU;
  CalcGravityEpEp_cuda <<<nblocks, nthreads>>> (ij_disp, walk,  dev_epi, dev_epj, dev_force);

  return 0;
}
PIKG::S32 RetrieveCalcGravityEpEp(const PIKG::S32 tag,
                       const PIKG::S32 n_walk,
                       const PIKG::S32 ni[],
                       ForceGrav *force[])
{
  int ni_tot = 0;
  for(int k=0; k<n_walk; k++){
    ni_tot += ni[k];
  }
  dev_force.d2h(ni_tot);

  int n_cnt = 0;
  std::cout << n_walk << std::endl;
  for(int iw=0; iw<n_walk; iw++){
    for(int i=0; i<ni[iw]; i++){
      force[iw][i].acc.x = (force[iw][i].acc.x+dev_force[n_cnt].acc.x);
      force[iw][i].acc.y = (force[iw][i].acc.y+dev_force[n_cnt].acc.y);
      force[iw][i].acc.z = (force[iw][i].acc.z+dev_force[n_cnt].acc.z);
      force[iw][i].pot = (force[iw][i].pot+dev_force[n_cnt].pot);
      std::cout << force[iw][i].acc << std::endl;
      n_cnt++;
    }
  }
  return 0;
}
