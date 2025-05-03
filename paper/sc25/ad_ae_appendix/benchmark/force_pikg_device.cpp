#include<particle_simulator.hpp>
#include "force_device.hpp"

#if defined(USE_GPU)
#include "cuda_pointer.h"
#endif
#include "kernel_pikg.hpp"

constexpr int N_EPJ_MAX = 128*1024;
EPJGrav* espj_buffer[N_WALK_LIMIT];
int      n_espj[N_WALK_LIMIT];


PS::S32 DispatchKernel(
                       const PS::S32          tag,
                       const PS::S32          n_walk,
                       const EPIGrav          *epi[],
                       const PS::S32          n_epi[],
                       const EPJGrav          *epj[],
                       const PS::S32          n_epj[],
                       const PS::SPJMonopole *spj[],
                       const PS::S32          n_spj[]){
    static_assert(sizeof(PS::SPJMonopole) == sizeof(EPJGrav));
#if 0
    DispatchCalcGravityEpEp(tag,n_walk,epi,n_epi,epj,n_epj,spj,n_spj,false);
#else
    static bool first = true;
    if(first){
      for(int w=0;w<N_WALK_LIMIT;w++){
	espj_buffer[w] = new EPJGrav[N_EPJ_MAX];
      }
      first = false;
    }
    const auto time_offset = PS::GetWtime();
    PS::S32 n = 0;
    PS::S32 n_max = 0;
    assert(n_walk > 0);
    for(int w=0; w<n_walk; w++){
      n_espj[w] = n_epj[w] + n_spj[w];
      n += n_espj[w];
      n_max = std::max(n_max,n_espj[w]);
    }
    assert(n_max > 0);
    #pragma omp parallel for
    for(int ii=0;ii<n_max*n_walk;ii++){
      const int w = ii / n_max;
      const int i = ii % n_max;
      const int offset = n_epj[w];
      const EPJGrav* pspj = (EPJGrav*)spj[w];
      if(i < n_epj[w]){
	espj_buffer[w][i] = epj[w][i];
      }else if(i < n_espj[w]){
	espj_buffer[w][i+offset] = pspj[i-offset];
      }
    }
    auto elapsed = PS::GetWtime() - time_offset;
    //std::cout << "copy_espj = " << elapsed << " ( " << 2/*readw/rite*/*4/*byte*/*4/*elements*/*n / elapsed * 1e-9 << " GB/s)" << std::endl;
    DispatchCalcGravityEpEp(tag,n_walk,epi,n_epi,(const EPJGrav**)espj_buffer,n_espj);
#endif
    return 0;
}
#if 0
PS::S32 DispatchKernelIndex(const PS::S32 tag,
                       const PS::S32 n_walk,
                       const EPIGrav ** epi,
                       const PS::S32 *  n_epi,
                       const PS::S32 ** id_epj,
                       const PS::S32 *  n_epj,
		       const PS::S32 ** id_spj,
                       const PS::S32  * n_spj,
                       const EPJGrav * epj,
		       const PS::S32   n_epj_tot,
                       const PS::SPJMonopole * spj,
		       const PS::S32   n_spj_tot,
		       const bool      send_flag){
    static_assert(sizeof(PS::SPJMonopole) == sizeof(EPJGrav));
    return DispatchCalcGravityEpEpIndex(tag,n_walk,epi,n_epi,id_epj,n_epj,id_spj,n_spj,epj,n_epj_tot,spj,n_spj_tot,send_flag);
}
#endif

PS::S32 RetrieveKernel(const PS::S32 tag,
                       const PS::S32 n_walk,
                       const PS::S32 ni[],
                       ForceGrav    *force[])
{
    RetrieveCalcGravityEpEp(tag,n_walk,ni,force);
    return 0;
}
