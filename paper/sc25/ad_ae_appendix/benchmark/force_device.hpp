#pragma once
#include <iostream>
#include <particle_simulator.hpp>
#include "user-defined.hpp"
PS::S32 DispatchKernel(const PS::S32 tag,
                       const PS::S32 n_walk,
                       const EPIGrav ** epi,
                       const PS::S32 *  n_epi,
                       const EPJGrav ** epj,
                       const PS::S32 *  n_epj,
                       const PS::SPJMonopole ** spj,
                       const PS::S32  * n_spj);
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
		       const bool      send_flag);
#endif

PS::S32 RetrieveKernel(const PS::S32 tag,
                       const PS::S32 n_walk,
                       const PS::S32 * ni,
                       ForceGrav      ** force);
