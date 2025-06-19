#include "user_defined.h"
#include "hydro_force.hpp"
#include "hydro_force_reference.hpp"

#include<iostream>
#include<fstream>
#include<string>
#include<cassert>
void ReadInputs(EP_hydro**& epi, int*& n_epi, EP_hydro*& epj_send, int& nsend_epj, EP_hydro**& epj, int**& id_epj, int*& n_epj, int& n_walk){
  std::ifstream ifs("input.txt");
  std::string dummy;

  ifs >> dummy >> nsend_epj;
  //std::cout << dummy << " " << nsend_epj << std::endl;
  assert(dummy == "nsend_epj");

  // read epj
  epj_send = new EP_hydro[nsend_epj];
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  for(int j=0;j<nsend_epj;j++){
    ifs >> dummy; assert(dummy == "EPJ");
    int idx, rank;
    ifs >> idx
        >> rank
	>> epj_send[j].pos.x >> epj_send[j].pos.y >> epj_send[j].pos.z 
	>> epj_send[j].mass
	>> epj_send[j].vel.x >> epj_send[j].vel.y >> epj_send[j].vel.z
	>> epj_send[j].eng 
	>> epj_send[j].h
	>> epj_send[j].dens 
	>> epj_send[j].pres
	>> epj_send[j].gradh
	>> epj_send[j].snds
	>> epj_send[j].BalSW
	>> epj_send[j].alpha;
  }
  ifs >> dummy >> n_walk;
  assert(dummy == "n_walk");
  assert(n_walk > 0);
  epi = new EP_hydro*[n_walk];
  n_epi = new int[n_walk];
  epj = new EP_hydro*[n_walk];
  id_epj = new int*[n_walk];
  n_epj = new int[n_walk];
  for(int w=0;w<n_walk;w++){
    ifs >> dummy >> n_epi[w];
    assert(dummy == "n_epi");
    epi[w] = new EP_hydro[n_epi[w]];
    ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    for(int i=0;i<n_epi[w];i++){
      ifs >> dummy;
      assert(dummy == "EPI");
      int rank;
      ifs >> rank
	  >> epi[w][i].pos.x >> epi[w][i].pos.y >> epi[w][i].pos.z 
	  >> epi[w][i].mass
	  >> epi[w][i].vel.x >> epi[w][i].vel.y >> epi[w][i].vel.z 
	  >> epi[w][i].eng
	  >> epi[w][i].h
	  >> epi[w][i].dens
	  >> epi[w][i].pres
	  >> epi[w][i].gradh
	  >> epi[w][i].snds
	  >> epi[w][i].BalSW
	  >> epi[w][i].alpha;
    }
    ifs >> dummy >> n_epj[w];
    assert(dummy == "n_epj");
    epj[w] = new EP_hydro[n_epj[w]];
    id_epj[w] = new int[n_epj[w]];
    for(int j=0;j<n_epj[w];j++){
      int index, walk;
      ifs >> dummy >> walk >> index >> id_epj[w][j];
      assert(dummy == "epj_index");
      assert(w == walk);
      assert(index == j);
      ifs >> dummy; assert(dummy == "EPJ");
      ifs >> dummy //idx
          >> dummy // rank
          >> epj[w][j].pos.x >> epj[w][j].pos.y >> epj[w][j].pos.z 
	  >> epj[w][j].mass
	  >> epj[w][j].vel.x >> epj[w][j].vel.y >> epj[w][j].vel.z
	  >> epj[w][j].eng 
	  >> epj[w][j].h
	  >> epj[w][j].dens 
	  >> epj[w][j].pres
	  >> epj[w][j].gradh
	  >> epj[w][j].snds
	  >> epj[w][j].BalSW
	  >> epj[w][j].alpha;
    }
  }
  return;
}

int main(int argc,char** argv){
    EP_hydro** epi;
    EP_hydro** epj;
    EP_hydro* epj_send;
    int** id_epj;
    int* n_epi;
    int* n_epj;
    int nsend_epj, n_walk;

    ReadInputs(epi,n_epi,epj_send,nsend_epj,epj,id_epj,n_epj,n_walk);

    Force_hydro** force = new Force_hydro*[n_walk];
    Force_hydro** force_ref = new Force_hydro*[n_walk];
    for(int w=0;w<n_walk;w++){
      force[w] = new Force_hydro[n_epi[w]];
      force_ref[w] = new Force_hydro[n_epi[w]];
      for(int i=0;i<n_epi[w];i++){
	force[w][i].acc = 0.0;
	force[w][i].eng_dot = 0.0;
	force[w][i].dt  = std::numeric_limits<float>::lowest();
	force_ref[w][i].acc = 0.0;
	force_ref[w][i].eng_dot = 0.0;
	force_ref[w][i].dt  = std::numeric_limits<float>::lowest();
      }
    }

    // check results
    auto ref_kernel = CalcHydroForceEpEpReference(1.0,0.3,1.0,1.0);
    for(int w=0;w<n_walk;w++){
      ref_kernel(epi[w], n_epi[w], epj[w], n_epj[w], force_ref[w]);
    }
#ifdef MULTI_WALK
    initialize_CalcHydroForceEpEp(1.0,0.3,1.0,1.0);
    DispatchCalcHydroForceEpEpIndex(0, n_walk, (const EP_hydro**)epi, n_epi, (const int**)id_epj, n_epj, (const EP_hydro*)epj_send, nsend_epj, true);
    DispatchCalcHydroForceEpEpIndex(0, n_walk, (const EP_hydro**)epi, n_epi, (const int**)id_epj, n_epj, (const EP_hydro*)epj_send, nsend_epj, false);
    RetrieveCalcHydroForceEpEp(0, n_walk, n_epi, force);
#else
    auto kernel = CalcHydroForceEpEp(1.0,0.3,1.0,1.0);
    for(int w=0;w<n_walk;w++){
      kernel(epi[w], n_epi[w], epj[w], n_epj[w], force[w]);
    }
#endif
    bool isOK = true;
    constexpr double abs_error = 1e-5;
    constexpr double rel_error = 1e-5;
    for(int w=0;w<n_walk;w++){
      for(int i=0;i<n_epi[w];i++){
	double abs_accx = std::abs(force[w][i].acc.x   - force_ref[w][i].acc.x);
	double abs_accy = std::abs(force[w][i].acc.y   - force_ref[w][i].acc.y);
	double abs_accz = std::abs(force[w][i].acc.z   - force_ref[w][i].acc.z);
	double abs_eng  = std::abs(force[w][i].eng_dot - force_ref[w][i].eng_dot);
	double abs_dt   = std::abs(force[w][i].dt      - force_ref[w][i].dt);
	double rel_accx = std::abs(force_ref[w][i].acc.x)   < abs_error ? 0.0 : std::abs(abs_accx / force_ref[w][i].acc.x);
	double rel_accy = std::abs(force_ref[w][i].acc.y)   < abs_error ? 0.0 : std::abs(abs_accy / force_ref[w][i].acc.y);
	double rel_accz = std::abs(force_ref[w][i].acc.z)   < abs_error ? 0.0 : std::abs(abs_accz / force_ref[w][i].acc.z);
	double rel_eng  = std::abs(force_ref[w][i].eng_dot) < abs_error ? 0.0 : std::abs(abs_eng  / force_ref[w][i].eng_dot);
	double rel_dt   = std::abs(force_ref[w][i].dt)      < abs_error ? 0.0 : std::abs(abs_dt   / force_ref[w][i].dt);
	if(rel_accx > rel_error){
          std::cout << "acc.x: " << w << " " << i << " " << force[w][i].acc.x << " " << force_ref[w][i].acc.x << " " << abs_accx << " " << rel_accx << std::endl;
	  isOK = false;
	}
	if(rel_accy > rel_error){
          std::cout << "acc.y: " << w << " " << i << " " << force[w][i].acc.y << " " << force_ref[w][i].acc.y << " " << abs_accy << " " << rel_accy << std::endl;
	  isOK = false;
	}
	if(rel_accz > rel_error){
          std::cout << "acc.z: " << w << " " << i << " " << force[w][i].acc.z << " " << force_ref[w][i].acc.z << " " << abs_accz << " " << rel_accz << std::endl;
	  isOK = false;
	}
	if(rel_eng > rel_error){
          std::cout << "eng_dot: " << w << " " << i << " " << force[w][i].eng_dot << " " << force_ref[w][i].eng_dot << " " << abs_eng << " " << rel_eng << std::endl;
	  isOK = false;
	}
	if(rel_dt > rel_error){
          std::cout << "dt: " << w << " " << i << " " << force[w][i].dt << " " << force_ref[w][i].dt << " " << abs_dt << " " << rel_dt << std::endl;
	  isOK = false;
        }
      }
    }
    if(isOK) std::cout << "Test passed!" << std::endl;
    else     std::cout << "TEST FAILED" << std::endl;

    for(int w=0;w<n_walk;w++){
      delete[] epi[w];
      delete[] epj[w];
      delete[] id_epj[w];
      delete[] force[w];
      delete[] force_ref[w];
    }
    delete[] epi;
    delete[] epj;
    delete[] force;
    delete[] force_ref;
    delete[] epj_send;

    return 0;
}
