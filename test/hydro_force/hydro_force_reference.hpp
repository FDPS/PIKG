#include<pikg_vector.hpp>
#include<cmath>
#include<limits>
#include<chrono>

struct CalcHydroForceEpEpReference{
PIKG::F32 g1u;
PIKG::F32 cfl;
PIKG::F32 a;
PIKG::F32 H;
CalcHydroForceEpEpReference(){}
CalcHydroForceEpEpReference(PIKG::F32 g1u,PIKG::F32 cfl,PIKG::F32 a,PIKG::F32 H):g1u(g1u),cfl(cfl),a(a),H(H){}
void initialize(PIKG::F32 g1u_,PIKG::F32 cfl_,PIKG::F32 a_,PIKG::F32 H_){
g1u = g1u_;
cfl = cfl_;
a = a_;
H = H_;
}
int kernel_id = 0;
void operator()(const EP_hydro* __restrict__ epi,const int ni,const EP_hydro* __restrict__ epj,const int nj,Force_hydro* __restrict__ force,const int kernel_select = 1){
static_assert(sizeof(EP_hydro) == 80,"check consistency of EPI member variable definition between PIKG source and original source");
static_assert(sizeof(EP_hydro) == 80,"check consistency of EPJ member variable definition between PIKG source and original source");
static_assert(sizeof(Force_hydro) == 20,"check consistency of FORCE member variable definition between PIKG source and original source");
if(kernel_select>=0) kernel_id = kernel_select;
if(kernel_id == 0){
std::cout << "ni: " << ni << " nj:" << nj << std::endl;
Force_hydro* force_tmp = new Force_hydro[ni];
std::chrono::system_clock::time_point  start, end;
double min_time = std::numeric_limits<double>::max();
{ // test Kernel_I1_J1
for(int i=0;i<ni;i++) force_tmp[i] = force[i];
start = std::chrono::system_clock::now();
Kernel_I1_J1(epi,ni,epj,nj,force_tmp);
end = std::chrono::system_clock::now();
double elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count();
std::cerr << "kerel 1: " << elapsed << " ns" << std::endl;
if(min_time > elapsed){
min_time = elapsed;
kernel_id = 1;
}
}
delete[] force_tmp;
} // if(kernel_id == 0)
if(kernel_id == 1) Kernel_I1_J1(epi,ni,epj,nj,force);
} // operator() definition 
void Kernel_I1_J1(const EP_hydro* __restrict__ epi,const PIKG::S32 ni,const EP_hydro* __restrict__ epj,const PIKG::S32 nj,Force_hydro* __restrict__ force){
PIKG::S32 i;
PIKG::S32 j;
PIKG::F32  __attribute__ ((aligned(64))) xiloc_tmp_x[ni];
PIKG::F32  __attribute__ ((aligned(64))) xiloc_tmp_y[ni];
PIKG::F32  __attribute__ ((aligned(64))) xiloc_tmp_z[ni];
PIKG::F32  __attribute__ ((aligned(64))) viloc_tmp_x[ni];
PIKG::F32  __attribute__ ((aligned(64))) viloc_tmp_y[ni];
PIKG::F32  __attribute__ ((aligned(64))) viloc_tmp_z[ni];
PIKG::F32  __attribute__ ((aligned(64))) miloc_tmp[ni];
PIKG::F32  __attribute__ ((aligned(64))) uiloc_tmp[ni];
PIKG::F32  __attribute__ ((aligned(64))) hiloc_tmp[ni];
PIKG::F32  __attribute__ ((aligned(64))) hi2loc_tmp[ni];
PIKG::F32  __attribute__ ((aligned(64))) hi4loc_tmp[ni];
PIKG::F32  __attribute__ ((aligned(64))) hi5loc_tmp[ni];
PIKG::F32  __attribute__ ((aligned(64))) rhoiloc_tmp[ni];
PIKG::F32  __attribute__ ((aligned(64))) Piloc_tmp[ni];
PIKG::F32  __attribute__ ((aligned(64))) filoc_tmp[ni];
PIKG::F32  __attribute__ ((aligned(64))) ciloc_tmp[ni];
PIKG::F32  __attribute__ ((aligned(64))) BalSWiloc_tmp[ni];
PIKG::F32  __attribute__ ((aligned(64))) ailoc_tmp[ni];
PIKG::F32  __attribute__ ((aligned(64))) epsloc_tmp[ni];
PIKG::F32  __attribute__ ((aligned(64))) xjloc_tmp_x[nj];
PIKG::F32  __attribute__ ((aligned(64))) xjloc_tmp_y[nj];
PIKG::F32  __attribute__ ((aligned(64))) xjloc_tmp_z[nj];
PIKG::F32  __attribute__ ((aligned(64))) vjloc_tmp_x[nj];
PIKG::F32  __attribute__ ((aligned(64))) vjloc_tmp_y[nj];
PIKG::F32  __attribute__ ((aligned(64))) vjloc_tmp_z[nj];
PIKG::F32  __attribute__ ((aligned(64))) mjloc_tmp[nj];
PIKG::F32  __attribute__ ((aligned(64))) ujloc_tmp[nj];
PIKG::F32  __attribute__ ((aligned(64))) hjloc_tmp[nj];
PIKG::F32  __attribute__ ((aligned(64))) hj2loc_tmp[nj];
PIKG::F32  __attribute__ ((aligned(64))) hj4loc_tmp[nj];
PIKG::F32  __attribute__ ((aligned(64))) hj5loc_tmp[nj];
PIKG::F32  __attribute__ ((aligned(64))) rhojloc_tmp[nj];
PIKG::F32  __attribute__ ((aligned(64))) Pjloc_tmp[nj];
PIKG::F32  __attribute__ ((aligned(64))) fjloc_tmp[nj];
PIKG::F32  __attribute__ ((aligned(64))) cjloc_tmp[nj];
PIKG::F32  __attribute__ ((aligned(64))) BalSWjloc_tmp[nj];
PIKG::F32  __attribute__ ((aligned(64))) ajloc_tmp[nj];
for(i = 0;i < ni;++i){
xiloc_tmp_x[i] = (to_f64(a)*(epi[i].pos.x-epi[0].pos.x));
} // loop of i
for(i = 0;i < ni;++i){
xiloc_tmp_y[i] = (to_f64(a)*(epi[i].pos.y-epi[0].pos.y));
} // loop of i
for(i = 0;i < ni;++i){
xiloc_tmp_z[i] = (to_f64(a)*(epi[i].pos.z-epi[0].pos.z));
} // loop of i
for(i = 0;i < ni;++i){
viloc_tmp_x[i] = (a*(H*to_f32(epi[i].pos.x)+epi[i].vel.x));
} // loop of i
for(i = 0;i < ni;++i){
viloc_tmp_y[i] = (a*(H*to_f32(epi[i].pos.y)+epi[i].vel.y));
} // loop of i
for(i = 0;i < ni;++i){
viloc_tmp_z[i] = (a*(H*to_f32(epi[i].pos.z)+epi[i].vel.z));
} // loop of i
for(i = 0;i < ni;++i){
miloc_tmp[i] = (1.0f/((to_f32(epi[i].mass)*epi[i].eng)*g1u));
} // loop of i
for(i = 0;i < ni;++i){
uiloc_tmp[i] = (epi[i].eng*g1u);
} // loop of i
for(i = 0;i < ni;++i){
hiloc_tmp[i] = (a*epi[i].h);
} // loop of i
for(i = 0;i < ni;++i){
hi2loc_tmp[i] = (1.0f/(((a*epi[i].h)*a)*epi[i].h));
} // loop of i
for(i = 0;i < ni;++i){
hi5loc_tmp[i] = (1.0f/(((((((((a*epi[i].h)*a)*epi[i].h)*a)*epi[i].h)*a)*epi[i].h)*a)*epi[i].h));
} // loop of i
for(i = 0;i < ni;++i){
rhoiloc_tmp[i] = epi[i].dens;
} // loop of i
for(i = 0;i < ni;++i){
Piloc_tmp[i] = (1.0f/epi[i].pres);
} // loop of i
for(i = 0;i < ni;++i){
filoc_tmp[i] = epi[i].gradh;
} // loop of i
for(i = 0;i < ni;++i){
ciloc_tmp[i] = epi[i].snds;
} // loop of i
for(i = 0;i < ni;++i){
BalSWiloc_tmp[i] = epi[i].BalSW;
} // loop of i
for(i = 0;i < ni;++i){
ailoc_tmp[i] = epi[i].alpha;
} // loop of i
for(i = 0;i < ni;++i){
epsloc_tmp[i] = (((((a*epi[i].h)*0.00001f)*a)*epi[i].h)*0.00001f);
} // loop of i
for(j = 0;j < nj;++j){
xjloc_tmp_x[j] = (to_f64(a)*(epj[j].pos.x-epi[0].pos.x));
} // loop of j
for(j = 0;j < nj;++j){
xjloc_tmp_y[j] = (to_f64(a)*(epj[j].pos.y-epi[0].pos.y));
} // loop of j
for(j = 0;j < nj;++j){
xjloc_tmp_z[j] = (to_f64(a)*(epj[j].pos.z-epi[0].pos.z));
} // loop of j
for(j = 0;j < nj;++j){
vjloc_tmp_x[j] = (a*(H*to_f32(epj[j].pos.x)+epj[j].vel.x));
} // loop of j
for(j = 0;j < nj;++j){
vjloc_tmp_y[j] = (a*(H*to_f32(epj[j].pos.y)+epj[j].vel.y));
} // loop of j
for(j = 0;j < nj;++j){
vjloc_tmp_z[j] = (a*(H*to_f32(epj[j].pos.z)+epj[j].vel.z));
} // loop of j
for(j = 0;j < nj;++j){
mjloc_tmp[j] = to_f32(epj[j].mass);
} // loop of j
for(j = 0;j < nj;++j){
ujloc_tmp[j] = ((to_f32(epj[j].mass)*epj[j].eng)*g1u);
} // loop of j
for(j = 0;j < nj;++j){
hjloc_tmp[j] = (a*epj[j].h);
} // loop of j
for(j = 0;j < nj;++j){
hj2loc_tmp[j] = (1.0f/(((a*epj[j].h)*a)*epj[j].h));
} // loop of j
for(j = 0;j < nj;++j){
hj5loc_tmp[j] = (1.0f/(((((((((a*epj[j].h)*a)*epj[j].h)*a)*epj[j].h)*a)*epj[j].h)*a)*epj[j].h));
} // loop of j
for(j = 0;j < nj;++j){
rhojloc_tmp[j] = epj[j].dens;
} // loop of j
for(j = 0;j < nj;++j){
Pjloc_tmp[j] = (1.0f/epj[j].pres);
} // loop of j
for(j = 0;j < nj;++j){
fjloc_tmp[j] = epj[j].gradh;
} // loop of j
for(j = 0;j < nj;++j){
cjloc_tmp[j] = epj[j].snds;
} // loop of j
for(j = 0;j < nj;++j){
BalSWjloc_tmp[j] = epj[j].BalSW;
} // loop of j
for(j = 0;j < nj;++j){
ajloc_tmp[j] = epj[j].alpha;
} // loop of j
PIKG::F32 tab0[] = {9.58004170197862318f,7.00522516417044578f,5.35895285590360211f,4.12866206443424342f,3.17062071823030587f,2.41185311812114876f,1.807527024614142441f,1.327023708411639595f,0.947923192455528625f,0.652953821363920891f,0.4282657430174724067f,0.262367050056875637f,0.1454154330327533217f,0.0687017377634103864f,0.0242156190165320471f,0.00416773812099807155f};
PIKG::F32 tab1[] = {-(59.8108501509125188f),-(31.05813328092411915f),-(22.44668696843435463f),-(17.2486364282689537f),-(13.58530941613384657f),-(10.80687593699297624f),-(8.60901668075986312f),-(6.82439607574598537f),-(5.35129661199654011f),-(4.123594784825375975f),-(3.096333696234045568f),-(2.238156835110497202f),-(1.527171556711805323f),-(0.948829213125291506f),-(0.495512143316267556f),-(0.1702541254620237798f)};
PIKG::F32 tab2[] = {472.125747428997911f,89.8785006474572409f,49.73973974443669843f,33.684795047245919f,24.97836215430330778f,19.4838847817621894f,15.6789152631370973f,12.86998444012059337f,10.69480846022754172f,8.9445317974717133f,7.4889013736202793f,6.240554988761701914f,5.13575827666686609f,4.121869456526907f,3.146437187997371147f,2.197916376996970456f};
PIKG::F32 tab3[] = {-(2817.13182179913717f),-(230.7305434099217958f),-(88.8341996044038539f),-(47.4683723462677079f),-(29.7371881212689023f),-(20.50676367606451472f),-(15.09880337604735571f),-(11.67177812519695984f),-(9.38065185944797533f),-(7.79538433501862086f),-(6.68263935674222611f),-(5.91510404149628297f),-(5.43624235906013173f),-(5.266341660526284156f),-(5.62445972826792737f),-(8.63456778215224677f)};
for(i = 0;i < ni;++i){
PIKG::F32 BalSWiloc;

BalSWiloc = BalSWiloc_tmp[i+0];
PIKG::F32 Piloc;

Piloc = Piloc_tmp[i+0];
PIKG::F32 ailoc;

ailoc = ailoc_tmp[i+0];
PIKG::F32 ciloc;

ciloc = ciloc_tmp[i+0];
PIKG::F32 epsloc;

epsloc = epsloc_tmp[i+0];
PIKG::F32 filoc;

filoc = filoc_tmp[i+0];
PIKG::F32 hi2loc;

hi2loc = hi2loc_tmp[i+0];
PIKG::F32 hi5loc;

hi5loc = hi5loc_tmp[i+0];
PIKG::F32 miloc;

miloc = miloc_tmp[i+0];
PIKG::F32 rhoiloc;

rhoiloc = rhoiloc_tmp[i+0];
PIKG::F32 uiloc;

uiloc = uiloc_tmp[i+0];
PIKG::F32vec viloc;

viloc.x = viloc_tmp_x[i+0];
viloc.y = viloc_tmp_y[i+0];
viloc.z = viloc_tmp_z[i+0];
PIKG::F32vec xiloc;

xiloc.x = xiloc_tmp_x[i+0];
xiloc.y = xiloc_tmp_y[i+0];
xiloc.z = xiloc_tmp_z[i+0];
PIKG::F32vec FORCE_acc;

FORCE_acc.x = 0.0f;
FORCE_acc.y = 0.0f;
FORCE_acc.z = 0.0f;
PIKG::F32 FORCE_dt;

FORCE_dt = std::numeric_limits<float>::lowest();
PIKG::F32 FORCE_eng_dot;

FORCE_eng_dot = 0.0f;
for(j = 0;j < nj;++j){
PIKG::F32 BalSWjloc;

BalSWjloc = BalSWjloc_tmp[j+0];
PIKG::F32 Pjloc;

Pjloc = Pjloc_tmp[j+0];
PIKG::F32 ajloc;

ajloc = ajloc_tmp[j+0];
PIKG::F32 cjloc;

cjloc = cjloc_tmp[j+0];
PIKG::F32 fjloc;

fjloc = fjloc_tmp[j+0];
PIKG::F32 hj2loc;

hj2loc = hj2loc_tmp[j+0];
PIKG::F32 hj5loc;

hj5loc = hj5loc_tmp[j+0];
PIKG::F32 mjloc;

mjloc = mjloc_tmp[j+0];
PIKG::F32 rhojloc;

rhojloc = rhojloc_tmp[j+0];
PIKG::F32 ujloc;

ujloc = ujloc_tmp[j+0];
PIKG::F32vec vjloc;

vjloc.x = vjloc_tmp_x[j+0];
vjloc.y = vjloc_tmp_y[j+0];
vjloc.z = vjloc_tmp_z[j+0];
PIKG::F32vec xjloc;

xjloc.x = xjloc_tmp_x[j+0];
xjloc.y = xjloc_tmp_y[j+0];
xjloc.z = xjloc_tmp_z[j+0];
PIKG::F32vec dr;

PIKG::F32vec dv;

PIKG::F32 __fkg_tmp1;

PIKG::F32 __fkg_tmp0;

PIKG::F32 r2;

PIKG::F32 __fkg_tmp3;

PIKG::F32 __fkg_tmp2;

PIKG::F32 dvdr;

PIKG::F32 r_inv;

PIKG::F32 wij;

PIKG::F32 __fkg_tmp4;

PIKG::F32 v_sig;

PIKG::F32 __fkg_tmp5;

PIKG::F32 aij;

PIKG::F32 __fkg_tmp9;

PIKG::F32 __fkg_tmp8;

PIKG::F32 __fkg_tmp7;

PIKG::F32 __fkg_tmp6;

PIKG::F32 __fkg_tmp10;

PIKG::F32 AV;

PIKG::F32 u_i;

PIKG::F32 du_i;

PIKG::F32 findex_i;

PIKG::U32 index_i;

PIKG::F32 c0_i;

PIKG::F32 c1_i;

PIKG::F32 c2_i;

PIKG::F32 c3_i;

PIKG::F32 __fkg_tmp12;

PIKG::F32 __fkg_tmp11;

PIKG::F32 tmp_i;

PIKG::F32 __fkg_tmp14;

PIKG::F32 __fkg_tmp13;

PIKG::F32 gradWi;

PIKG::F32 u_j;

PIKG::F32 du_j;

PIKG::F32 findex_j;

PIKG::U32 index_j;

PIKG::F32 c0_j;

PIKG::F32 c1_j;

PIKG::F32 c2_j;

PIKG::F32 c3_j;

PIKG::F32 __fkg_tmp16;

PIKG::F32 __fkg_tmp15;

PIKG::F32 tmp_j;

PIKG::F32 __fkg_tmp18;

PIKG::F32 __fkg_tmp17;

PIKG::F32 gradWj;

PIKG::F32 __fkg_tmp19;

PIKG::F32 gradWij;

PIKG::F32 fij;

PIKG::F32 __fkg_tmp20;

PIKG::F32 fji;

PIKG::F32 tmp0;

PIKG::F32 tmp1;

PIKG::F32 __fkg_tmp21;

PIKG::F32 tmp2;

PIKG::F32 tmp3;

PIKG::F32 __fkg_tmp23;

PIKG::F32 __fkg_tmp22;

PIKG::F32 __fkg_tmp24;

PIKG::F32 tmp4;

dr.x = (xiloc.x-xjloc.x);
dr.y = (xiloc.y-xjloc.y);
dr.z = (xiloc.z-xjloc.z);
dv.x = (viloc.x-vjloc.x);
dv.y = (viloc.y-vjloc.y);
dv.z = (viloc.z-vjloc.z);
__fkg_tmp1 = (dr.x*dr.x+epsloc);
__fkg_tmp0 = (dr.y*dr.y+__fkg_tmp1);
r2 = (dr.z*dr.z+__fkg_tmp0);
__fkg_tmp3 = (dr.y*dv.y);
__fkg_tmp2 = (dr.x*dv.x+__fkg_tmp3);
dvdr = (dr.z*dv.z+__fkg_tmp2);
r_inv = rsqrt(r2);
wij = min(0.0f,(dvdr*r_inv));
__fkg_tmp4 = (ciloc+cjloc);
v_sig = (__fkg_tmp4 - 3.0f*wij);
FORCE_dt = max(FORCE_dt,v_sig);
__fkg_tmp5 = (ailoc+ajloc);
aij = (0.25f*__fkg_tmp5);
__fkg_tmp9 = -(aij);
__fkg_tmp8 = (__fkg_tmp9*v_sig);
__fkg_tmp7 = (__fkg_tmp8*wij);
__fkg_tmp6 = (__fkg_tmp7*inv((rhoiloc+rhojloc)));
__fkg_tmp10 = (BalSWiloc+BalSWjloc);
AV = (__fkg_tmp6*__fkg_tmp10);
u_i = (r2*hi2loc);
du_i = min(u_i,1.0f);
findex_i = (du_i*16.0f);
index_i = to_uint(findex_i);
index_i = min(index_i,15);
du_i = (du_i - 0.0625f*to_float(index_i));
c0_i = table(tab0,index_i);
c1_i = table(tab1,index_i);
c2_i = table(tab2,index_i);
c3_i = table(tab3,index_i);
__fkg_tmp12 = (du_i*c3_i+c2_i);
__fkg_tmp11 = (du_i*__fkg_tmp12+c1_i);
tmp_i = (du_i*__fkg_tmp11+c0_i);
__fkg_tmp14 = (tmp_i*tmp_i);
__fkg_tmp13 = -(__fkg_tmp14);
gradWi = (__fkg_tmp13*hi5loc);
u_j = (r2*hj2loc);
du_j = min(u_j,1.0f);
findex_j = (du_j*16.0f);
index_j = to_uint(findex_j);
index_j = min(index_j,15);
du_j = (du_j - 0.0625f*to_float(index_j));
c0_j = table(tab0,index_j);
c1_j = table(tab1,index_j);
c2_j = table(tab2,index_j);
c3_j = table(tab3,index_j);
__fkg_tmp16 = (du_j*c3_j+c2_j);
__fkg_tmp15 = (du_j*__fkg_tmp16+c1_j);
tmp_j = (du_j*__fkg_tmp15+c0_j);
__fkg_tmp18 = (tmp_j*tmp_j);
__fkg_tmp17 = -(__fkg_tmp18);
gradWj = (__fkg_tmp17*hj5loc);
__fkg_tmp19 = (gradWi+gradWj);
gradWij = (0.5f*__fkg_tmp19);
fij = (1.0f - filoc*inv(ujloc));
__fkg_tmp20 = (1.0f - fjloc*miloc);
fji = (__fkg_tmp20*Pjloc);
tmp0 = (uiloc*ujloc);
tmp1 = (mjloc*AV);
__fkg_tmp21 = (fij*Piloc);
tmp2 = (__fkg_tmp21*gradWi);
tmp3 = (tmp1*gradWij);
__fkg_tmp23 = (tmp0*tmp2);
__fkg_tmp22 = (0.5f*tmp3+__fkg_tmp23);
FORCE_eng_dot = (__fkg_tmp22*dvdr+FORCE_eng_dot);
__fkg_tmp24 = (fji*gradWj+tmp2);
tmp4 = (tmp0*__fkg_tmp24+tmp3);
FORCE_acc.x = (FORCE_acc.x - tmp4*dr.x);
FORCE_acc.y = (FORCE_acc.y - tmp4*dr.y);
FORCE_acc.z = (FORCE_acc.z - tmp4*dr.z);
} // loop of j

force[i+0].acc.x = (force[i+0].acc.x+FORCE_acc.x);
force[i+0].acc.y = (force[i+0].acc.y+FORCE_acc.y);
force[i+0].acc.z = (force[i+0].acc.z+FORCE_acc.z);
force[i+0].dt = max(FORCE_dt,force[i+0].dt);
force[i+0].eng_dot = (force[i+0].eng_dot+FORCE_eng_dot);
} // loop of i
} // Kernel_I1_J1 definition 
PIKG::F64 rsqrt(PIKG::F64 op){ return 1.0/std::sqrt(op); }
PIKG::F64 sqrt(PIKG::F64 op){ return std::sqrt(op); }
PIKG::F64 inv(PIKG::F64 op){ return 1.0/op; }
PIKG::F64 max(PIKG::F64 a,PIKG::F64 b){ return std::max(a,b);}
PIKG::F64 min(PIKG::F64 a,PIKG::F64 b){ return std::min(a,b);}
PIKG::F32 rsqrt(PIKG::F32 op){ return 1.f/std::sqrt(op); }
PIKG::F32 sqrt(PIKG::F32 op){ return std::sqrt(op); }
PIKG::F32 inv(PIKG::F32 op){ return 1.f/op; }
PIKG::S64 max(PIKG::S64 a,PIKG::S64 b){ return std::max(a,b);}
PIKG::S64 min(PIKG::S64 a,PIKG::S64 b){ return std::min(a,b);}
PIKG::S32 max(PIKG::S32 a,PIKG::S32 b){ return std::max(a,b);}
PIKG::S32 min(PIKG::S32 a,PIKG::S32 b){ return std::min(a,b);}
PIKG::F64 table(PIKG::F64 tab[],PIKG::S64 i){ return tab[i]; }
PIKG::F32 table(PIKG::F32 tab[],PIKG::S32 i){ return tab[i]; }
PIKG::F64 to_float(PIKG::U64 op){return (PIKG::F64)op;}
PIKG::F32 to_float(PIKG::U32 op){return (PIKG::F32)op;}
PIKG::F64 to_float(PIKG::S64 op){return (PIKG::F64)op;}
PIKG::F32 to_float(PIKG::S32 op){return (PIKG::F32)op;}
PIKG::S64   to_int(PIKG::F64 op){return (PIKG::S64)op;}
PIKG::S32   to_int(PIKG::F32 op){return (PIKG::S32)op;}
PIKG::U64  to_uint(PIKG::F64 op){return (PIKG::U64)op;}
PIKG::U32  to_uint(PIKG::F32 op){return (PIKG::U32)op;}
template<typename T> PIKG::F64 to_f64(const T& op){return (PIKG::F64)op;}
template<typename T> PIKG::F32 to_f32(const T& op){return (PIKG::F32)op;}
template<typename T> PIKG::S64 to_s64(const T& op){return (PIKG::S64)op;}
template<typename T> PIKG::S32 to_s32(const T& op){return (PIKG::S32)op;}
template<typename T> PIKG::U64 to_u64(const T& op){return (PIKG::U64)op;}
template<typename T> PIKG::U32 to_u32(const T& op){return (PIKG::U32)op;}
};// kernel functor definition 
