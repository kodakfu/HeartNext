//MyEmail:xsjshr1108@163.com
#include "VectorWrapper.cuh"
#include<string>
#include<map>
#include<tuple>
extern "C"{
	#include<lua.h>
	#include <lualib.h>
    #include <lauxlib.h>
}
#include<stdio.h>

using std::string;
using std::map;
static map<string, double> ParamSpec={
	{"Cai" , 0.00007}
	,{"CaSR" , 1.3}
	,{"CaSS" , 0.00007}
	,{"Nai" , 7.67}
	,{"Ki" , 138.3}
	,{"sm" , 0.0}
	,{"sh" , 0.75}
	,{"sj" , 0.75}
	,{"sxr1" , 0.0}
	,{"sxr2" , 1}
	,{"sxs" , 0.0}
	,{"sr" , 0.0}
	,{"ss" , 1.0}
	,{"sd" , 0.0}
	,{"sf" , 1.0}
	,{"sf2" , 1.0}
	,{"sfcass" , 1.0}
	,{"sOO" , 0.0}
	,{"sRR" , 1}
	,{"sml" , 0.0}
	,{"shl" , 0.75}
};
using namespace thrust;

void call_kernel06(double* Cai_buffer, double* CaSR_buffer, double* CaSS_buffer, double* Nai_buffer,
				   double* Ki_buffer,  double* sm_buffer, double* sh_buffer, double* sj_buffer,
				   double* sxr1_buffer, double* sxr2_buffer, double* sxs_buffer, double* sr_buffer,
				   double* ss_buffer, double* sd_buffer, double* sf_buffer, double* sf2_buffer,
				   double* sfcass_buffer, double* sOO_buffer, double* sRR_buffer, double* sml_buffer,
				   double* shl_buffer, double* volt_buffer, double* istim_buffer, double* itot_buffer,
				   int* type_buffer, int cellNumber, int threadCount);

void UpdateVentricle06model_ischemia_1a_1b_infarct
(
 DoubleGPUVector& Cai_buffer, DoubleGPUVector& CaSR_buffer, DoubleGPUVector& CaSS_buffer, DoubleGPUVector& Nai_buffer,
 DoubleGPUVector& Ki_buffer, DoubleGPUVector& sm_buffer, DoubleGPUVector& sh_buffer, DoubleGPUVector& sj_buffer,
 DoubleGPUVector& sxr1_buffer, DoubleGPUVector& sxr2_buffer, DoubleGPUVector& sxs_buffer, DoubleGPUVector& sr_buffer,
 DoubleGPUVector& ss_buffer, DoubleGPUVector& sd_buffer, DoubleGPUVector& sf_buffer, DoubleGPUVector& sf2_buffer,
 DoubleGPUVector& sfcass_buffer, DoubleGPUVector& sOO_buffer, DoubleGPUVector& sRR_buffer, DoubleGPUVector& sml_buffer,
 DoubleGPUVector& shl_buffer, DoubleGPUVector& volt_buffer, DoubleGPUVector& istim_buffer, DoubleGPUVector& itot_buffer,
 IntGPUVector& type_buffer, int cellNumber, int threadCount
)
{
	call_kernel06( Cai_buffer.raw_ptr(),  CaSR_buffer.raw_ptr(),  CaSS_buffer.raw_ptr(),  Nai_buffer.raw_ptr(),
				   Ki_buffer.raw_ptr(),   sm_buffer.raw_ptr(),  sh_buffer.raw_ptr(),  sj_buffer.raw_ptr(),
				   sxr1_buffer.raw_ptr(),  sxr2_buffer.raw_ptr(), sxs_buffer.raw_ptr(),  sr_buffer.raw_ptr(),
				   ss_buffer.raw_ptr(),  sd_buffer.raw_ptr(),  sf_buffer.raw_ptr(),  sf2_buffer.raw_ptr(),
				   sfcass_buffer.raw_ptr(),  sOO_buffer.raw_ptr(),  sRR_buffer.raw_ptr(),  sml_buffer.raw_ptr(),
				   shl_buffer.raw_ptr(),  volt_buffer.raw_ptr(),  istim_buffer.raw_ptr(),  itot_buffer.raw_ptr(),
				   type_buffer.raw_ptr(), cellNumber, threadCount);
}
static int return06ParamSpec_ischemia_1a_1b_infarct(lua_State* s){
	lua_newtable(s);
	Helper::map2table(s, ParamSpec);
	return 1;
}
static int UpdateVentricleCurrent_06model_ischemia_1a_1b_infarct(lua_State* s)
{
	// 1 table contain {"Cai"->{}, "CaSR"->{}, ...}
	if(lua_gettop(s) < 1 || !lua_istable(s, -1)){
		luaL_error(s, " wrong param number or param not a table!");
		return 0;
	}
	util::Lua vm(s);
	vm.load("Cai"); auto& Cai_buffer = *dynamic_cast<DoubleGPUVector*>( vm.object(-1) ); vm.pop();
	vm.load("CaSR"); auto& CaSR_buffer = *dynamic_cast<DoubleGPUVector*>( vm.object(-1) ); vm.pop();
	vm.load("CaSS"); auto& CaSS_buffer = *dynamic_cast<DoubleGPUVector*>( vm.object(-1) ); vm.pop();
	vm.load("Nai"); auto& Nai_buffer = *dynamic_cast<DoubleGPUVector*>( vm.object(-1) ); vm.pop();
	vm.load("Ki"); auto& Ki_buffer = *dynamic_cast<DoubleGPUVector*>( vm.object(-1) ); vm.pop();
	vm.load("sm"); auto& sm_buffer = *dynamic_cast<DoubleGPUVector*>( vm.object(-1) ); vm.pop();
	vm.load("sh"); auto& sh_buffer = *dynamic_cast<DoubleGPUVector*>( vm.object(-1) ); vm.pop();
	vm.load("sj"); auto& sj_buffer = *dynamic_cast<DoubleGPUVector*>( vm.object(-1) ); vm.pop();
	vm.load("sxr1"); auto& sxr1_buffer = *dynamic_cast<DoubleGPUVector*>( vm.object(-1) ); vm.pop();
	vm.load("sxr2"); auto& sxr2_buffer = *dynamic_cast<DoubleGPUVector*>( vm.object(-1) ); vm.pop();
	vm.load("sxs"); auto& sxs_buffer = *dynamic_cast<DoubleGPUVector*>( vm.object(-1) ); vm.pop();
	vm.load("sr"); auto& sr_buffer = *dynamic_cast<DoubleGPUVector*>( vm.object(-1) ); vm.pop();
	vm.load("ss"); auto& ss_buffer = *dynamic_cast<DoubleGPUVector*>( vm.object(-1) ); vm.pop();
	vm.load("sd"); auto& sd_buffer = *dynamic_cast<DoubleGPUVector*>( vm.object(-1) ); vm.pop();
	vm.load("sf"); auto& sf_buffer = *dynamic_cast<DoubleGPUVector*>( vm.object(-1) ); vm.pop();
	vm.load("sf2"); auto& sf2_buffer = *dynamic_cast<DoubleGPUVector*>( vm.object(-1) ); vm.pop();
	vm.load("sfcass"); auto& sfcass_buffer = *dynamic_cast<DoubleGPUVector*>( vm.object(-1) ); vm.pop();
	vm.load("sOO"); auto& sOO_buffer = *dynamic_cast<DoubleGPUVector*>( vm.object(-1) ); vm.pop();
	vm.load("sRR"); auto& sRR_buffer = *dynamic_cast<DoubleGPUVector*>( vm.object(-1) ); vm.pop();
	vm.load("sml"); auto& sml_buffer = *dynamic_cast<DoubleGPUVector*>( vm.object(-1) ); vm.pop();
	vm.load("shl"); auto& shl_buffer = *dynamic_cast<DoubleGPUVector*>( vm.object(-1) ); vm.pop();
	vm.load("volt"); auto& volt_buffer = *dynamic_cast<DoubleGPUVector*>( vm.object(-1) ); vm.pop();
	vm.load("istim"); auto& istim_buffer = *dynamic_cast<DoubleGPUVector*>( vm.object(-1) ); vm.pop();
	vm.load("itot"); auto& itot_buffer = *dynamic_cast<DoubleGPUVector*>( vm.object(-1) ); vm.pop();
	vm.load("type"); auto& type_buffer = *dynamic_cast<IntGPUVector*>( vm.object(-1) ); vm.pop();
	vm.load("cellNumber"); int cellNumber = lua_tointeger(s, -1); vm.pop();
	int threadCount = lua_tointeger(s, 1);
	//*
	call_kernel06( Cai_buffer.raw_ptr(),  CaSR_buffer.raw_ptr(),  CaSS_buffer.raw_ptr(),  Nai_buffer.raw_ptr(),
				   Ki_buffer.raw_ptr(),   sm_buffer.raw_ptr(),  sh_buffer.raw_ptr(),  sj_buffer.raw_ptr(),
				   sxr1_buffer.raw_ptr(),  sxr2_buffer.raw_ptr(), sxs_buffer.raw_ptr(),  sr_buffer.raw_ptr(),
				   ss_buffer.raw_ptr(),  sd_buffer.raw_ptr(),  sf_buffer.raw_ptr(),  sf2_buffer.raw_ptr(),
				   sfcass_buffer.raw_ptr(),  sOO_buffer.raw_ptr(),  sRR_buffer.raw_ptr(),  sml_buffer.raw_ptr(),
				   shl_buffer.raw_ptr(),  volt_buffer.raw_ptr(),  istim_buffer.raw_ptr(),  itot_buffer.raw_ptr(),
				   type_buffer.raw_ptr(), cellNumber, threadCount);
	//*/
	return 0;
}
extern "C" int load06model_ischemia_1a_1b_infarct(lua_State* s){
	static util::Lua vm(s);
	lua_pushcfunction(s, return06ParamSpec_ischemia_1a_1b_infarct);
	lua_setglobal(s, "return06ParamSpec_ischemia_1a_1b_infarct");
	lua_pushcfunction(s, UpdateVentricleCurrent_06model_ischemia_1a_1b_infarct);
	lua_setglobal(s, "UpdateVentricleCurrent_06model_ischemia_1a_1b_infarct");
	//vm.export_function("UpdateVentricle06Model", UpdateVentricle06Model);
	return 0;
}
