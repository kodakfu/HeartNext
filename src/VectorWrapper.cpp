#include"VectorWrapper.h"
extern "C" {
#include <lualib.h>
#include <lauxlib.h>
};

#include<stdio.h>
#include<string>
#include<algorithm>
#include<sstream>
#include<map>
#include<set>
#include<iostream>
#include<fstream>
#include<iterator>
#include <arpa/inet.h>
//#include<functional>
using std::string;

int unpack4byte(string data)
{
	if(data.size() != 4){
		std::cout << "wrong size to unpack4byte " << data.size() << std::endl;
		throw("unpack4byte need 4 bytes");
	} 
	uint32_t* p = (uint32_t*)&data[0];
	uint32_t rs = ntohl(*p);
	return rs;
}
string pack4byte(int data)
{
	uint32_t  net = htonl(data);
	char* iter = (char*)&net;
	string rs(iter, iter+4);
	return rs;
}
template<typename T>
void doWriteFile(std::ofstream& file, T& data, std::false_type)
{
	file.write( (char*)&data[0], sizeof(data[0]) * data.size() );
}
template<typename T>
void doWriteFile(std::ofstream& file, T& data, std::true_type)
{
	typename vector_trait<T>::ConvertableType host = data;
	file.write( (char*)&host[0], sizeof(host[0]) * host.size() );
}
template<typename Container>
void VectorWrapper<Container>::write_to_file(string name, bool append)
{
	using std::ofstream;
	std::ios_base::openmode mode = ofstream::out;
	if(append) mode |= ofstream::app;
	ofstream file(name.c_str(), mode);
	doWriteFile(file, data, typename vector_trait<Container>::IsDevice());
}
template<typename T>
void VectorWrapper<T>::fillValueAtIndexs(double v, IntHostVector& indexs)
{
	for(int i : indexs.data){
		data[i] = v;
	}
}
void passByRef(DoubleHostVector& a){
	a.data[0] = 100;
	printf("in passByRef\n");
	a.dump();
	//printf("%d, %d\n", a->size(), b);
}
template<typename Device, typename Host>
Device* Host2Device(Host& host){
	return new Device(host);
}
template<typename T, typename P>
std::pair<int,int> PosNegSurroundIndex(const T& rawData, const T& gIndex, int index, int stride,
									   int lower, int upper, int nilVoxel, const P& ignore)
{
	int pos = index + stride;
	int neg = index - stride;
	int vpos = nilVoxel, vneg = nilVoxel;
	if(pos >= lower && pos < upper) {vpos = rawData[pos]; pos = gIndex[pos];}
	else pos = -1;
	if(neg >= lower && neg < upper) {vneg = rawData[neg]; neg = gIndex[neg];}
	else neg = -1;
	index = gIndex[index];
	if(vpos == nilVoxel || ignore.count(vpos) == 1) pos = -1;
	if(vneg == nilVoxel || ignore.count(vneg) == 1) neg = -1;
	if(pos < 0 && neg >= 0) pos = neg;
	if(pos < 0) pos = index;
	if(neg < 0) neg = pos;
	return std::make_pair(pos, neg);
}

std::map<int, std::set<int> > GetIgnoreTable(lua_State* state, int index){
	std::map<int, std::set<int>> rs;
	Helper::for_each_kv(state, index, [&rs](const char* k, const char* v){
			//std::cout << k << "->" << v << std::endl;
			int ik = std::stoi(k);
			string sv(v);
			std::istringstream stream(sv);
			std::set<int> set;
			int iv;
			printf("ignore table for %d:",ik);
			while(stream>>iv){ set.insert(iv); printf("%d,", iv);}
			printf("\n");
			rs[ik] = set;
		});
	return rs;
}
template<typename T>
void dumpMap(std::map<int, T>& map){
	using std::cout;
	for(auto &pair: map){
		cout << pair.first << ":";
		std::for_each(pair.second.begin(), pair.second.end(), [](int v){cout << v << ',';});
		cout << std::endl;
	}
}
template<typename T, typename I,typename V>
bool judge(T& srd, I& tofile, V& volume)
{
	return volume[tofile[srd.first]] > 15 || volume[tofile[srd.second]] > 15;
}
// dimx, dimy, dimz, data(host int), nilVoxel, ignore table:{1:{12,3}, 2:{3,3}}
extern "C" int parse_tissue(lua_State* state){
	util::Lua vm(state);
	if(lua_gettop(state) != 6)
		luaL_error(state, " wrong param number!");
	const int dimx = vm.tonumber(1);
	const int dimy = vm.tonumber(2);
	const int dimz = vm.tonumber(3);
	IntHostVector& data = *(IntHostVector*)vm.object(4);
	host_vector_int& rawData = data.data;
	const int nilVoxel = vm.tonumber(5);
	auto ignore =  GetIgnoreTable(state, 6);
	//dumpMap(ignore);
	std::map<string, IntHostVector> surround;
	IntHostVector& xplus = surround["xpos"];
	IntHostVector& xminus = surround["xneg"];
	IntHostVector& yplus = surround["ypos"];
	IntHostVector& yminus = surround["yneg"];
	IntHostVector& zplus = surround["zpos"];
	IntHostVector& zminus = surround["zneg"];
	IntHostVector globalIndex, indexInTissueFile;
	
	int totalCellNumber = 0;
	std::map<int, IntHostVector> local2global;
	
	if(dimx * dimy * dimz != rawData.size())
		luaL_error(state, "wrong dim");
	std::map<int, int> voxelNumbers;
	std::cout << "start first pass\n";
	for(auto v : rawData) voxelNumbers[v]+=1;
	globalIndex.data.reserve(rawData.size());
	indexInTissueFile.data.reserve(rawData.size() - voxelNumbers[nilVoxel]);
	for(auto kv : voxelNumbers){
		if(kv.first == nilVoxel) continue;
		local2global[kv.first].data.reserve(kv.second);
	}
	for(int rawIndex = 0; rawIndex != rawData.size(); ++rawIndex){
		int voxel = rawData[rawIndex];
		if(voxel == nilVoxel){
			globalIndex.data.push_back(-1); continue;
		}
		local2global[voxel].data.push_back(totalCellNumber);
		globalIndex.data.push_back(totalCellNumber);
		indexInTissueFile.data.push_back(rawIndex);
		++totalCellNumber;
	}
	std::cout << "first pass done \n";
	//int isolated = 0;
	for(int z = 0; z < dimz; ++z)
		for(int y = 0; y < dimy; ++y)
			for(int x=0; x < dimx; ++x){
				int rawIndex = x + y * dimx + z * dimx * dimy;
				int voxel = rawData[rawIndex];
				if(voxel == nilVoxel) continue;
				if(ignore.count(voxel) == 0) luaL_error(state, "missing ignore");
				auto& theIgnore = ignore[voxel];
				int planSize = dimx * dimy;
				int planStart = z * planSize;
				int rowStart = planStart + y * dimx;
				auto pn = PosNegSurroundIndex(rawData, globalIndex.data, rawIndex, 1,
											  rowStart, rowStart + dimx,  nilVoxel, theIgnore);
				xplus.data.push_back(pn.first);
				xminus.data.push_back(pn.second);
				
				pn = PosNegSurroundIndex(rawData, globalIndex.data, rawIndex, dimx,
										 planStart, planStart + planSize, nilVoxel, theIgnore);
				yplus.data.push_back(pn.first);
				yminus.data.push_back(pn.second);
				
				pn = PosNegSurroundIndex(rawData, globalIndex.data, rawIndex, dimx * dimy,
										 0, planSize * dimz, nilVoxel, theIgnore);
				zplus.data.push_back(pn.first);
				zminus.data.push_back(pn.second);
			}
	//std::cout << "isolated:"<<isolated;
	std::cout << "second pass done \n";
	vm.table();
	vm.number(totalCellNumber);
	vm.save("totalCellNumber", -2);
	Helper::map2table(state, local2global);
	vm.table();
	Helper::map2table(state, surround);
	vm.save("surround", -2);
	Helper::pushHelper(state, globalIndex);
	vm.save("globalIndex", -2);
	Helper::pushHelper(state, indexInTissueFile);
	vm.save("indexInTissueFile", -2);
	Helper::pushHelper(state, dimx); vm.save("dimx", -2);
	Helper::pushHelper(state, dimy); vm.save("dimy", -2);
	Helper::pushHelper(state, dimz); vm.save("dimz", -2);
	return 1;
}
void do_MapDoubleDataByIndex(double* src, int* index, double* dst, int number, bool get);
void MapDoubleDataByIndex( DoubleGPUVector& src, IntGPUVector& index, DoubleGPUVector& dst, bool get){
	if(src.size() != index.size()){
		printf("damn! src.size() != index.size()\n");
		return;
	}
	const int number = src.size();
	do_MapDoubleDataByIndex(src.raw_ptr(), index.raw_ptr(), dst.raw_ptr(), number, get);
}
void copyDoubleDeviceVector(DoubleGPUVector& src, DoubleGPUVector& dst)
{
	dst.data = src.data;
}
void wait_gpu_work_done(int index)
{
	cudaSetDevice(index);
	cudaDeviceSynchronize();
}
void select_GPU(int index){
	cudaSetDevice(index);
}
int GPU_number(){
	int rs = 0;
	cudaGetDeviceCount(&rs);
	return rs;
}
void do_CalcVolt(double* new_volt, double* volt, double* current, int* xpos, int* xneg,
				 int* ypos, int* yneg, int* zpos, int* zneg,
				 double DD, double dt, int number);
int CalcVolt(lua_State* s)
{
	util::Lua vm(s);
	// volt obj, current obj, 1 table contain surround, dd, dt
	auto& volt = *dynamic_cast<DoubleGPUVector*>( vm.object(1) ); 
	auto& current = *dynamic_cast<DoubleGPUVector*>( vm.object(2) );
	vm.load("xpos", 3); auto& xpos = *dynamic_cast<IntGPUVector*>( vm.object(-1) ); vm.pop();
	vm.load("xneg", 3); auto& xneg = *dynamic_cast<IntGPUVector*>( vm.object(-1) ); vm.pop();
	vm.load("ypos", 3); auto& ypos = *dynamic_cast<IntGPUVector*>( vm.object(-1) ); vm.pop();
	vm.load("yneg", 3); auto& yneg = *dynamic_cast<IntGPUVector*>( vm.object(-1) ); vm.pop();
	vm.load("zpos", 3); auto& zpos = *dynamic_cast<IntGPUVector*>( vm.object(-1) ); vm.pop();
	vm.load("zneg", 3); auto& zneg = *dynamic_cast<IntGPUVector*>( vm.object(-1) ); vm.pop();
	double DD = vm.tonumber(4);
	double dt = vm.tonumber(5);
	auto& new_volt = *dynamic_cast<DoubleGPUVector*>( vm.object(6) );
	const int number = volt.size();
	//printf("dump surround in CalcVolt\n");
	//xpos.dump(); xneg.dump(); ypos.dump(); yneg.dump(); zpos.dump(); zneg.dump();
	if(number != current.size() || number != new_volt.size() ){
		printf("damn! size dismatch\n"); return 0;
	}
	do_CalcVolt(new_volt.raw_ptr(), volt.raw_ptr(), current.raw_ptr(),
				xpos.raw_ptr(), xneg.raw_ptr(),
				ypos.raw_ptr(), yneg.raw_ptr(),
				zpos.raw_ptr(), zneg.raw_ptr(),
				DD, dt, number);
	return 0;
}
static void LoadIntFromByteFile(string name, IntHostVector& dst)
{
	using std::ifstream;
	ifstream file(name);
	if(!file){
		std::cout << "open " << name << " fail\n";
		return;
	}
	file.seekg(0, file.end);
	size_t length = static_cast<size_t>(file.tellg());
	file.seekg(0, file.beg);
	if(length > 0){
		using std::vector;
		using std::istream_iterator;
		vector<char> bytes;
		bytes.resize(length);
		dst.data.resize(bytes.size());
		file.read( &bytes[0], length);
		std::copy(bytes.begin(), bytes.end(), dst.data.begin());
	}
	
}
static void LoadDoubleFromFile(string name, DoubleGPUVector& dst)
{
	using std::ifstream;
	ifstream file(name);
	if(!file){
		std::cout << "open " << name << " fail\n";
		return;
	}
	file.seekg(0, file.end);
	size_t bytesNumber = static_cast<size_t>(file.tellg());
	file.seekg(0, file.beg);
	size_t n = bytesNumber / sizeof(double);
	if(n > 0){
		host_vector_double buffer(n);
		file.read( (char*)&buffer[0], bytesNumber);
		dst.data = buffer;
	}
	else{
		std::cout << "wrong double number " << n << " from file " << name << std::endl;
		return;
	}
}
extern "C" int load(lua_State* state){
	static util::Lua vm(state);
	lua_pushcfunction(state, parse_tissue);
	lua_setglobal(state, "parse_tissue");
	lua_pushcfunction(state, CalcVolt);
	lua_setglobal(state, "CalcVolt");
	DoubleHostVector::export_me(vm);
	DoubleGPUVector::export_me(vm);
	IntHostVector::export_me(vm);
	IntGPUVector::export_me(vm);
	vm.export_function("passByRef", passByRef);
	DoubleGPUVector* (*doubleH2D)(DoubleHostVector&) = Host2Device;
	vm.export_function("doubleH2D", doubleH2D);
	DoubleHostVector* (*doubleD2H)(DoubleGPUVector&) = Host2Device;
	vm.export_function("doubleD2H", doubleD2H);
	vm.export_function("MapDoubleDataByIndex", MapDoubleDataByIndex);
	vm.export_function("LoadIntFromByteFile", LoadIntFromByteFile);
	vm.export_function("LoadDoubleFromFile", LoadDoubleFromFile);
	vm.export_function("unpack4byte", unpack4byte);
	vm.export_function("pack4byte", pack4byte);
	vm.export_function("select_GPU", select_GPU);
	vm.export_function("GPU_number", GPU_number);
	vm.export_function("wait_gpu_work_done", wait_gpu_work_done);
	vm.export_function("copyDoubleDeviceVector", copyDoubleDeviceVector);
	return 0;
}
