#ifndef VECTORWRAPPER_H
#define VECTORWRAPPER_H

#include "Lua.hh"
#include "base.h"
#include <thrust/device_ptr.h>
#include <iostream>
template<typename Container>
struct vector_trait{
};
template<>
struct vector_trait<host_vector_double>{
	static const char* const name(){
		return  "DoubleHostVector";
	}
	typedef std::false_type IsDevice;
	typedef device_vector_double ConvertableType;
};
template<>
struct vector_trait<device_vector_double>{
	static const char* const name(){
		return  "DoubleDeviceVector";
	}
	typedef std::true_type IsDevice;
	typedef host_vector_double ConvertableType;
};
template<>
struct vector_trait<host_vector_int>{
	static const char* const name(){
		return  "IntHostVector";
	}
	typedef std::false_type IsDevice;
	typedef device_vector_int ConvertableType;
};
template<>
struct vector_trait<device_vector_int>{
	static const char* const name(){
		return  "IntDeviceVector";
	}
	typedef std::true_type IsDevice;
	typedef host_vector_int ConvertableType;
};
template<class T>
void fill_thrust_vector(T& data, lua_Integer n, typename T::value_type v, std::false_type){
	data.resize(n, v);
}

void fill_n_v(device_vector_double& , int, typename device_vector_double::value_type);
void fill_n_v(device_vector_int& , int, typename device_vector_int::value_type);
template<class T>
void fill_thrust_vector(T& data, lua_Integer n, typename T::value_type v, std::true_type){
	fill_n_v(data, n, v);
}
template<typename Container> class VectorWrapper;
typedef VectorWrapper< host_vector_double > DoubleHostVector;
typedef VectorWrapper< device_vector_double> DoubleGPUVector;
typedef VectorWrapper< host_vector_int > IntHostVector;
typedef VectorWrapper< device_vector_int > IntGPUVector;

template<typename Container>
class VectorWrapper: public util::LuaClass{
	
	public:
	Container data;
	typedef typename Container::value_type value_type;
	VectorWrapper(){}
    //VectorWrapper(int size, int def):data(size, def){}
	template<typename Other>
	VectorWrapper(const VectorWrapper<Other>& v):data(v.data){}
	int size(){return data.size();}
	void fill_n(lua_Integer n, typename Container::value_type v){
		fill_thrust_vector(data, n, v, typename vector_trait<Container>::IsDevice() );
	}
	void fillValueAtIndexs(double v, IntHostVector& indexs);
	void assign(int index, typename Container::value_type v){
		data[index] = v;
	}
	void push_back(value_type v){ data.push_back(v); }
	void write_to_file(std::string name, bool append);
	void dump(){
		
		for(auto v : data) std::cout << v << ',';
		std::cout << std::endl;
	}
	
	value_type* raw_ptr(){ return thrust::raw_pointer_cast(&data[0]);}
	value_type at(lua_Integer index){return data[index];}
	static VectorWrapper* construct(){
		return new VectorWrapper();
	}
	static VectorWrapper* construct(lua_Integer size, typename Container::value_type def){
		auto v = new VectorWrapper;
		v->fill_n(size, def);
		return v;
	}
	template<typename Other>
	static VectorWrapper* construct( VectorWrapper<Other>& v){
		return new VectorWrapper(v);
	}
	void deconstruct(){
		delete this;
	}
	static void export_me(util::Lua& vm){
		vm.export_class<VectorWrapper>();
	}
	static void export_class(util::Lua& vm){
		VectorWrapper* (*construct0)() = &VectorWrapper::construct;
		vm.export_function("new0", construct0);
		
		typedef VectorWrapper< typename vector_trait<Container>::ConvertableType > OtherVector;
		VectorWrapper* (*construct1)( OtherVector&) = &VectorWrapper::construct;
		vm.export_function("new1", construct1);
		VectorWrapper* (*construct2)(lua_Integer,value_type) = &VectorWrapper::construct;
		vm.export_function("new2", construct2);
		vm.export_method("deconstruct", &VectorWrapper::deconstruct);
		vm.export_method("size", &VectorWrapper::size);
		vm.export_method("at", &VectorWrapper::at);
		vm.export_method("set", &VectorWrapper::assign);
		vm.export_method("dump", &VectorWrapper::dump);
		vm.export_method("write_to_file", &VectorWrapper::write_to_file);
		vm.export_method("fillValueAtIndexs", &VectorWrapper::fillValueAtIndexs);
		//vm.export_method("push_back", &VectorWrapper::push_back);
	}
	static const std::string class_name() {
        return vector_trait<Container>::name();
    }

    virtual const  std::string obj_class_name() const  {
        return vector_trait<Container>::name();
    }
};
namespace Helper{
	template<typename T>
		void pushHelper(lua_State* s, T& v){
		static_assert(std::is_convertible<T*, util::LuaClass*>::value,
					  "LuaClass * required!");
		util::Lua vm(s);
		vm.object(static_cast<util::LuaClass*>(new T(v)), T::class_name());
	}
	void pushHelper(lua_State* s, int v){
		lua_pushinteger(s,v);
	}
	void pushHelper(lua_State* s, double v){
		lua_pushnumber(s, v);
	}
	void pushHelper(lua_State* s, const std::string& v){
		lua_pushstring(s, v.c_str());
	}

	template<typename T>
		void map2table(lua_State* s, T& map){
		for(auto &pair : map){
			pushHelper(s, pair.first);
			pushHelper(s, pair.second);
			lua_settable(s, -3);
		}
	}
	void for_each_kv(lua_State* state, int index, std::function<void(const char*,const char*)> functor){
		lua_pushnil(state);
		while(lua_next(state, index)){
			lua_pushvalue(state,-2);
			functor(lua_tostring(state, -1), lua_tostring(state, -2));
			lua_pop(state, 2);
		}
	}

};

#endif
