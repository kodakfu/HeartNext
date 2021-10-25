#include"Lua.hh"
extern "C" int load(lua_State*);
extern "C" int load06Model(lua_State*);
extern "C" int loadPurkinjeModel(lua_State*);
extern "C" int loadBlockKrCalModel(lua_State* s);
int main(int argc, char** argv)
{
	util::Lua vm;
	lua_State* s = vm.getState();
	load(vm.getState());
	load06Model(vm.getState());
	loadPurkinjeModel(vm.getState());
	loadBlockKrCalModel(vm.getState());
	vm.table();
	for(int i = 2; i < argc; ++i){
		lua_pushinteger(s, i - 1);
		lua_pushstring(s, argv[i]);
		lua_settable(s, -3);
	}
	lua_setglobal(s, "arg");
	vm.file(argv[1]);
	return 0;
}
