#!/usr/local/bin/lua
package.path = '/usr/share/lua/5.1/?.lua;' .. package.path
package.cpath = '/usr/lib/x86_64-linux-gnu/lua/5.1/?.so;' .. package.cpath
dofile('core.lua')
global_gpu_index = 0
json = require('cjson')
socket = require('socket')
server = assert(socket.bind('127.0.0.1', 3333))
--server:listen()
client = server:accept()
client:settimeout(nil, 'b')
function receiveCheck(sock, n)
   local data, error_msg = sock:receive(n)
   if not data then error(error_msg) end
   print( 'receiveCheck, ', string.len(data) )
   return data
end
function sendCheck(sock, data)
   local send_n, error_msg, last_i = sock:send(data)
   if send_n ~= string.len(data) then
	  error(error_msg)
   end
end
function wrapData(data)
   local n = string.len(data)
   --if n > 512 then error('data too long') end
   return pack4byte(n) .. data
end
function reparse_fill_surround(ignoreTable)
   after = parse_tissue(tissueSpec.dimx, tissueSpec.dimy, tissueSpec.dimz,
							 tissueSpec.rawVoxels, tissueSpec.nilVoxel, ignoreTable)
   for name, v in pairs(after.surround) do
	  Gsurround[name]:deconstruct()
	  Gsurround[name] = IntDeviceVector.new1(v)
   end
end
function openTissue(fname)
   tissueSpec = dofile(fname)
   print('parse_tissue')
   allTissues = parse_tissue(tissueSpec.dimx, tissueSpec.dimy, tissueSpec.dimz,
							 tissueSpec.rawVoxels, tissueSpec.nilVoxel, tissueSpec.ignoreTable)
   print("parse_tissue done")
   select_GPU(global_gpu_index)
   Gvolt = DoubleDeviceVector.new2(allTissues.totalCellNumber, 0)
   GNewVolt = DoubleDeviceVector.new2(allTissues.totalCellNumber, 0)
   Gcurrent = DoubleDeviceVector.new2(allTissues.totalCellNumber, 0)

   VoxelVolt = DoubleDeviceVector.new2(allTissues.dimx * allTissues.dimy * allTissues.dimz, -90)
   VoxelVoltIndex = IntDeviceVector.new1(allTissues.indexInTissueFile)
   Gsurround = {}
   for name, v in pairs(allTissues.surround) do
	  Gsurround[name] = IntDeviceVector.new1(v)
   end
   DD = tissueSpec.DD
   dt = 0.02
end
function initCellData()
   select_GPU_for_cell(allTissues, tissueSpec.allVoxels)
   print("begin fill cell parameters")
   if tissueSpec.VentricleType and tissueSpec.VentricleVoxel then
	  print('fill ventricle')
	  fillCellUse06VentricleModel(allTissues, tissueSpec.VentricleVoxel, tissueSpec.VentricleType)
	  --fillCellUseBlockKrCalVentricleModel(allTissues, tissueSpec.VentricleVoxel, tissueSpec.VentricleType)
   end
   if tissueSpec.PurkinjeVoxel then
	  fillCellUsePurkinjeModel(allTissues, tissueSpec.PurkinjeVoxel)
   end
   --tissueSpec.Init(allTissues)
   print("done")
   -- the global data struct
   print("construct global data")
   select_GPU(global_gpu_index)
   for i, voxel in ipairs(tissueSpec.allVoxels) do
	  if allTissues[voxel] then
		 print('map', voxel, 'volt')
		 local var =  allTissues[voxel].params
		 local src = var.volt
		 if var.volt_in_global_gpu then
			src = var.volt_in_global_gpu
		 end
		 MapDoubleDataByIndex(src, var.local2global, Gvolt, false)
	  end
   end
   
   print("done")
   

end
function travelNeedSerializParam(dirName, cb)
   for _, voxel in ipairs(tissueSpec.allVoxels) do
	  if allTissues[voxel] then
		 local params = allTissues[voxel].params
		 local nameTemplat = dirName..'%d.%s'
		 for k,v in pairs(params) do
			if type(v) == 'table' and v.needSerialize then
			   local dstName = string.format(nameTemplat, voxel, k)
			   cb(dstName, v)
			end
		 end
	  end
   end
end
function saveSnap(dirName)
   Gvolt:write_to_file(dirName..'Gvolt', false)
   GNewVolt:write_to_file(dirName..'GNewVolt', false)
   Gcurrent:write_to_file(dirName..'Gcurrent', false)
   travelNeedSerializParam(dirName,
						   function(saveName, v)
							  print('save', saveName)
							  v:write_to_file(saveName, false)
						   end
   )
end
function loadSnap(dirName)
   LoadDoubleFromFile(dirName..'Gvolt', Gvolt)
   LoadDoubleFromFile(dirName..'GNewVolt', GNewVolt)
   LoadDoubleFromFile(dirName..'Gcurrent', Gcurrent)
   travelNeedSerializParam(dirName,
						   function(loadName, v)
							  print('load', loadName)
							  LoadDoubleFromFile(loadName, v)
						   end
   )
end
function step(istepStart, istepEnd)
   --print('step')
   for istep = istepStart, istepEnd do
	  for _, info in ipairs(allTissues.stimulate) do
		 local dst = allTissues[info.voxel].params.istim
		 if istep == math.floor(info.startTime/dt) then
			dst:fillValueAtIndexs(info.stimulate, info.index)
			local msg = string.format("stimulate %d at %f cells %d", info.voxel, info.startTime, info.index:size())
			print(msg)
		 elseif istep == math.floor(info.endTime/dt) then
			dst:fillValueAtIndexs(0, info.index)
		 end
	  end
	  for _, voxel in ipairs(tissueSpec.allVoxels) do
		 if allTissues[voxel] then
			select_GPU(allTissues[voxel].gpu_index)
			local params = allTissues[voxel].params
			params.updator(352, params)
		 end
	  end
	  wait_gpu_work_done(0)
	  wait_gpu_work_done(1)
	  for _, voxel in ipairs(tissueSpec.allVoxels) do
		 tissue = allTissues[voxel]
		 if tissue then
			local params = tissue.params
			select_GPU(global_gpu_index)
			local srcItot = params.itot
			if params.itot_in_global_gpu then
			   copyDoubleDeviceVector(params.itot, params.itot_in_global_gpu)
			   srcItot = params.itot_in_global_gpu
			end
			MapDoubleDataByIndex(srcItot, params.local2global, Gcurrent, false)
		 end
	  end
	  select_GPU(global_gpu_index)
	  CalcVolt(Gvolt, Gcurrent, Gsurround, DD, dt, GNewVolt)
	  for _, voxel in ipairs(tissueSpec.allVoxels) do
		 if allTissues[voxel] then
			local params = allTissues[voxel].params
			local dstVolt = params.volt
			if params.volt_in_global_gpu then
			   dstVolt = params.volt_in_global_gpu
			end
			select_GPU(global_gpu_index)
			MapDoubleDataByIndex( dstVolt, params.local2global,  GNewVolt, true)
			if params.volt_in_global_gpu then
			   copyDoubleDeviceVector(params.volt_in_global_gpu, params.volt )
			end
		 end
	  end
	  Gvolt, GNewVolt = GNewVolt, Gvolt
   end
   
end
function write_volt(fname)
   --print('write_volt:'..fname)
   MapDoubleDataByIndex(Gvolt, VoxelVoltIndex, VoxelVolt, false)
   VoxelVolt:write_to_file(fname , false)
end
function dump_table(data)
   print('{')
   for k,v in pairs(data) do
	  if type(v) == 'table' then
		 print(k)
		 dump_table(v)
	  else print(k,v) end
   end
   print('}')
end
function pointInBalls(balls, x, y, z)
   for _, ball in ipairs(balls) do
   --local ball = balls[2]
	  local local_x = x - ball[1]
	  local local_y = y - ball[2]
	  local local_z = z - ball[3]
	  local r = ball[4]
	  if local_x * local_x + local_y * local_y + local_z * local_z < r * r then return true end
   end
   return false
end
function setRunTimeParams(caseSpec)
   print('get caseSpec:'..caseSpec)
   caseSpec = json.decode(caseSpec)
   --DD = caseSpec.DD
   allTissues.stimulate = {}
   for _, v in ipairs(caseSpec.stimulates) do
	  print('regions size',table.maxn(v.regions) )
	  local a = makeStimulate(allTissues, v.voxel, v.strength, v.startTime, v.endTime,
							  function(x,y,z)
								 return pointInBalls(v.regions, x, y, z)
							  end
	  )
	  print('stimulate number',a.index:size())
	  table.insert(allTissues.stimulate, a)
   end
end
function initAllParams()
   print('call initAllParams')
   initCellData()
end
while 1 do
   msg_size = unpack4byte(receiveCheck(client, 4))
   command = receiveCheck(client, msg_size)
   print(msg_size, command)
   result = assert(loadstring(command))()
   --[[]
   result = 'done'
   sleep = math.random(1, 4)
   oscmd = 'sleep '..tostring(sleep)
   print('execute', oscmd)
   os.execute(oscmd)
   --]]
   data = wrapData(tostring(result))
   sendCheck(client, data)
end
   
