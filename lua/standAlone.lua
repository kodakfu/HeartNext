dofile('../lua/core.lua')

if not arg[1] then
   print("not given input file")
   return
end
print('input file:', arg[1])
global_gpu_index = 0
tissueSpec = dofile(arg[1])
print('parse_tissue')
allTissues = parse_tissue(tissueSpec.dimx, tissueSpec.dimy, tissueSpec.dimz,
						  tissueSpec.rawVoxels, tissueSpec.nilVoxel, tissueSpec.ignoreTable)
select_GPU_for_cell(allTissues, tissueSpec.allVoxels)
print("parse_tissue done")
print("begin fill cell parameters")
if tissueSpec.VentricleType and tissueSpec.VentricleVoxel then
   print('fill ventricle')
   fillCellUse06VentricleModel(allTissues, tissueSpec.VentricleVoxel, tissueSpec.VentricleType)
end
if tissueSpec.PurkinjeVoxel then
   fillCellUsePurkinjeModel(allTissues, tissueSpec.PurkinjeVoxel)
end
tissueSpec.Init(allTissues)
print("done")

print("construct global data")

select_GPU(global_gpu_index)
Gvolt = DoubleDeviceVector.new2(allTissues.totalCellNumber, 0)
for i, voxel in ipairs(tissueSpec.allVoxels) do
   if allTissues[voxel] then
	  local var =  allTissues[voxel].params
	  local src = var.volt
	  if var.volt_in_global_gpu then
		 src = var.volt_in_global_gpu
	  end
	  MapDoubleDataByIndex(src, var.local2global, Gvolt, false)
   end
end
GNewVolt = DoubleDeviceVector.new2(allTissues.totalCellNumber, 0)
Gcurrent = DoubleDeviceVector.new2(allTissues.totalCellNumber, 0)

VoxelVolt = DoubleDeviceVector.new2(allTissues.dimx * allTissues.dimy * allTissues.dimz, -100)
VoxelVoltIndex = IntDeviceVector.new1(allTissues.indexInTissueFile)
Gsurround = {}
for name, v in pairs(allTissues.surround) do
   Gsurround[name] = IntDeviceVector.new1(v)
   --v:write_to_file(name..'.index',false)
end
--
--
print("done")
DD = tissueSpec.DD
dt = 0.02
simulateTime = tissueSpec.simulateTime
print("begin stimulate")
endI = math.floor(simulateTime / dt)
oldTime = os.time()
warpSize = 256
timeStimulateBegin = os.time()
for step = 0, endI do
   for _, info in ipairs(allTissues.stimulate) do
	  select_GPU(allTissues[info.voxel].gpu_index)
	  local dst = allTissues[info.voxel].params.istim
	  if step == math.floor(info.startTime/dt) then
		 dst:fillValueAtIndexs(info.stimulate, info.index)
		 print("stimulate voxel ", info.voxel, " at ", info.startTime, 'number', info.index:size())
	  elseif step == math.floor(info.endTime/dt) then
		 dst:fillValueAtIndexs(0, info.index)
	  end
   end
   for _, voxel in ipairs(tissueSpec.allVoxels) do
	  if allTissues[voxel] then
		 local srcGPU = allTissues[voxel].gpu_index
		 select_GPU(srcGPU)
		 local params = allTissues[voxel].params
		 params.updator(warpSize, params)
	  end
   end
   -- this should modify to enumerate all gpu
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
		 MapDoubleDataByIndex(srcItot, params.local2global,  Gcurrent, false)
	  end
   end
   select_GPU(global_gpu_index)
   CalcVolt(Gvolt, Gcurrent, Gsurround, DD, dt, GNewVolt)
   --Gcurrent:write_to_file('allitot', false)
   --GNewVolt:write_to_file('allcell', false)
   for _, voxel in ipairs(tissueSpec.allVoxels) do
	  if allTissues[voxel] then
		 local srcGPU = allTissues[voxel].gpu_index
		 
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
   --warpSize = warpSize + 32
   if step % 50 == 0 then
   --if step then
	  
	  newTime = os.time()
	  print(os.difftime(newTime, oldTime))
	  oldTime = newTime
	  select_GPU(global_gpu_index)
	  MapDoubleDataByIndex(Gvolt, VoxelVoltIndex, VoxelVolt, false)
	  local fname = string.format(tissueSpec.voltFileName, math.floor(step/50) )
	  print('write to', fname)
	  if tissueSpec.writeSeperateFile then
		 VoxelVolt:write_to_file(fname, false)
		 --GNewVolt:write_to_file(fname, false)
	  else
		 VoxelVolt:write_to_file(fname, true)
	  end
   end
end
print('total time', os.difftime(os.time(), timeStimulateBegin))
