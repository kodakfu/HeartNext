function mostly_equal_sum_indexs(data, sum)
   local dp = {}
   local pre = {}
   for i, v in ipairs(data) do
	  table.insert(dp, 0)
	  table.insert(pre, i)
   end
   -- construct dp and pre
   for i = 1, table.maxn(data) do
	  dp[i] = data[i]
	  pre[i] = i
	  for j = 1, i - 1 do
		 local v = dp[j] + data[i]
		 if math.abs(v - sum) < math.abs(dp[i] - sum) then
			dp[i] = v
			pre[i] = j
		 end
	  end
   end
   local startIndex = 1
   local nearest_dist = math.abs(dp[startIndex] - sum)
   for i = 2, table.maxn(dp) do
	  local dist_i = math.abs(dp[i] - sum)
	  if nearest_dist >  dist_i then
		 nearest_dist = dist_i
		 startIndex = i
	  end
   end
   local preIndex = pre[startIndex]
   local rs = {startIndex}
   while startIndex ~= preIndex do
	  table.insert(rs, preIndex)
	  startIndex = preIndex
	  preIndex = pre[startIndex]
   end
   return rs
end

function select_GPU_for_cell(allTissues, allVoxels)
   -- WARNING! can only handle 2 GPU
   local gpu_number = 2
   local avg_cell_per_gpu = math.floor(allTissues.totalCellNumber / gpu_number)
   local sum = 0
   local gpu_index = 0
   local tissue_size_array = {}
   local tissue_array = {}
   for i, v in ipairs(allVoxels) do
	  table.insert(tissue_size_array, allTissues[v]:size())
	  table.insert(tissue_array, allTissues[v])
   end
   local part1 = mostly_equal_sum_indexs(tissue_size_array, avg_cell_per_gpu)
   for i,v in ipairs(part1) do
	  local tissue = table.remove(tissue_array, v)
	  print('make',allVoxels[v], 'GPU0')
	  tissue.gpu_index = 0
   end
   for i,v in ipairs(tissue_array) do
	  print('make this GPU1')
	  v.gpu_index = 1
   end
end
function fillCellUse06VentricleModel(tissues, voxels, CellVoxel)
   local name2Type = {endo=1, mcell=2, epi=3}
   local VoxelToCellType = {}
   for k, vs in pairs(CellVoxel) do
	  for _, v in ipairs(vs) do
		 VoxelToCellType[v] = k
	  end
   end
   local init_v = -86.2
   for i, v in ipairs(voxels) do
	  local theTissue = tissues[v]
	  if theTissue then
		 select_GPU(theTissue.gpu_index)
		 local cellNumber = theTissue:size()
		 theTissue.params = {}
		 if global_gpu_index ~= theTissue.gpu_index then
			select_GPU(global_gpu_index)
			theTissue.params.itot_in_global_gpu = DoubleDeviceVector.new2(cellNumber, 0)
			theTissue.params.volt_in_global_gpu = DoubleDeviceVector.new2(cellNumber, init_v)
		 end
		 theTissue.params.local2global = IntDeviceVector.new1(theTissue)
		 
		 select_GPU(theTissue.gpu_index)
		 theTissue.params.cellNumber = cellNumber
		 local cellType = name2Type[VoxelToCellType[v]]
		 if not cellType then error('nill cell type') end
		 theTissue.params.type = IntDeviceVector.new2(cellNumber, cellType)
		 theTissue.params.volt = DoubleDeviceVector.new2(cellNumber, init_v)
		 theTissue.params.volt.needSerialize = true
		 theTissue.params.itot = DoubleDeviceVector.new2(cellNumber, 0)
		 theTissue.params.itot.needSerialize = true
		 theTissue.params.istim = DoubleDeviceVector.new2(cellNumber, 0)
		 theTissue.params.istim.needSerialize = true
		 theTissue.params.updator = UpdateVentricleCurrent_06Model
		 for name, init in pairs(return06ParamSpec()) do
			if name == 'CaSR' and VoxelToCellType[v] == 'mcell' then
			   print("init = init * 3 for voxel", v)
			   init = init * 3
			end
			theTissue.params[name] = DoubleDeviceVector.new2(cellNumber, init)
			theTissue.params[name].needSerialize = true
		 end
	  end
   end
end
function fillCellUseBlockKrCalVentricleModel(tissues, voxels, CellVoxel)
   local name2Type = {endo=1, mcell=2, epi=3}
   local VoxelToCellType = {}
   for k, vs in pairs(CellVoxel) do
	  for _, v in ipairs(vs) do
		 VoxelToCellType[v] = k
	  end
   end
   local init_v = -86.2
   for i, v in ipairs(voxels) do
	  local theTissue = tissues[v]
	  if theTissue then
		 select_GPU(theTissue.gpu_index)
		 local cellNumber = theTissue:size()
		 theTissue.params = {}
		 if global_gpu_index ~= theTissue.gpu_index then
			select_GPU(global_gpu_index)
			theTissue.params.itot_in_global_gpu = DoubleDeviceVector.new2(cellNumber, 0)
			theTissue.params.volt_in_global_gpu = DoubleDeviceVector.new2(cellNumber, init_v)
		 end
		 theTissue.params.local2global = IntDeviceVector.new1(theTissue)
		 
		 select_GPU(theTissue.gpu_index)
		 theTissue.params.cellNumber = cellNumber
		 local cellType = name2Type[VoxelToCellType[v]]
		 if not cellType then error('nill cell type') end
		 theTissue.params.block_kr = 0.43
		 theTissue.params.block_cal = 0.85
		 theTissue.params.type = IntDeviceVector.new2(cellNumber, cellType)
		 theTissue.params.volt = DoubleDeviceVector.new2(cellNumber, init_v)
		 theTissue.params.volt.needSerialize = true
		 theTissue.params.itot = DoubleDeviceVector.new2(cellNumber, 0)
		 theTissue.params.itot.needSerialize = true
		 theTissue.params.istim = DoubleDeviceVector.new2(cellNumber, 0)
		 theTissue.params.istim.needSerialize = true
		 theTissue.params.updator = UpdateVentricleCurrent_BlockKrCal 
		 for name, init in pairs(returnBlockKrCalParamSpec()) do
			theTissue.params[name] = DoubleDeviceVector.new2(cellNumber, init)
			theTissue.params[name].needSerialize = true
		 end
	  end
   end
end


function fillCellUsePurkinjeModel(tissues, voxels)
   local init_v = -74.78905
   for i,v in ipairs(voxels) do
	  local tissue = tissues[v]
	  if tissue then
		 select_GPU(tissue.gpu_index)
		 local number = tissue:size()
		 tissue.params={}
		 if global_gpu_index ~= tissue.gpu_index then
			select_GPU(global_gpu_index)
			tissue.params.itot_in_global_gpu = DoubleDeviceVector.new2(number, 0)
			tissue.params.volt_in_global_gpu = DoubleDeviceVector.new2(number, init_v)
		 end
		 tissue.params.local2global = IntDeviceVector.new1(tissue)
		 select_GPU(tissue.gpu_index)
		 
		 tissue.params.cellNumber = number
		 tissue.params.volt = DoubleDeviceVector.new2(number, init_v)
		 tissue.params.volt.needSerialize = true
		 tissue.params.itot = DoubleDeviceVector.new2(number, 0)
		 tissue.params.itot.needSerialize = true
		 tissue.params.istim = DoubleDeviceVector.new2(number, 0)
		 tissue.params.istim.needSerialize = true
		 tissue.params.updator = UpdatePurkinjeCurrent
		 for name, init in pairs(returnPurkinjeParamSpec()) do
			tissue.params[name] = DoubleDeviceVector.new2(number, init)
			tissue.params[name].needSerialize = true
		 end
	  end
   end
end

function filterVoxelByPosition(localIndex, indexInTissueFile, judge, dimx, dimy)
   local trs = {}
   local last = localIndex:size() - 1
   for i = 0, last do
	  local linearIndex = indexInTissueFile:at(localIndex:at(i))
	  local z = math.floor(linearIndex / (dimx * dimy) )
	  local mid = linearIndex - z * (dimx * dimy)
	  local y = math.floor(mid / dimx)
	  local x = mid % dimx
	  if judge(x,y,z) then
		 table.insert(trs, i)
	  end
   end
   local rs = IntHostVector.new2(table.maxn(trs), 0)
   for k, v in ipairs(trs) do rs:set(k-1, v) end
   return rs
end
function makeStimulate(tissues, voxel, strength, s, e, canIt)
   local rs = {
	  voxel = voxel,
	  stimulate = strength,
	  startTime = s, endTime = e,
	  index = filterVoxelByPosition(tissues[voxel], tissues.indexInTissueFile, canIt, tissues.dimx, tissues.dimy)
   }
   return rs
end

