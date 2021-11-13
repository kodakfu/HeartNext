dofile('core.lua')--引入该文件中函数，包含makeStimulate
--无浦肯野纤维细胞
--读取切片数据，设置切片数据的详细信息
local raw = IntHostVector.new0()--用来存储组织体数据
tissueFileName = "fang_ischemia1a_1b_infarct.dat"--体数据文件名
LoadIntFromByteFile(tissueFileName, raw)--读取体数据
dimx = 600
dimy = 600
dimz = 1
nilVoxel = 0--无细胞的体素值
ignoreTable = {
     [20]='',[40]='',[64]=''
}--每个细胞与哪些细胞不传导,[20]='21 22'
cellVoxel = {64,20,40}--所有细胞的体素
voxel2cellType={[64]=6,[20]=9,[40]=12}--细胞体素到细胞的类型的映射
--解析组织数据
parsedTissue = parse_tissue(dimx, dimy, dimz, raw, nilVoxel, ignoreTable)
--设置初值
for i, cell in ipairs(cellVoxel) do
   local tissue = parsedTissue[cell]
   local ncell = tissue:size()
   local init_v = -87.5--bai_change -86.2
   local cellType = voxel2cellType[cell]
   --模型通用的参数
   tissue.params = {
	  local2global = IntDeviceVector.new1(tissue),
	  cellNumber = ncell,
	  type = IntDeviceVector.new2(ncell, cellType),
	  volt = DoubleDeviceVector.new2(ncell, init_v),
	  itot = DoubleDeviceVector.new2(ncell, 0),
	  istim = DoubleDeviceVector.new2(ncell, 0),
	  updator = UpdateVentricleCurrent_06model_ischemia_1a_1b_infarct,
	  resetor = function(theTissue, needResetIndex)
		    theTissue.params.volt:fillValueAtIndexs(init_v, needResetIndex)
		    for name, init in pairs(return06ParamSpec_ischemia_1a_1b_infarct()) do
		       theTissue.params[name]:fillValueAtIndexs(init, needResetIndex)
		    end
	end,
   }
   local params = tissue.params
   params.volt.needSerialize = true
   params.itot.needSerialize = true
   params.istim.needSerialize = true
   --特定模型参数
   for name, init in pairs(return06ParamSpec_ischemia_1a_1b_infarct()) do
	  tissue.params[name] = DoubleDeviceVector.new2(ncell, init)
	  tissue.params[name].needSerialize = true
   end
end
--设置刺激
parsedTissue.stimulate = {}
function pointInBalls(balls, x, y, z)
   for _, ball in ipairs(balls) do
	  local local_x = x - ball[1]
	  local local_y = y - ball[2]
	  local local_z = z - ball[3]
	  local r = ball[4]
	  if local_x * local_x + local_y * local_y + local_z * local_z < r * r then return true end
   end
   return false
end

k=1000
for pixel = 1, 30, 1 do
	s1_name=string.format("s11_%s",tostring(pixel))
    parsedTissue.stimulate[s1_name] = makeStimulate(parsedTissue, 20, -120, (pixel-1)*k,  (pixel-1)*k+1,
            function(x, y, z)
                if x<=3  then
                    return true 
                else
                    return false
                end
            end
    )
s1_name=string.format("s12_%s",tostring(pixel))
    parsedTissue.stimulate[s1_name] = makeStimulate(parsedTissue, 40, -120, (pixel-1)*k,  (pixel-1)*k+1,
            function(x, y, z)
                if x<=3  then
                    return true 
                else
                    return false
                end
            end
    )
s1_name=string.format("s13_%s",tostring(pixel))
    parsedTissue.stimulate[s1_name] = makeStimulate(parsedTissue, 64, -120, (pixel-1)*k,  (pixel-1)*k+1,
            function(x, y, z)
                if x<=3  then
                    return true 
                else
                    return false
                end
            end
    )
end

k=250
for pixel = 1, 40, 1 do
	s1_name=string.format("s21_%s",tostring(pixel))
    parsedTissue.stimulate[s1_name] = makeStimulate(parsedTissue, 20, -120, (pixel-1)*k+30000,  (pixel-1)*k+30000+1,
            function(x, y, z)
                if x<=3 then
                    return true 
                else
                    return false
                end
            end
    )
s1_name=string.format("s22_%s",tostring(pixel))
    parsedTissue.stimulate[s1_name] = makeStimulate(parsedTissue, 40, -120, (pixel-1)*k+30000,  (pixel-1)*k+30000+1,
            function(x, y, z)
                if x<=3 then
                    return true 
                else
                    return false
                end
            end
    )
s1_name=string.format("s23_%s",tostring(pixel))
    parsedTissue.stimulate[s1_name] = makeStimulate(parsedTissue, 64, -120, (pixel-1)*k+30000,  (pixel-1)*k+30000+1,
            function(x, y, z)
                if x<=3 then
                    return true 
                else
                    return false
                end
            end
    )
end



function resetCell(tissue, cell, callback)
   local needReset = filterVoxelByPosition(tissue[cell], tissue.indexInTissueFile, callback, dimx, dimy)
   tissue[cell].params.resetor(tissue[cell], needReset)
end
--分配global相关的数据
Gvolt = DoubleDeviceVector.new2(parsedTissue.totalCellNumber, 0)
GNewVolt = DoubleDeviceVector.new2(parsedTissue.totalCellNumber, 0)
Gcurrent = DoubleDeviceVector.new2(parsedTissue.totalCellNumber, 0)

VoxelVolt = DoubleDeviceVector.new2(dimx * dimy * dimz, -100)--全体体素（包括空体素）
VoxelVoltIndex = IntDeviceVector.new1(parsedTissue.indexInTissueFile)--细胞体素到全体体素映射

Gsurround = {}--存储在GPU上的surround
for name,v in pairs(parsedTissue.surround) do
   Gsurround[name] = IntDeviceVector.new1(v)
end
--开始仿真迭代
DD = 0.154/(0.25*0.25)--扩散系数
simulateTime = 41000
dt = 0.02
endStep = math.floor(simulateTime/dt)
oldTime = os.time()
saveAtStep = math.floor(500/dt)
toLoad = false
startStep = 0
--加载状态
if toLoad then
   print('load------------------------')
   local dir2load = '/content/graphic/data2/snap_ischemia1a_1b_infarct_fangxing/'
   loadData(dir2load, parsedTissue, cellVoxel)
   startStep = savedStep
end
--为Gvolt设置初值
function localVolt2Global()
   for i, cell in ipairs(cellVoxel) do
      local params = parsedTissue[cell].params--这种细胞的参数
      MapDoubleDataByIndex(params.volt, params.local2global, Gvolt, false)--将local电压映射到global电压
   end

end
localVolt2Global()

for step = startStep, endStep do
   --保存状态
   if saveAtStep == step then
	  print('save---------------------')
	  --saveData('/content/graphic/data2/snap_ischemia1a_1b_infarct_fangxing/', step, parsedTissue, cellVoxel)--liang_change
   end
   --加刺激
   for name, info in pairs(parsedTissue.stimulate) do
	  local dst = parsedTissue[info.voxel].params.istim
	  if step == math.floor(info.startTime/dt) then
		 print("apply stimulate",name)
		 dst:fillValueAtIndexs(info.stimulate, info.index)--在index位置赋值
	  elseif step == math.floor(info.endTime/dt) then
		 print("end stimulate",name)
		 dst:fillValueAtIndexs(0, info.index)
	  end
   end
   --切波
   if 43100000 == step * dt then
      
      function decider(x,y,z)
        if y <300 then
	-- if y <300  and x > 300 then
	    return true
	 else
	    return false
	 end
      end
      local cellNeedReset = 20
      resetCell(parsedTissue, cellNeedReset, decider)
      localVolt2Global()
   end
   --更新电流
   for i, voxel in ipairs(cellVoxel) do
	  local params = parsedTissue[voxel].params
	  params.updator(256, params)--调用更新电流函数
	  MapDoubleDataByIndex(params.itot, params.local2global, Gcurrent, false)--映射到全局电流
   end
   CalcVolt(Gvolt, Gcurrent, Gsurround, DD, dt, GNewVolt)--计算新电压
   --将电压映射回local
   for i, voxel in ipairs(cellVoxel) do
	  local params = parsedTissue[voxel].params
	  MapDoubleDataByIndex(params.volt, params.local2global, GNewVolt, true)--true是反向
   end
   
   Gvolt, GNewVolt = GNewVolt, Gvolt--更新旧电压
   --写数据
   if step % 500 == 0 and step>=30000*50   then --step>900*4*50 
	  newTime = os.time()
	  print('cost time',os.difftime(newTime, oldTime))--距离上次写数据花了多少时间
	  oldTime = newTime
	  
	  MapDoubleDataByIndex(Gvolt, VoxelVoltIndex, VoxelVolt, false)--将细胞电压映射到全体体素电压
	  fname = string.format("/content/graphic/data2/volts_ischemia1a_1b_infarct_fangxing/%d.volt", math.floor(step/50))--构造文件名  
	  print("write to", fname)
	  VoxelVolt:write_to_file(fname, false)--覆盖(false)写文件
   end
   
end
