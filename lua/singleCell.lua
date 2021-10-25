cellType = 2
type2name = {'endo','mcell','epi'}
for cellType = 1,3 do
   print( 'for', type2name[cellType])
   cell = {
	  
   }
   cell.params = {
	  volt = DoubleDeviceVector.new2(1, -86.2),
	  itot = DoubleDeviceVector.new2(1, 0),
	  istim = DoubleDeviceVector.new2(1, 0),
	  type = IntDeviceVector.new2(1, cellType),
	  cellNumber = 1,
   }
   for name, init in pairs(return06ParamSpec()) do
	  cell.params[name]= DoubleDeviceVector.new2(1, init)
   end

   dt = 0.02
   simulateTime = 800 --ms
   endI = math.floor(simulateTime/dt)
   warpSize = 32
   time = 0
   for step = 0, endI do
	  if time > 0 and time < 1 then
		 cell.params.istim:set(0, -52)
	  else
		 cell.params.istim:set(0, 0)
	  end
	  UpdateVentricleCurrent_06Model(32, cell.params)
	  local volt = cell.params.volt:at(0)
	  volt = volt + dt * (-cell.params.itot:at(0))
	  cell.params.volt:set(0,volt)
	  if step % 50 == 0 then
		 cell.params.volt:write_to_file(type2name[cellType], true)
	  end
	  time = time + dt
   end


end

