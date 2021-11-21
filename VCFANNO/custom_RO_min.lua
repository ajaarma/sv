function overlap_bp(as, ae)
	--print("Inside overlap_bp")
	--print(start)
	--print(stop)
	--print(as)
	--print(ae)
	local start = start+1
	local stop = stop -1
	local diff = math.min(ae,stop) - math.max(as,start)
	diff = diff+1
	--print(result)
	--print("End of overlap_bp\n")
	return string.format("%.4f",diff)
end

function query_overlap(qas, qae)
	--print("Inside query_overlap")
	local qstart = start+1
	local qstop = stop -1
	--
       	diff = math.min(qae,qstop) - math.max(qas,qstart)+1
	if diff	~=0 
	then 
       		query_ov = diff/((qstop-qstart)+1)
	else
		query_ov = "NA"
	end 
	
	--print(result)
	--print(query_ov)
	--print("\n")
	return string.format("%.4f",query_ov)
end

function db_overlap(das, dae)
	--print("Inside db_overlap")
	local qstart = start+1
	local qstop = stop-1
       	diff = math.min(dae,qstop) - math.max(das,qstart)+1
	if diff ~=0
	then 
		db_ov = diff/((dae-das)+1)
	else
		db_ov = "NA"
	end
	--print(result)
	--print(db_ov)
	--print("\n")
	return string.format("%.4f",db_ov)
end

function printOps(ss,se,flag)
	--print("Inside printOps")
	--print(ss)
	--print(se)
	local ss_t = {}
	local se_t = {}

	index = 1
	for ele in string.gmatch(ss,"[^,]+") do
		ss_t[index] = ele
		index = index+1
	end

	index = 1
	for ele in string.gmatch(se,"[^,]+") do
		se_t[index] = ele
		index = index+1
	end
	
	local result = {}
	local ovlap_query = {}
	local ovlap_db = {}

	print("General: ",start,stop,vals,ref,alt)

	for i =1,table.getn(ss_t) do 
		local sst = tonumber(ss_t[i])
		local set = tonumber(se_t[i])
		
		--if flag==0 then
		bp_ovlap = overlap_bp(sst, set)
		--elseif flag==1 then
		q_ovlap = tonumber(query_overlap(sst, set))
		--elseif flag==2 then
		d_ovlap = tonumber(db_overlap(sst, set))
		--end
		--print(ss)
		--print(se)
		--print(sst,set,bp_ovlap,q_ovlap,d_ovlap)
		--print(d_ovlap)

		if q_ovlap ~=nil and d_ovlap ~= nil and q_ovlap >=0.70 and d_ovlap >=0.70 then
			if flag==0 then
				--print("Inside cut-off-0",start,stop)
				--print(ss)
				--print(se)
				--print(sst,set,bp_ovlap,q_ovlap,d_ovlap)
				result[i] = bp_ovlap
			elseif flag==1 then
				--print(q_ovlap)
				result[i] = q_ovlap
			elseif flag==2 then
				--print(d_ovlap)
				result[i] = d_ovlap
			end
		else
			--print("Outside the cut-off: ",start,stop,ss,se,bp_ovlap,q_ovlap,d_ovlap)
		--	print(stop)
		--	print(ss)
		--	print(se)
		--	print(bp_ovlap)
		--	print(q_ovlap)
		--	print(d_ovlap)
		--	print("\n")
			result[i]="NA"
		end
		--overlap_bp[i] = ovlap
	end

	--print(result)
	--print(table.concat(result,","))
	--print("\n")
	--result[1]="NA"
	return table.concat(result,",")

end
