local basis_index = {
	[0b000] = 1, -- 1
	[0b001] = 2, -- x
	[0b010] = 3, -- y
	[0b100] = 4, -- z
	[0b011] = 5, -- xy
	[0b101] = 6, -- xz
	[0b110] = 7, -- yz
	[0b111] = 8 -- i
}
local basis_list = {0b000, 0b001, 0b010, 0b100, 0b011, 0b101, 0b110, 0b111}

local basis_tokens = {
	"one", "x", "y", "z", "xy", "xz", "yz", "im"
}

local EPSLOOKUP = {
	X = {0, 1, 0, 0, 1, 1, 0, 1},
	Y = {0, 0, 1, 0, 1, 0, 1, 1},
	Z = {0, 0, 0, 1, 0, 1, 1, 1}
}

local EPSTBL = {
	[0] = {
		[0] = 1,
		[1] = 1,
	},
	[1] = {
		[0] = 1,
		[1] = -1
	}
}

local function eps3(i, j)
	return EPSTBL[EPSLOOKUP.X[i]][EPSLOOKUP.Z[j]] * EPSTBL[EPSLOOKUP.X[i]][EPSLOOKUP.Y[j]] * EPSTBL[EPSLOOKUP.Y[i]][EPSLOOKUP.Z[j]]
end

local three
local Three
Three = {
	unary = {
		reversion = function(t)
			return three(t[1], t[2], t[3], t[4], -t[5], -t[6], -t[7], -t[8])
		end,
		components = function(t)
			return t[1], t[2], t[3], t[4], t[5], t[6], t[7], t[8]
		end,
		vector = function(t)
			return t[2], t[3], t[4]
		end,
		pseudovector = function(t)
			return t[4], t[5], t[6]
		end,
		pseudoscalar = function(t)
			return t[8]
		end,
		scalar = function(t)
			return t[1]
		end,
		bar = function(t)
			return three(t[1], -t[2], -t[3], -t[4], -t[5], -t[6], -t[7], t[8])
		end,
		sq = function(t)
			return t * t
		end,
		-- dual in the sense that t = t_dual * I
		dual = function(t)
			return three(-t[8], -t[7], t[6], -t[5], -t[4], t[3], -t[2], -t[1])
		end,
		step = function(t)
			if t[8] ~= 0 then
				return 3
			elseif t[7] ~= 0 or t[6] ~= 0 or t[5] ~= 0 then
				return 2
			elseif t[4] ~= 0 or t[3] ~= 0 or t[2] ~= 0 then
				return 1
			end
			return 0
		end,
		sqrt = function(t)
			
		end,
		magnitude = function(t)
			local t_bar = t.bar
			local b = (t * t_bar)
			
			return math.sqrt(math.sqrt(b[1] * b[1] + b[8] * b[8]))
		end,
		unit = function(t)
			-- not 100% certain how to interpret the direction of a multivector with no magnitude
			-- naively discard it
			local mag = t.magnitude
			
			if mag == 0 then
				return Three.lib.zero
			end
			
			return (1 / t.magnitude) * t
		end,
		inverse = function(t)
			-- inverse function expanded for higher performance. not verified closely but experimentally accurate
			-- defined in accordance to the resource specified via method #inv2
			-- floating point error lives here
			
			local b_0 = t[1] * t[1] - t[2] * t[2] - t[3] * t[3] - t[4] * t[4] + t[5] * t[5] + t[6] * t[6] + t[7] * t[7] - t[8] * t[8]
			local b_7 = 2 * (-t[4] * t[5] + t[3] * t[6] - t[2] * t[7] + t[1] * t[8])
			
			local s = b_0 * b_0 + b_7 * b_7
			
			if s == 0 then
				return Three.lib.zero
			end
			
			return (1 / s) * three(b_0 * t[1] + b_7 * t[8],
				-b_0 * t[2] - b_7 * t[7],
				-b_0 * t[3] + b_7 * t[6],
				-b_0 * t[4] - b_7 * t[5],
				-b_0 * t[5] + b_7 * t[4],
				-b_0 * t[6] - b_7 * t[3],
				-b_0 * t[7] + b_7 * t[2],
				b_0 * t[8] - b_7 * t[1])
		end,
		inv2 = function(t)
			-- defined verbatim from https://arxiv.org/pdf/1712.05204.pdf (p. 6)
			local t_bar = t.bar
			local b = t * t_bar
			
			local s = (b[1] * b[1] + b[8] * b[8])
			
			if s == 0 then
				return Three.lib.zero
			end

			return (1 / (b[1] * b[1] + b[8] * b[8])) * t_bar * three(b[1], 0, 0, 0, 0, 0, 0, -b[8])
		end,
	},
	lib = {
		v3 = function(arg1, arg2, arg3)
			if typeof(arg1) == "Vector3" then
				arg1, arg2, arg3 = arg1.X, arg1.Y, arg1.Z
			end
			return three(0, arg1, arg2, arg3, 0, 0, 0, 0)
		end,
		inner = function(t1, t2)
			return 0.5 * (t1 * t2 + t2 * t1)
		end,
		outer = function(t1, t2)
			return 0.5 * (t1 * t2 - t2 * t1)
		end,
		regressive = function(t1, t2)
			return Three.lib.inner(t1.dual, t2)
		end,
		project = function(t1, t2)
			return (t1 .. t2) * t2.inverse
		end,
		reject = function(t1, t2)
			return (t1 ^ t2) * t2.inverse
		end,
		reflect = function(t1, t2)
			return -t2 * t1 * t2.inverse
		end,
		lerp = function(t1, t2, alpha, plane)
			
		end,
		rotateTo = function(t1, t2, theta, plane)
			if not plane then
				plane = t1 ^ t2
			end
			
			
		end,
		axisAngle = function(t1, axis, theta)
			local plane = axis.dual
			
			
		end,
		rotateOn = function(plane, t1, t2, theta)
			
		end,
		-- this is the exponential function for bivectors (plane segments) only (assumed t1.scalar = t1.vector = t1.pseudoscalar = 0)
		biExp = function(t1, theta)
			return math.cos(theta) + math.sin(theta) * t1
		end,
	},
	__unm = function(t)
		return three(-t[1], -t[2], -t[3], -t[4], -t[5], -t[6], -t[7], -t[8])
	end,
	__add = function(t1, t2)
		if typeof(t1) == "number" then
			return three(t1 + t2[1], t2[2], t2[3], t2[4], t2[5], t2[6], t2[7], t2[8])
		elseif typeof(t2) == "number" then
			return three(t2 + t1[1], t1[2], t1[3], t1[4], t1[5], t1[6], t1[7], t1[8])
		end
		return three(t1[1] + t2[1], t1[2] + t2[2], t1[3] + t2[3], t1[4] + t2[4], t1[5] + t2[5], t1[6] + t2[6], t1[7] + t2[7], t1[8] + t2[8])
	end,
	__sub = function(t1, t2)
		if typeof(t1) == "number" then
			return three(t1 - t2[1], t2[2], t2[3], t2[4], t2[5], t2[6], t2[7], t2[8])
		elseif typeof(t2) == "number" then
			return three(t1[1] - t2, t1[2], t1[3], t1[4], t1[5], t1[6], t1[7], t1[8])
		end
		return three(t1[1] - t2[1], t1[2] - t2[2], t1[3] - t2[3], t1[4] - t2[4], t1[5] - t2[5], t1[6] - t2[6], t1[7] - t2[7], t1[8] - t2[8])
	end,
	__mul = function(t1, t2)
		if typeof(t1) == "number" then
			return three(t1 * t2[1], t1 * t2[2], t1 * t2[3], t1 * t2[4], t1 * t2[5], t1 * t2[6], t1 * t2[7], t1 * t2[8])
		elseif typeof(t2) == "number" then
			return three(t2 * t1[1], t2 * t1[2], t2 * t1[3], t2 * t1[4], t2 * t1[5], t2 * t1[6], t2 * t1[7], t2 * t1[8])
		end
		if typeof(t2) == "Vector3" then
			return t1 * three(0, t2.X, t2.Y, t2.Z)
		elseif typeof(t2) == "Vector2" then
			return t1 * three(0, t2.X, t2.Y)
		end
		local t_prod = {0, 0, 0, 0, 0, 0, 0, 0}
		
		for i = 1, 8 do
			if t1[i] == 0 then
				continue
			end
			for j = 1, 8 do
				if t2[j] == 0 then
					continue
				end
--				print("tick")
				t_prod[basis_index[bit32.bxor(basis_list[i], basis_list[j])]] += EPSTBL[EPSLOOKUP.X[j]][EPSLOOKUP.Z[i]] * EPSTBL[EPSLOOKUP.X[j]][EPSLOOKUP.Y[i]] * EPSTBL[EPSLOOKUP.Y[j]][EPSLOOKUP.Z[i]] * t1[i] * t2[j]
			end
		end
		
		return three(t_prod[1], t_prod[2], t_prod[3], t_prod[4], t_prod[5], t_prod[6], t_prod[7], t_prod[8])
	end,
	__div = function(t1, t2)
		if typeof(t1) == "number" then
			return t1 * t2.Inverse
		elseif typeof(t2) == "number" then
			return t1 * (1 / t2)
		end
		if typeof(t2) == "Vector3" then
			--t2 = three(0, t2.X, t2.Y, t2.Z).Inverse
			local sqmag = t2.X * t2.X + t2.Y * t2.Y + t2.Z * t2.Z
			t2 = three(0, t2.X / sqmag, t2.Y / sqmag, t2.Z / sqmag)
			return t1 * t2
		elseif typeof(t2) == "Vector2" then
			local sqmag = t2.X * t2.X + t2.Y * t2.Y
			t2 = three(0, t2.X / sqmag, t2.Y / sqmag)
			return t1 * t2
		end
		
		return t1 * t2.Inverse
	end,
	__pow = function(t1, t2)
		return Three.lib.outer(t1, t2)
	end,
	__concat = function(t1, t2)
		return Three.lib.inner(t1, t2)
	end,
	__mod = function(t1, t2)
		return Three.lib.regressive(t1, t2)
	end,
	__eq = function(t1, t2)
		if typeof(t1) == "number" then
			for i = 2, 8 do
				if t2[i] ~= 0 then
					return false
				end
			end
			return t2[1] == t1
		elseif typeof(t1) == "Vector3" then
			if t2[1] ~= 0 then
				return false
			end
			for i = 5, 8 do
				if t2[i] ~= 0 then
					return false
				end
			end
			return t2[2] == t1.X and t2[3] == t1.Y and t2[4] == t1.Z
		end
		if typeof(t2) == "number" then
			for i = 2, 8 do
				if t1[i] ~= 0 then
					return false
				end
			end
			return t1[1] == t2
		elseif typeof(t2) == "Vector3" then
			if t1[1] ~= 0 then
				return false
			end
			for i = 5, 8 do
				if t1[i] ~= 0 then
					return false
				end
			end
			return t1[2] == t2.X and t1[3] == t2.Y and t1[4] == t2.Z
		end
		for i = 1, 8 do
			if t1[i] ~= t2[i] then
				return false
			end
		end
		return true
	end,
	__tostring = function(t)
		local strt = {}
		for i = 1, 8 do
			if t[i] ~= 0 then
				strt[#strt + 1] = string.format("%s%.3f%s", t[i] < 0 and "－" or (#strt > 0 and "＋" or ""), math.abs(t[i]), basis_tokens[i])
			end
		end
		if #strt == 0 then
			return "𝟎"
		end
		return table.concat(strt)
	end,
}
Three.unary.Magnitude = Three.unary.magnitude
Three.unary.Unit = Three.unary.unit
Three.unary.Inverse = Three.unary.inverse
Three.unary.Reversion = Three.unary.reversion
Three.unary.Bar = Three.unary.bar
Three.__index = function(t, k)
	if Three.unary[k] then
		return Three.unary[k](t)
	end
	if Three.lib[k] then
		return Three.lib[k]
	end
	return Three[k]
end
three = function(a, x, y, z, xy, xz, yz, xyz)
	if typeof(a) == "Vector3" then
		x = a.X
		y = a.Y
		z = a.Z
		a = 0
	end
	local v = setmetatable({a or 0, x or 0, y or 0, z or 0, xy or 0, xz or 0, yz or 0, xyz or 0}, Three)
	
	return v
end

Three.lib.im = three(0, 0, 0, 0, 0, 0, 0, 1)
Three.lib.iminv = three(0, 0, 0, 0, 0, 0, 0, -1)
Three.lib.zero = three()

return {
	three = setmetatable({}, {
		__index = function(t, k)
			return Three.lib[k]
		end,
		__call = function(_, ...)
			return three(...)
		end,
	})
}