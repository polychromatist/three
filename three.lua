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

local EPSLOOKUP_X = {0, 1, 0, 0, 1, 1, 0, 1}
local EPSLOOKUP_Y = {0, 0, 1, 0, 1, 0, 1, 1}
local EPSLOOKUP_Z = {0, 0, 0, 1, 0, 1, 1, 1}

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

local PHI = {
	nil,
	{4, 6, 5, 3},
	{2, 5, 7, 4},
	{3, 7, 6, 2},
	{2, 3, 6, 7},
	{4, 2, 7, 5},
	{3, 4, 5, 6}
}

local function eps3(i, j)
	return EPSTBL[EPSLOOKUP_X[i]][EPSLOOKUP_Z[j]] * EPSTBL[EPSLOOKUP_X[i]][EPSLOOKUP_Y[j]] * EPSTBL[EPSLOOKUP_Y[i]][EPSLOOKUP_Z[j]]
end

local bxor = bit32.bxor

local three
local Three
Three = {
	unary = {
		reversion = function(t)
			return three(t[1], t[2], t[3], t[4], -t[5], -t[6], -t[7], -t[8])
		end,
		gradeinvol = function(t)
			return three(t[1], -t[2], -t[3], -t[4], t[5], t[6], t[7], -t[8])
		end,
		-- this value is not defined for general multivectors. this is defined for k-blades (products of linearly independent vectors)
		grade = function(t)
			if t[8] > 1e-6 or t[8] < -1e-6 then
				return 3
			elseif t[7] > 1e-6 or t[7] < -1e-6 or t[6] > 1e-6 or t[6] < -1e-6 or t[5] > 1e-6 or t[5] < -1e-6 then
				return 2
			elseif t[4] > 1e-6 or t[4] < -1e-6 or t[3] > 1e-6 or t[3] < -1e-6 or t[2] > 1e-6 or t[2] < -1e-6 then
				return 1
			end
			return 0
		end,
		components = function(t)
			return t[1], t[2], t[3], t[4], t[5], t[6], t[7], t[8]
		end,
		vector = function(t)
			return Vector3.new(t[2], t[3], t[4])
		end,
		pseudovector = function(t)
			return three(0, 0, 0, 0, t[5], t[6], t[7], 0)
		end,
		pseudoscalar = function(t)
			return three(0, 0, 0, 0, 0, 0, 0, t[8])
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
			return three(t[8], t[7], -t[6], t[5], -t[4], t[3], -t[2], -t[1])
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
		-- can return 0 for multivectors composed of a null versor
		-- another formulation: math.sqrt(t.reversion:scalarprod(t)) per Leo Dorst
		-- Applications of Geometric Algebra in Computer Science and Engineering, Ch. 2 The Inner Products of Geometric Algebra, p. 40
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
		vec = function(arg1, arg2, arg3, arg4)
			if typeof(arg1) == "Vector3" then
				arg4 = arg2
				arg1, arg2, arg3 = arg1.X, arg1.Y, arg1.Z
			end
			return three(0, arg1, arg2, arg3, 0, 0, 0, arg4 or 0)
		end,
		bivec = function(arg1, arg2, arg3, arg4)
			if typeof(arg1) == "Vector3" then
				arg4 = arg2
				arg1, arg2, arg3 = arg1.X, arg1.Y, arg1.Z
			end
			return three(arg4 or 0, 0, 0, 0, arg1, arg2, arg3, 0)
		end,
		-- this product also has the identity:
		-- 2 * a:inner(b) = (a + b).sq - a.sq - b.sq
		inner = function(t1, t2)
			return 0.5 * (t1 * t2 + t2 * t1)
		end,
		outer = function(t1, t2)
			return 0.5 * (t1 * t2 - t2 * t1)
		end,
		-- direct implementation of exterior product
		exterior = function(t1, t2)
			local prod = three()
			
			local t1r, t2s, term_rs
			for r = 0, 3 do
				t1r = three(t1:gradeproject(r))
				for s = 0, 3 do
					t2s = three(t2:gradeproject(s))
					
					term_rs = (t1r * t2s):gradeproject(r + s)
					
					if term_rs then
						prod += term_rs
					end
					
					term_rs = nil
				end
			end
			
			return prod
			
			--[=[
			-- not sure why i can take the highest grade of multivectors and just use that but ok
			
			return (t1 * t2):gradeproject(t1.grade + t2.grade) or Three.lib.zero]=]
		end,
		-- the scalar part of the geometric product
		-- note that sometimes "the scalar product" refers to the regular inner product
		-- also note the following identity:
		-- t1:scalarprod(t2) + t1:fatdot(t2) = t1:left(t2) + t1:right(t2)
		scalarprod = function(t1, t2)
			-- equivalent to (t1 * t2):gradeproject(0)
			local prod = 0
			
			for i = 1, 4 do
				prod += t1[i] * t2[i] - t1[i + 4] * t2[i + 4]
			end
			
			return prod
		end,
		fatdot = function(t1, t2)
			local prod = three()

			local t1r, t2s, term_rs
			for r = 0, 3 do
				t1r = three(t1:gradeproject(r))
				for s = 0, 3 do
					t2s = three(t2:gradeproject(s))

					term_rs = (t1r * t2s):gradeproject(math.abs(s - r))

					if term_rs then
						prod += term_rs
					end

					term_rs = nil
				end
			end

			return prod
		end,
		left = function(t1, t2)
			local prod = three()

			local t1r, t2s, term_rs
			for r = 0, 3 do
				t1r = three(t1:gradeproject(r))
				for s = 0, 3 do
					t2s = three(t2:gradeproject(s))

					term_rs = (t1r * t2s):gradeproject(s - r)

					if term_rs then
						prod += term_rs
					end

					term_rs = nil
				end
			end

			return prod
		end,
		right = function(t1, t2)
			local prod = three()

			local t1r, t2s, term_rs
			for r = 0, 3 do
				t1r = three(t1:gradeproject(r))
				for s = 0, 3 do
					t2s = three(t2:gradeproject(s))

					term_rs = (t1r * t2s):gradeproject(r - s)

					if term_rs then
						prod += term_rs
					end

					term_rs = nil
				end
			end

			return prod
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
		-- get a multivector's k-th grade component
		gradeproject = function(t1, k)
			if k == 0 then
				return t1.scalar
			elseif k == 1 then
				return t1.vector
			elseif k == 2 then
				return t1.pseudovector
			elseif k == 3 then
				return t1.pseudoscalar
			end
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
		-- for multivectors with many nonzero components: more scalable implementation of geometric product
		-- 100000 repeated multiplications of the multivectors three(1, 2, 3, 4, 5, 6, 7, 8) * three(8, 7, 6, 5, 4, 3, 2, 1):
		-- conventional implementation: 1.35s runtime
		-- geo2: 0.22s runtime
		-- less or equally efficient in general for sparse multivectors, like pure vectors, pure pseudovectors or versors (Aone + Bim)
		geo2 = function(a, b)
			local c = {0, 0, 0, 0, 0, 0, 0, 0}
			
			a[6] = -a[6]
			b[6] = -b[6]
			
			for i = 1, 4 do
				c[1] += a[i] * b[i] - a[i + 4] * b[i + 4]
			end
			do
				local phi_i1
				for i = 2, 4 do
					phi_i1 = PHI[i]
					c[i] = a[1] * b[i] + a[i] * b[1] + a[phi_i1[1]] * b[phi_i1[2]] - a[phi_i1[2]] * b[phi_i1[1]] + a[phi_i1[3]] * b[phi_i1[4]] - a[phi_i1[4]] * b[phi_i1[3]] + a[8] * b[9 - i] + a[9 - i] * b[8] 
				end
				for i = 5, 7 do
					phi_i1 = PHI[i]
					-- note: if a == b, the phi terms are cancelling
					-- therefore we are left with the scalar and the pseudoscalar contribution only
					-- could be interesting in a sqrt implementation
					if i ~= 6 then
						c[i] = a[1] * b[i] + a[i] * b[1] + a[phi_i1[1]] * b[phi_i1[2]] - a[phi_i1[2]] * b[phi_i1[1]] + a[phi_i1[3]] * b[phi_i1[4]] - a[phi_i1[4]] * b[phi_i1[3]] - a[8] * b[9 - i] + a[9 - i] * b[8]
					else
						c[i] = -(a[1] * b[i] + a[i] * b[1] + a[phi_i1[1]] * b[phi_i1[2]] - a[phi_i1[2]] * b[phi_i1[1]] + a[phi_i1[3]] * b[phi_i1[4]] - a[phi_i1[4]] * b[phi_i1[3]] - a[8] * b[9 - i] + a[9 - i] * b[8])
					end 
				end
			end
			-- without lookup table
			--[=[
			do
				local pi1x, pi1y, pi2x, pi2y
				for i = 2, 4 do
					pi1x, pi1y = (i % 3) + 2, ((1 - i) % 3) + 2
					pi2x, pi2y = 9 - pi1x, 9 - pi1y
					c[i] = a[1] * b[i] + a[i] * b[1] + a[pi1x] * b[pi1y] - a[pi1y] * b[pi1x] + a[pi2x] * b[pi2y] - a[pi2y] * b[pi2x] + a[8] * b[9 - i] + a[9 - i] * b[8] 
				end
				for i = 5, 7 do
					pi1x, pi1y = ((5 - i) % 3) + 2, ((6 - i) % 3) + 2
					pi2x, pi2y = 9 - pi1y, 9 - pi1x
					if i ~= 6 then
						c[i] = a[1] * b[i] + a[i] * b[1] + a[pi1x] * b[pi1y] - a[pi1y] * b[pi1x] + a[pi2x] * b[pi2y] - a[pi2y] * b[pi2x] - a[8] * b[9 - i] + a[9 - i] * b[8]
					else
						
						c[i] = -(a[1] * b[i] + a[i] * b[1] + a[pi1x] * b[pi1y] - a[pi1y] * b[pi1x] + a[pi2x] * b[pi2y] - a[pi2y] * b[pi2x] - a[8] * b[9 - i] + a[9 - i] * b[8])
					end
				end
			end]=]
			--c[5] = -c[5]
			--c[7] = -c[7]
			for i = 1, 8 do
				c[8] += a[i] * b[9 - i]
			end

			a[6] = -a[6]
			b[6] = -b[6]
			
			return setmetatable(c, Three)
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
		elseif typeof(t2) == "Vector3" then
			return three(t2.X + t1[1], t2.Y + t1[2], t2.Z + t1[3], t1[4], t1[5], t1[6], t1[7], t1[8])
		elseif typeof(t2) == "Vector2" then
			return three(t2.X + t1[1], t2.Y + t1[2], t1[3], t1[4], t1[5], t1[6], t1[7], t1[8])
		end
		return three(t1[1] + t2[1], t1[2] + t2[2], t1[3] + t2[3], t1[4] + t2[4], t1[5] + t2[5], t1[6] + t2[6], t1[7] + t2[7], t1[8] + t2[8])
	end,
	__sub = function(t1, t2)
		if typeof(t1) == "number" then
			return three(t1 - t2[1], t2[2], t2[3], t2[4], t2[5], t2[6], t2[7], t2[8])
		elseif typeof(t2) == "number" then
			return three(t1[1] - t2, t1[2], t1[3], t1[4], t1[5], t1[6], t1[7], t1[8])
		elseif typeof(t2) == "Vector3" then
			return three(t1[1] - t2.X, t1[2] - t2.Y, t1[3] - t2.Z, t1[4], t1[5], t1[6], t1[7], t1[8])
		elseif typeof(t2) == "Vector2" then
			return three(t1[1] - t2.X, t1[2] - t2.Y, t1[3], t1[4], t1[5], t1[6], t1[7], t1[8])
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
				t_prod[basis_index[bxor(basis_list[i], basis_list[j])]] += EPSTBL[EPSLOOKUP_X[j]][EPSLOOKUP_Z[i]] * EPSTBL[EPSLOOKUP_X[j]][EPSLOOKUP_Y[i]] * EPSTBL[EPSLOOKUP_Y[j]][EPSLOOKUP_Z[i]] * t1[i] * t2[j]
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
Three.lib.join = Three.lib.regressive
Three.lib.meet = Three.lib.regressive
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
	if getmetatable(a) == Three then
		return a
	end
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
Three.lib.one = three(1, 0, 0, 0, 0, 0, 0, 0)

--[=[
return {
	three = setmetatable({}, {
		__index = function(t, k)
			return Three.lib[k]
		end,
		__call = function(_, ...)
			return three(...)
		end,
	})
}]=]
return setmetatable({}, {
	__index = function(t, k)
		return Three.lib[k]
	end,
	__call = function(_, ...)
		return three(...)
	end,
})