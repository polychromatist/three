# three
Roblox Lua library for 3D multivectors, the geometric product and associated operations

### what/why
The space of 3D multivectors is a space of regular 3D vectors and an assortment of objects, all of which become transformations under the **geometric product**, a useful way to multiply vectors. Projections, rotations, reflections are defined with this product in a way that is practical, but also illuminating, suggestive, and privy to tinkering.

The language of geometric algebra is quickly taking the stage as a universal representation of transformations in space.

Please take some time to understand these captivating objects: https://www.youtube.com/watch?v=60z_hpEAtD8

### usage
Place contents of `three.lua` into ReplicatedStorage\
`local three = require(game:GetService("ReplicatedStorage").three)`
To specify a multivector directly:
`local m1 = three(a, x, y, z, xy, xz, yz, im)`
- `a` describes the scalar component `m1.scalar` (returns a number)
- `x, y, z` describes the vector component `m1.vector` (returns a 3-tuple of numbers)
- `xy, xz, yz` describes the pseudovector component `m1.pseudovector` (similarly a 3-tuple)
- `im` describes the pseudoscalar component `m1.pseudoscalar` (number)

_If you have two multivectors `m1`, `m2`, you can..._
- Add / subtract them
- Take their geometric product: `m1 * m2` (m2 can also be a scalar, Vector3, or Vector2)
- Take their inner product: `m1 .. m2` OR `m1:inner(m2)` (wrapper for `0.5 * (m1 * m2 + m2 * m1)`)
- Take their outer (wedge) product: `m1 ^ m2` OR `m1:outer(m2)` (wrapper for `0.5 * (m1 * m2 - m2 * m1)`)
- Take a quotient (if m2 is invertible): `m1 / m2` (wrapper for `m1 * m2.inverse`)

_If you have one multivector `m1`, you can..._
- Take its reverse conjugate `m1.reversion`, which is the multivector such that, if `m1 = v1 * v2 * ... * v_n` for vectors v_i, then `m1.reversion = v_n * v_(n-1) * ... * v_1`
- Take its Clifford conjugate `m1.bar`, which is useful to provide a value analogous to magnitude (by `m1 * m1.bar`)
- Take its "magnitude" `m1.magnitude`, which is vector magnitude if the multivector is a vector, but is defined by `sqrt(sqrt(b.scalar ^ 2 + b.pseudoscalar ^ 2))`, where `b = m1 * m1.bar`
- Take its "unit" multivector `m1.unit`, simply `m1 / m1.magnitude`
- Take its inverse `m1.inverse`, defined by `m1 * m1.inverse = id = unit scalar = three(1)`\
See an interesting paper on generalized multivector inverse: https://arxiv.org/abs/1712.05204v2 "INVERSE OF MULTIVECTOR: BEYOND P+Q=5 THRESHOLD" (A. ACUS AND A. DARGYS)

### not confirmed
- Dual `m1.dual`, which I define by `m1 = m1.dual * I`, where `I = three(0,0,0,0,0,0,0,1) = three.im`
- Reflections `m1:reflect(m2)`, regressive product `m1 % m2` OR `m1:regressive(m2)`\

### not yet implemented
- rotations, multivector SQRT\
Motivating paper on generalized multivector SQRT: https://arxiv.org/abs/2003.06873v1 "Square root of a multivector of Clifford algebras in 3D: A game with signs" (A. Acus, A. Dargys)

### Warning: May Not Be Suitable for Production
- This code is superbly undertested: I have only verified it with bouts of trial-and-error. Use at your own risk.
- There is a concern with floating-point error here in that some operations that should yield a zero value in a component of a multivector, such as inverse, yields an epsilon value. This is resolvable by, for example, implementing a symbolic equation solver that would destroy cancelling terms before they are evaluated. (toeventuallydo)
