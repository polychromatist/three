# three
Roblox Lua library for 3D multivectors, the geometric product and associated operations

### Warning: May Not Be Suitable for Production
- This code is superbly undertested: I have only verified it with bouts of trial-and-error. Use at your own risk.
- There is a concern with floating-point error here in that some operations that should yield a zero value in a component of a multivector, such as inverse, yields an epsilon value. This is resolvable by, for example, implementing a symbolic equation solver that would destroy cancelling terms before they are evaluated. (toeventuallydo)
