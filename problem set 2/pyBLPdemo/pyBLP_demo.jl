using Conda
for p ∈ ["numpy" "scipy" "sympy" "patsy" "pyhdfe"]
    Conda.add(p)
end
Conda.pip_interop(true)
Conda.pip("install", "pyblp")
Conda.pip_interop(false)

using PyCall
pyblp = pyimport("pyblp")
