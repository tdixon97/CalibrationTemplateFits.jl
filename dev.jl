using Pkg
Pkg.instantiate()
Pkg.test()
using Aqua
Aqua.test_all(MyPackage)
