using Aqua
@testset "Aqua" begin
    Aqua.test_all(Biomodelling; ambiguities = false, deps_compat = (check_extras = false,), piracies = false, stale_deps = (ignore = [:Printf, :LinearAlgebra],))
end
