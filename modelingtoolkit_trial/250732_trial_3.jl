@mtkmodel FOLUnconnectedFunction begin
    @parameters begin
        τ # parameters
    end
    @variables begin
        x(t) # dependent variables
        f(t)
        RHS(t)
    end
    @equations begin
        RHS ~ f
        D(x) ~ (RHS - x) / τ
    end
end
@mtkmodel FOLConnected begin
    @components begin
        fol_1 = FOLUnconnectedFunction(; τ = 2.0, x = -0.5)
        fol_2 = FOLUnconnectedFunction(; τ = 4.0, x = 1.0)
    end
    @equations begin
        fol_1.f ~ 1.5
        fol_2.f ~ fol_1.x
    end
end

using BenchmarkTools

@mtkcompile connected = FOLConnected()
prob = ODEProblem(connected, [], (0.0, 10.0))
prob_an = ODEProblem(connected, [], (0.0, 10.0); jac = true)
prob_sparse = ODEProblem(connected, [], (0.0, 10.0); jac = true, sparse = true)

@btime solve(prob_an, Rodas4());
