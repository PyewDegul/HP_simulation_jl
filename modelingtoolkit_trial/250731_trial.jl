using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

@mtkmodel FOL begin
    @parameters begin
        τ = 3.0 # parameters and their values
    end
    @variables begin
        x(t) = 0.0 # dependent variables and their initial conditions
        RHS(t)
    end
    @equations begin
        RHS ~ (1 - x) / τ
        D(x) ~ RHS
    end
end


using OrdinaryDiffEq
using Plots

@mtkcompile fol = FOL()
prob = ODEProblem(fol, [], (0.0, 10.0))
sol = solve(prob)
plot(sol, idxs = [fol.x, fol.RHS])
