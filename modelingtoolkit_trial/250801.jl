using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

@mtkmodel MDK begin
    @parameters begin
        F = 100
        d = 1
        k = 1000
    end
    @variables begin
        x(t) = 0.0
        ẋ(t) = F/d
        ẍ(t) = 0.0
    end
    @equations begin
        ẋ ~ D(x)
        ẍ ~ D(ẋ)
        d*ẋ + k*x^1.5 ~ F
    end
end

using OrdinaryDiffEq
@mtkcompile mdk = MDK()
prob = ODEProblem(mdk, [], (0.0, 0.01), [])
sol = solve(prob)

using Plots
plot(sol; idxs=ẍ, xlabel="time [s]", ylabel="ẍ [m/s^2]")