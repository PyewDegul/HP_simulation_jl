using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

value_vector = randn(10)
f_fun(t) = t >= 10 ? value_vector[end] : value_vector[Int(floor(t)) + 1]
@register_symbolic f_fun(t)

@mtkmodel FOLExternalFunction begin
    @parameters begin
        τ = 0.75 # parameters and their values
    end
    @variables begin
        x(t) = 0.0 # dependent variables and their initial conditions
        f(t)
    end
    @equations begin
        f ~ f_fun(t)
        D(x) ~ (f - x) / τ
    end
end

@mtkcompile fol_external_f = FOLExternalFunction()

prob = ODEProblem(fol_external_f, [], (0.0, 10.0))
sol = solve(prob)
plot(sol, idxs = [fol_external_f.x, fol_external_f.f])