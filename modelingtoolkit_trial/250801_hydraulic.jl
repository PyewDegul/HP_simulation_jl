using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using DifferentialEquations
using Symbolics
using Plots

# parameters -------
pars = @parameters begin
    r0 = 876 #kg/m^3
    beta = 1.2e9 #Pa
    A = 0.01 #m²
    x0 = 1.0 #m
    M = 10_000 #kg
    g = 9.807 #m/s²
    amp = 5e-2 #m
    f = 15 #Hz
end

dt = 1e-4 #s
t_end = 0.2 #s
time = 0:dt:t_end

x_fun(t,amp,f) = amp*sin(2π*t*f) + x0
xdot_expr = expand_derivatives(D(x_fun(t,amp,f)))  # dx/dt
xdot_fun = build_function(xdot_expr, t, amp, f; expression=false)

vars = @variables begin
    x(t) = x0
    xdot(t)
    xddot(t)
    p(t) = M*g / A #Pa
    mdot(t)
    r(t)
    rdot(t)
end

function get_base_equations(density_type)
    eqs = [
        D(x) ~ xdot
        D(xdot) ~ xddot
        D(r) ~ rdot

        r ~ r0 * (1 + p / beta)
        mdot ~ rdot * x* A + (density_type) * xdot * A
        M*xddot ~ p*A - M*g
    ]
    return eqs
end
eqs_mdot1 = [
    get_base_equations(r0)...
    mdot ~ xdot_fun(t, amp, f)*A*r
]
eqs_mdot2 = [
    get_base_equations(r)...
    mdot ~ xdot_fun(t, amp, f)*A*r
]
eqs_x = [
    get_base_equations(r)...
    x ~ x_fun(t,amp,f)
]

@mtkbuild odesys_x = ODESystem(eqs_x, t, vars, pars)
@mtkbuild odesys_mdot1 = ODESystem(eqs_mdot1, t, vars, pars)
@mtkbuild odesys_mdot2 = ODESystem(eqs_mdot2, t, vars, pars)

prob_x = ODEProblem(odesys_x, [], (0, t_end))

u0 = [sol_x[s][1] for s in unknowns(odesys_mdot1)]
prob_mdot1 = ODEProblem(odesys_mdot1, u0, (0, t_end))
prob_mdot2 = ODEProblem(odesys_mdot2, u0, (0, t_end))

sol_x = solve(prob_x; saveat=time)
sol_mdot1 = solve(prob_mdot1; initializealg=NoInit());
sol_mdot2 = solve(prob_mdot2; initializealg=NoInit());

plot(sol_x; idxs=x, label="solution", ylabel="x")
plot!(sol_mdot1; idxs=x, label="case 1: r₀")
plot!(sol_mdot2; idxs=x, label="case 2: r")