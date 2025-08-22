using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

@connector MechanicalPort begin
    v(t)
    f(t), [connect = Flow]
end

@mtkmodel Mass begin
    @parameters begin
        m = 10
    end
    @variables begin
        v(t)
        f(t)
    end
    @components begin
        port = MechanicalPort()
    end
    @equations begin
        # connectors
        port.v ~ v
        port.f ~ f

        # physics
        f ~ m*D(v)
    end
end

@mtkmodel Damper begin
    @parameters begin
        d = 1
    end

    @variables begin
        v(t)
        f(t)
    end

    @components begin
        port_a = MechanicalPort()
        port_b = MechanicalPort()
    end
    @equations begin
        # connectors
        (port_a.v - port_b.v) ~ v
        port_a.f ~ +f
        port_b.f ~ -f

        # physics
        f ~ d*v
    end
end

@mtkmodel Spring begin
    @parameters begin
        k = 100
    end

    @variables begin
        x(t)
        v(t)
        f(t)
    end

    @components begin
        port_a = MechanicalPort()
        port_b = MechanicalPort()
    end

    @equations begin
        # Derivatives
        D(x) ~ v

        # connectors
        (port_a.v - port_b.v) ~ v
        port_a.f ~ +f
        port_b.f ~ -f

        # physics
        f ~ k*x
    end
end

@mtkmodel Reference begin
    @parameters begin

    end
    @variables begin
        f(t)
    end
    @components begin
        port = MechanicalPort(;f,v=0)
    end
    @equations begin
        # connectors
        port.v ~ 0
        port.f ~ -f
    end
end

@mtkmodel ConstantForce begin
    @parameters begin
        f
    end
    @variables begin
        v(t)
    end
    @components begin
        port = MechanicalPort()
    end
    @equations begin
        # connectors
        port.v ~ v
        port.f ~ -f
    end
end

@mtkmodel System_msd_ref begin
    @parameters begin
        v=0
        x=0
        m=100
        d=10
        k=1000
        f=1
    end
    @components begin
        mass = Mass(;v,m)
        damper = Damper(;v, d)
        spring = Spring(;v, k, x)
        ref = Reference()
        force = ConstantForce(;v,f)
    end
    @equations begin
        connect(mass.port, damper.port_a, spring.port_a, force.port)
        connect(damper.port_b, spring.port_b, ref.port)
    end
end

@mtkmodel MassSpringDamper begin
    @parameters begin
        m
        k
        d
        v
        x
    end
    @components begin
        port_a = MechanicalPort()
        port_b = MechanicalPort()
        mass = Mass(;v,m)
        damper = Damper(;v, d)
        spring = Spring(;v, k, x)
    end
    @equations begin
        connect(mass.port, damper.port_a, spring.port_a, port_a)
        connect(damper.port_b, spring.port_b, port_b)
    end
end

@mtkmodel System_3msd begin
    @parameters begin
        v = 0
        x = 0
    end

    @components begin
        msd1 = MassSpringDamper(m=10, d=1, k=1000, v, x)
        msd2 = MassSpringDamper(m=20, d=2, k=2000, v, x)
        msd3 = MassSpringDamper(m=30, d=3, k=3000, v, x)
        ref = Reference()
        force = ConstantForce(;f=1, v)
    end

    @equations begin
        connect(force.port, msd1.port_a)
        connect(msd1.port_b, msd2.port_a)
        connect(msd2.port_b, msd3.port_a)
        connect(msd3.port_b, ref.port)
    end
end

@mtkbuild sys = System_3msd()

prob = ODEProblem(sys, [], (0, 2))
sol = solve(prob)
plot(sol; idxs=[sys.msd1.spring.x, sys.msd2.spring.x, sys.msd3.spring.x])