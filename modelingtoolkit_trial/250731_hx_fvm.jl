using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using CoolProp

@variables t
@parameters UA V

@connector FluidPort begin        # Across / Stream / Flow 지정
    p::Real
    h::Real
    @variables ṁ(t)::Real [connect=Flow]
end

@connector HeatPort begin
    T::Real
    @variables Q̇(t)::Real [connect=Flow]   # → 연결된 두 Port Q̇ 합 = 0
end

temperature_ph(p,h) = PropsSI("T","P",p,"H",h,"R134a")
@register_symbolic temperature_ph(p,h)


@component FVCell1D begin
    in::FluidPort
    out::FluidPort
    wall::HeatPort
    @parameters V UA                       # 제어부피, 열전달계수
    @variables rho(t) h(t)
    T = temperature_ph(out.p,h)            # Gibbs 관계 포함
    eqs = quote
        Differential(rho)*V              ~  in.ṁ - out.ṁ
        rho*V*Differential(h) + h*V*Differential(rho)
                                        ~  in.ṁ*in.h - out.ṁ*h + wall.Q̇
        wall.Q̇                        ~  UA*(T - wall.T)
        out.h                          ~  h
    end
end

@component Source begin      # 간단한 정압/정엔탈피 소스
    port::FluidPort
    @parameters p0 h0 ṁ0
    eqs = quote
        port.p  ~ p0
        port.h  ~ h0
        port.ṁ  ~ -ṁ0        # 양수를 ‘시스템 → 소스’ 로 정의했을 때
    end
end

@component Sink begin
    port::FluidPort
    @parameters p_set
    eqs = quote
        port.p  ~ p_set
        port.ṁ  ~ 0          # 자유유량(압력만 고정)
    end
end