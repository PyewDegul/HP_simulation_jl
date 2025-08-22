using ModelingToolkit
using ModelingToolkitStandardLibrary
using ModelingToolkit: t_nounits as t, D_nounits as D

@connector FluidPort begin
    p(t)                      # 압력 [Pa]
    ṁ(t),                     [connect = Flow]  # 질량유량 (+: 포트 → 컴포넌트)
    h(t)                      # 엔탈피 [J/kg]
end

############################
# 1. Source:  Pin, hin 지정
############################
@mtkmodel PHSource begin
    @parameters Pin    # Pa
    @parameters hin    # J/kg
    port = FluidPort()
    @equations begin
        port.p  ~ Pin
        port.h  ~ hin
    end
end

############################
# 2. Sink:  Pout만 고정
############################
@mtkmodel PSink begin
    @parameters Pout = 1.0e6       # Pa
    port = FluidPort()
    @equations begin
        port.p  ~ Pout
        # port.h ~ 
        # h 는 유입 유체 엔탈피와 같아야 하므로 별도 방정식 X
        # ṁ 는 연결된 상대가 결정
    end
end

############################
# 3. Valve:  mdot, Hdot 연산
############################

@mtkmodel SimpleValve begin
    @parameters Kv = 1.0e-7        # 밸브 계수 [kg/(s·Pa)]
    @variables  del_p(t) ṁ(t)
    port_a = FluidPort(); 
    port_b = FluidPort()
    @equations begin
        del_p ~ port_a.p - port_b.p
        ṁ  ~  Kv * del_p            # 단순 선형 관계
        port_a.ṁ ~  ṁ
        port_b.ṁ ~ -ṁ
        port_a.h  ~ port_b.h       # 단열·엔탈피 손실 무시
    end
end

############################
# 4. System
############################
@mtkmodel total_model begin
    @parameters begin
    end
    @components begin
        source = PHSource(Pin = 1.0e6, hin = 1000.0)  # 예시 값
        sink = PSink(Pout = 1.0e5)  # 예시 값
        valve = SimpleValve(Kv = 1.0e-7)
    end
    @equations begin
        connect(source.port, valve.port_a)
        connect(valve.port_b, sink.port)
    end
end

@mtkbuild total_sys  = total_model()