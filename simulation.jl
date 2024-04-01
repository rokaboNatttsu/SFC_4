using StatsPlots

T = 150

#   パラメータ設定
α1, α2, α3, α4, α5, α6 = 0.9, 0.1, 0.1, 0.02, 0.5, 0.2
β1, β2 = 0.05, 0.5  #   0.05, 0.5
γ1, γ2 = 0.015, 0.02
δ1, δ2 = fill(0.02, T), 0.3  #   0.02, 0.3
δ1[end-40:end] .= 0.02
ϵ1, ϵ2, ϵ3, ϵ4 = 0.5, 0.7, 0.8, 0.05    #   0.5, 0.7, 0.8, 0.05
ζ1, ζ2, ζ3 = 1.0, fill(300.0, T), 0.023   #   1.0, 300.0, 0.023     (ζ2[1] = k[1])
for t = 2:T
    ζ2[t] = ζ2[t-1]*(1 + ζ3*abs(randn()))
end
θ1, θ2 = 0.4, 0.1   #   0.4, 0.1
ι1, ι2, ι3, ι4, ι5, ι6, ι7 = 0.1, 1.0, 0.9, 0.5, 10.0, 0.5, 0.5
λe = 0.5
μ1, μ2, μ3 = 0.3, 0.1, 0.5
τ1, τ2 = 0.2, 0.2
uT = 0.8
G0 = 100.0
r = 0.01

#   配列定義
G, GD, g, gD = zeros(T), zeros(T), zeros(T), zeros(T)
Wf, Wg, Wb, W, wg = zeros(T), zeros(T), zeros(T), zeros(T), zeros(T)
Cw, Ci, C, cw, ci, c, CwD, CiD = zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T)
K, k = zeros(T), zeros(T)
Mw, Mi, Mf, M = zeros(T), zeros(T), zeros(T), zeros(T)
Ei, E, Eb, ei, e, eb = zeros(T), zeros(T), zeros(T), zeros(T), fill(100.0, T), zeros(T)
Lw, Lf, L, LfD = zeros(T), zeros(T), zeros(T), zeros(T)
Hw, Hb, H = zeros(T), zeros(T), zeros(T)
NWw, NWi, NWf, NWb, NWg = zeros(T), zeros(T), zeros(T), zeros(T), zeros(T)
NLw, NLi, NLf, NLb, NLg = zeros(T), zeros(T), zeros(T), zeros(T), zeros(T)
Tef, Tei, Tew, Te = zeros(T), zeros(T), zeros(T), zeros(T)
Tii, Tiw, Ti = zeros(T), zeros(T), zeros(T)
Tew, Tei, Tef, Te = zeros(T), zeros(T), zeros(T), zeros(T)
Tff, Tfb, Tf = zeros(T), zeros(T), zeros(T)
Πi, Π, Πf, Πb = zeros(T), zeros(T), zeros(T), zeros(T)
I, i = zeros(T), zeros(T)
p, pe = zeros(T), zeros(T)
u, Δu, ue = zeros(T), zeros(T), zeros(T)
Δk, ΔK = zeros(T), zeros(T), zeros(T), zeros(T)
ΔMw, ΔMi, ΔMf, ΔM = zeros(T), zeros(T), zeros(T), zeros(T)
ΔLw, ΔLf, ΔL = zeros(T), zeros(T), zeros(T)
Δei, Δeb, ΔEi, ΔEb, ΔE = zeros(T), zeros(T), zeros(T), zeros(T), zeros(T)
ΔHw, ΔHb, ΔH = zeros(T), zeros(T), zeros(T)
Δp, Δpe = zeros(T), zeros(T), zeros(T)
Le, Ebe = zeros(T), zeros(T)
We, Tiwe, Tewe = zeros(T), zeros(T), zeros(T)
Πie, Πbe, Ee, Tiie, Teie = zeros(T), zeros(T), zeros(T), zeros(T), zeros(T)
NWie, ΔNWi = zeros(T), zeros(T)
ΔW, ΔTiw, ΔTew = zeros(T), zeros(T), zeros(T)
ΔΠi, ΔΠb, ΔTii, ΔTei = zeros(T), zeros(T), zeros(T), zeros(T)
ΔCi, ΔCw = zeros(T), zeros(T)
EiT, EbT = zeros(T), zeros(T)

#   初期値設定
K[1] = 300
Mw[1], Mi[1], Mf[1], M[1] = 400.0, 50.0, 150.0, 600.0
Ei[1], E[1], Eb[1] = 50.0, 100.0, 50.0
ei[1], e[1], eb[1] = 50.0, 100.0, 50.0
Lw[1], Lf[1], L[1] = 200.0, 0.0, 200.0
Hw[1], Hb[1], H[1] = 20.0, 580.0, 600.0
G[1], gD[1] = G0, G0
Wf[1], Wg[1], W[1], wg[1] = 1.0*G0, 0.5*G0, 1.5*G0, 0.5*G0
p[1], pe[1] = 1.0, 1.0
k[1] = K[1]/p[1]
CwD[1], CiD[1] = 1.0*G0, 0.1*G0
Cw[1], Ci[1] = CwD[1]/max(1, (CwD[1] + CiD[1])/(p[1]*ζ1*k[1])), CiD[1]/max(1, (CwD[1] + CiD[1])/(p[1]*ζ1*k[1]))
cw[1], ci[1] = Cw[1]/p[1], Ci[1]/p[1]
C[1], c[1] = Cw[1] + Ci[1], cw[1] + ci[1]
u[1] = min(1, (C[1] + Ci[1])/(p[1] * ζ1 * k[1]))
ue[1] = u[1]
Π[1] = C[1] + I[1] + G[1] - Wf[1] - Tef[1] - Tff[1] - r*Lf[1]
Πi[1] = max(0, θ1*(Π[1]-I[1])*Ei[1]/E[1] + θ2*(Mf[1]-Lf[1]))
Πb[1] = max(0, θ1*(Π[1]-I[1])*Eb[1]/E[1] + θ2*(Mf[1]-Lf[1]))
Πf[1] = Π[1] - Πi[1] - Πb[1]
Tiw[1], Tii[1], Tei[1], Tew[1] = τ1*W[1], τ1*Πi[1], γ2*(Mi[1] + Ei[1]), γ2*(Mw[1] + Hw[1])
Le[1] = L[1]
Πie[1], Πbe[1], Ee[1] = Πi[1], Πb[1], E[1]
We[1], Tiwe[1], Tewe[1] = W[1], Tiw[1], Tew[1]
Tiie[1], Teie[1] = Tii[1], Tei[1]
EiT[1], EbT[1] = Ei[1], Eb[1]

#   会計的整合性を満たすように、残りの変数の値を決める
#   初期値の計算においては、使わない会計恒等式も多い
NLw[1] = -Cw[1] + W[1] - Tiw[1] - Tew[1] - r*L[1]
ΔMw[1] = NLw[1] # ΔLw[1] = ΔHw[1] = 0ということにする
NLi[1] = Πi[1] - Ci[1] - Tii[1] - Tei[1]  # i[1] = 0を仮定
ΔMi[1] = NLi[1]  # Δei[1] = 0ということにする
NLf[1] = Πf[1]  # i[1] = 0を仮定
ΔMf[1] = NLf[1] # ΔLf[1] = 0を仮定
Ti[1] = Tiw[1] + Tii[1]
Te[1] = Tew[1] + Tei[1] + Tef[1]
ΔM[1] = ΔMw[1] + ΔMi[1] + ΔMf[1]
ΔL[1] = ΔLw[1] + ΔLf[1]
NLb[1] = -Tfb[1] + Πb[1] + r*L[1]
ΔHb[1] = NLb[1] + ΔM[1] - ΔL[1] - pe[1]*Δeb[1]
ΔH[1] = ΔHw[1] + ΔHb[1]
NLg[1] = -G[1] - Wg[1] + Ti[1] + Te[1] + Tf[1]

NWw[1] = Mw[1] - Lw[1] + Hw[1]
NWi[1] = Mi[1] + Ei[1]
NWf[1] = K[1] + Mf[1] - E[1] - Lf[1]
NWb[1] = -M[1] + Eb[1] + L[1] + Hb[1]
NWg[1] = -H[1]

#   会計的整合性の確認
println("------initial------")
println("-G-Wg+Ti+Te+Tf+ΔH=", -G[1]-Wg[1]+Ti[1]+Te[1]+Tf[1]+ΔH[1])
println("NLw+NLi+NLf+NLb+NLg=", NLw[1]+NLi[1]+NLf[1]+NLb[1]+NLg[1])
println("-M+Mw+Mi+Mf=", -M[1]+Mw[1]+Mi[1]+Mf[1])
println("-E+Ei+Eb=", -E[1]+Ei[1]+Eb[1])
println("-L+Lw+Lf=", -L[1]+Lw[1]+Lf[1])
println("-H+Hw+Hf=", -H[1]+Hw[1]+Hb[1])
println("NWw+NWi+NWf+NWb+NWg-K=", NWw[1]+NWi[1]+NWf[1]+NWb[1]+NWg[1]-K[1])

function btw(A, B, C)
    if A > C
        return nothing
    else
        return min(max(A, B), C)
    end
end

function run()
    #   シミュレーション実行
    for t = 2:T
        ue[t] = ((1 - λe)*ue[t-1] + λe*u[t-1])*u[t-1]/(u[t-1]-Δu[t-1])
        p[t] = p[t-1]*exp(μ3*(min(1, max(ue[t], (CwD[t-1] + CiD[t-1] + GD[t-1])/(p[t-1]*ζ2[t-1])) - uT)))
        Δp[t] = p[t] - p[t-1]
        Wf[t] = (1.0 - ϵ1)*Wf[t-1] + ϵ1*max(0.0, (Wf[t-1] + Π[t-1] - I[t-1])*ϵ3)
        Wb[t] = (1.0 - ϵ1)*Wb[t-1] + ϵ1*max(0.0, ϵ2*(Πb[t-1] + r*L[t-1]) + ϵ4*(Hb[t-1] - M[t-1]))
        wg[t] = wg[t-1]*(1.0 + δ1[t] - δ2*(p[t]-p[t-1])/p[t-1])
        Wg[t] = p[t]*wg[t-1]
        W[t] = Wf[t] + Wg[t] + Wb[t]
        ΔW[t] = W[t] - W[t-1]
        gD[t] = gD[t-1]*(1.0 + δ1[t] - δ2*(p[t]-p[t-1])/p[t-1])
        GD[t] = p[t]*gD[t]
        We[t] = ((1 - λe)*We[t-1] + λe*W[t-1])*W[t-1]/(W[t-1]-ΔW[t-1])
        Tiwe[t] = ((1 - λe)*Tiwe[t-1] + λe*Tiw[t-1])*Tiw[t-1]/(Tiw[t-1]-ΔTiw[t-1])
        Tewe[t] = ((1 - λe)*Tewe[t-1] + λe*Tew[t-1])*Tew[t-1]/(Tew[t-1]-ΔTew[t-1])
        CwD[t] = (1 - α5)*CwD[t-1] + α5*max(0, α1*(We[t]-Tiwe[t]-Tewe[t]-r*Lw[t-1]) + α2*(Mw[t-1] + Hw[t-1] - Lw[t-1]))*Cw[t-1]/(Cw[t-1]-ΔCw[t-1])
        Πie[t] = ((1 - λe)*Πie[t-1] + λe*Πi[t-1])*Πi[t-1]/(Πi[t-1]-ΔΠi[t-1])
        Tiie[t] = ((1 - λe)*Tiie[t-1] + λe*Tii[t-1])*Tii[t-1]/(Tii[t-1]-ΔTii[t-1])
        Teie[t] = ((1 - λe)*Teie[t-1] + λe*Tei[t-1])*Tei[t-1]/(Tei[t-1]-ΔTei[t-1])
        CiD[t] = (1 - α5)*CiD[t-1] + α5*max(0, α3*(Πie[t]-Tiie[t]-Teie[t]) + (α4 + α6*Ei[t-1]/E[t-1])*Mi[t-1])*Ci[t-1]/(Ci[t-1]-ΔCi[t-1])
        Cw[t] = CwD[t]/max(1, (CwD[t] + CiD[t] + GD[t])/(p[t]*min(ζ1*k[t-1], ζ2[t])))
        ΔCw[t] = Cw[t] - Cw[t-1]
        Ci[t] = CiD[t]/max(1, (CwD[t] + CiD[t] + GD[t])/(p[t]*min(ζ1*k[t-1], ζ2[t])))
        ΔCi[t] = Ci[t] - Ci[t-1]
        C[t] = Cw[t] + Ci[t]
        c[t] = C[t]/p[t]
        G[t] = GD[t]/max(1, (CwD[t] + CiD[t] + GD[t])/(p[t]*min(ζ1*k[t-1], ζ2[t])))
        g[t] = G[t]/p[t]
        u[t] = (c[t] + g[t])/(ζ1*k[t-1])
        Δu[t] = u[t] - u[t-1]
        Tef[t] = γ1*K[t-1]
        Tei[t] = γ2*(Mi[t-1] + Ei[t-1])
        ΔTei[t] = Tei[t] - Tei[t-1]
        Tew[t] = γ2*(Mw[t-1] + Hw[t-1])
        Te[t] = Tef[t] + Tei[t] + Tew[t]
        i[t] = btw(0, (u[t-1] - uT)*k[t-1] + β1*k[t-1], max(0, β2*(Mf[t-1] - Lf[t-1])/p[t]))
        I[t] = p[t]*i[t]
        k[t] = (1 - β1)k[t-1] + i[t]
        Δk[t] = k[t] - k[t-1]
        K[t] = p[t-1]*k[t-1] + Δp[t]*k[t-1] + p[t]*Δk[t]
        ΔK[t] = K[t] - K[t-1]
        Tii[t] = τ1*Πi[t-1]
        ΔTii[t] = Tii[t] - Tii[t-1]
        Tiw[t] = τ1*W[t-1]
        ΔTiw[t] = Tiw[t] - Tiw[t-1]
        Ti[t] = Tii[t] + Tiw[t]
        Tff[t] = τ2*(C[t] + I[t] + G[t] - Wf[t] - Tef[t] - r*Lf[t-1])
        Π[t] = C[t] + I[t] + G[t] - Wf[t] - Tef[t] - Tff[t] - r*Lf[t-1]
        Πi[t] = max(0, (θ1*(Π[t] - I[t]) + θ2*(Mf[t-1] - Lf[t-1]))*Ei[t-1]/E[t-1])
        ΔΠi[t] = Πi[t] - Πi[t-1]
        Πb[t] = max(0, (θ1*(Π[t] - I[t]) + θ2*(Mf[t-1] - Lf[t-1]))*Eb[t-1]/E[t-1])
        ΔΠb[t] = Πb[t] - Πb[t-1]
        Πf[t] = Π[t] - Πi[t] - Πb[t]
        Tfb[t] = max(0, τ2*(Πb[t] + r*L[t-1] - Wb[t] + Δpe[t]*eb[t-1]))
        Tf[t] = Tff[t] + Tfb[t]
        NLw[t] = -Cw[t] + W[t] - Tiw[t] - Tew[t] - r*Lw[t-1]
        NLi[t] = -Ci[t] - Tii[t] - Tei[t] + Πi[t]
        NLf[t] = -I[t] + Πf[t]
        NLb[t] = -Wb[t] - Tfb[t] + Πb[t] + r*L[t-1]
        NLg[t] = -G[t] - Wg[t] + Ti[t] + Te[t] + Tf[t]
        Hw[t] = ι1*Cw[t]
        ΔHw[t] = Hw[t] - Hw[t-1]
        Lw[t] = ι2*W[t]
        ΔLw[t] = Lw[t] - Lw[t-1]
        ΔMw[t] = NLw[t] - ΔHw[t] + ΔLw[t]
        Mw[t] = Mw[t-1] + ΔMw[t]
        if Πi[t-1]-ΔΠi[t-1] == 0.0
            Πie[t] = 0.0
        else
            Πie[t] = ((1 - λe)*Πie[t-1] + λe*Πi[t-1])*Πi[t-1]/(Πi[t-1]-ΔΠi[t-1])
        end
        if Πb[t-1]-ΔΠb[t-1] == 0.0
            Πbe[t] = 0.0
        else
            Πbe[t] = ((1 - λe)*Πbe[t-1] + λe*Πb[t-1])*Πb[t-1]/(Πb[t-1]-ΔΠb[t-1])
        end
        Ee[t] = ((1 - λe)*Ee[t-1] + λe*E[t-1])*E[t-1]/(E[t-1]-ΔE[t-1])
        Le[t] = ((1 - λe)*Le[t-1] + λe*L[t-1])*L[t-1]/(L[t-1]-ΔL[t-1])
        NWie[t] = ((1 - λe)*NWie[t-1] + λe*NWi[t-1])*NWi[t-1]/(NWi[t-1]-ΔNWi[t-1])
        ΔNWi[t] = NWi[t] - NWi[t-1]
        EiT[t], EbT[t] = ι3*NWie[t], (1.0 - ι7)*Eb[t-1] + ι7*((Πie[t] + Πbe[t])/Ee[t]*Le[t]/r)
        #TODO   もうちょっとどうにかならんか？ei,eb,peをめぐる行動方程式
        pe[t] = (EiT[t] + EbT[t])/e[t]
        Δpe[t] = pe[t] - pe[t-1]
        E[t] = pe[t]*e[t]
        ΔE[t] = E[t] - E[t-1]
        NWw[t] = NWw[t-1] + NLw[t]
        NWi[t] = NWi[t-1] + NLi[t] + Δpe[t]*ei[t-1]
        NWf[t] = NWf[t-1] + NLf[t] - ΔE[t] + ΔK[t]
        NWb[t] = NWb[t-1] + NLb[t] + Δpe[t]*eb[t-1]
        NWg[t] = NWg[t-1] + NLg[t]
        Ei[t] = min(ι3*NWi[t], E[t])
        ei[t] = Ei[t]/pe[t]
        ΔEi[t] = Ei[t] - Ei[t-1]
        Δei[t] = ei[t] - ei[t-1]
        eb[t] = e[t] - ei[t]
        Δeb[t] = eb[t] - eb[t-1]
        Eb[t] = pe[t]*eb[t]
        ΔEb[t] = Eb[t] - Eb[t-1]
        ΔMi[t] = NLi[t] - pe[t]*Δei[t]
        Mi[t] = NWi[t] - Ei[t]
        LfD[t] = (1.0 - ι4)*LfD[t-1] + ι4 * max(0.0, I[t] + Wf[t] + Tef[t] + Tff[t] + r*Lf[t-1] - Mf[t-1] - NLf[t])
        Lf[t] = btw(0.0, LfD[t], max(0.0, ι5*NLf[t]))
        ΔLf[t] = Lf[t] - Lf[t-1]
        ΔMf[t] = ΔLf[t] + NLf[t]
        Mf[t] = Mf[t-1] + ΔMf[t]
        M[t] = Mw[t] + Mi[t] + Mf[t]
        ΔM[t] = ΔMw[t] + ΔMi[t] + ΔMf[t]
        L[t] = Lw[t] + Lf[t]
        ΔL[t] = ΔLw[t] + ΔLf[t]
        ΔHb[t] = NLb[t] + ΔM[t] - ΔL[t] - pe[t]*Δeb[t]
        Hb[t] = Hb[t-1] + ΔHb[t]
        ΔH[t] = ΔHw[t] + ΔHb[t]
        H[t] = H[t-1] + ΔH[t]
    end

    #   会計的整合性の確認
    println("--------final--------")
    println("-G-Wg+Ti+Te+Tf+ΔH=",-G[end]-Wg[end]+Ti[end]+Te[end]+Tf[end]+ΔH[end])
    println("NLw+NLi+NLf+NLb+NLg=",NLw[end]+NLi[end]+NLf[end]+NLb[end]+NLg[end])
    println("Δei+Δeb=",Δei[end]+Δeb[end])
    println("L-L_{-1}-ΔL=",L[end]-L[end-1]-ΔL[end])
    println("E-Ei-Eb=",E[end]-Ei[end]-Eb[end])
    println("NWw+NWi+NWf+NWb+NWg-K=",NWw[end]+NWi[end]+NWf[end]+NWb[end]+NWg[end]-K[end])
    println("-NWw+Mw-Lw+Hw=",-NWw[end]+Mw[end]-Lw[end]+Hw[end])
    println("-NWi+Mi+Ei=",-NWi[end]+Mi[end]+Ei[end])
    println("-NWf+K+Mf-E-Lf=",-NWf[end]+K[end]+Mf[end]-E[end]-Lf[end])
    println("-NWb-M+Eb+L+Hb=",-NWb[end]-M[end]+Eb[end]+L[end]+Hb[end])
    println("-NWg-H=",-NWg[end]-H[end])
end

run()

function plot_transition()
    plot(p[end-50:end]./p[end-50], label="p")
    savefig("figs/p.png")

    plot(pe[end-50:end]./p[end-50], label="pe")
    savefig("figs/pe.png")

    plot(NLw[end-50:end]./G[end-50], label="NLw")
    plot!(NLi[end-50:end]./G[end-50], label="NLi")
    plot!(NLf[end-50:end]./G[end-50], label="NLf")
    plot!(NLb[end-50:end]./G[end-50], label="NLb")
    plot!(NLg[end-50:end]./G[end-50], label="NLg")
    savefig("figs/NL.png")

    plot(u[end-50:end], label="u")
    savefig("figs/u.png")

    plot(W[end-50:end]./G[end-50], label="W")
    plot!(Wf[end-50:end]./G[end-50], label="Wf")
    plot!(Wb[end-50:end]./G[end-50], label="Wb")
    plot!(Wg[end-50:end]./G[end-50], label="Wg")
    savefig("figs/W.png")

    plot(Π[end-50:end]./G[end-50], label="Π")
    plot!(Πf[end-50:end]./G[end-50], label="Πf")
    plot!(Πb[end-50:end]./G[end-50], label="Πb")
    plot!(Πi[end-50:end]./G[end-50], label="Πi")
    savefig("figs/Π.png")

    plot(k[end-50:end]./g[end-50], label="k")
    plot!(i[end-50:end]./g[end-50], label="i")
    savefig("figs/k_and_i.png")

    plot(C[end-50:end]./G[end-50], label="C")
    plot!(Ci[end-50:end]./G[end-50], label="Ci")
    plot!(Cw[end-50:end]./G[end-50], label="Cw")
    savefig("figs/nominal_C.png")

    plot(ζ2[end-50:end], label="ζ2")
    savefig("figs/ζ2.png")

    plot(C[end-50:end]./G[end-50], label="C")
    plot!(I[end-50:end]./G[end-50], label="I")
    plot!(G[end-50:end]./G[end-50], label="G")
    savefig("figs/nominal_Y.png")

    plot(c[end-50:end]./g[end-50], label="c")
    plot!(i[end-50:end]./g[end-50], label="i")
    plot!(g[end-50:end]./g[end-50], label="g")
    savefig("figs/y.png")

    plot(M[end-50:end]./G[end-50], label="M")
    plot!(Mw[end-50:end]./G[end-50], label="Mw")
    plot!(Mi[end-50:end]./G[end-50], label="Mi")
    plot!(Mf[end-50:end]./G[end-50], label="Mf")
    savefig("figs/M.png")

    plot(L[end-50:end]./G[end-50], label="L")
    plot!(Lw[end-50:end]./G[end-50], label="Lw")
    plot!(Lf[end-50:end]./G[end-50], label="Lf")
    savefig("figs/L.png")

    plot(NWw[end-50:end]./G[end-50], label="NWw")
    plot!(NWi[end-50:end]./G[end-50], label="NWi")
    plot!(NWf[end-50:end]./G[end-50], label="NWf")
    plot!(NWb[end-50:end]./G[end-50], label="NWb")
    plot!(NWg[end-50:end]./G[end-50], label="NWg")
    savefig("figs/NW.png")

    plot(Ei[end-50:end]./G[end-50], label="Ei")
    plot!(E[end-50:end]./G[end-50], label="E")
    plot!(Eb[end-50:end]./G[end-50], label="Eb")
    savefig("figs/nominal_E.png")

    plot(ei[end-50:end], label="ei")
    plot!(eb[end-50:end], label="eb")
    savefig("figs/e.png")

    plot(Hw[end-50:end]./G[end-50], label="Hw")
    plot!(Hb[end-50:end]./G[end-50], label="Hb")
    plot!(H[end-50:end]./G[end-50], label="H")
    savefig("figs/H.png")

    plot(Tiw[end-50:end]./G[end-50], label="Tiw")
    plot!(Tii[end-50:end]./G[end-50], label="Tii")
    plot!(Ti[end-50:end]./G[end-50], label="Ti")
    savefig("figs/Ti.png")

    plot(Tew[end-50:end]./G[end-50], label="Tew")
    plot!(Tei[end-50:end]./G[end-50], label="Tei")
    plot!(Tef[end-50:end]./G[end-50], label="Tef")
    plot!(Te[end-50:end]./G[end-50], label="Te")
    savefig("figs/Te.png")

    plot(Tfb[end-50:end]./G[end-50], label="Tfb")
    plot!(Tff[end-50:end]./G[end-50], label="Tff")
    plot!(Tf[end-50:end]./G[end-50], label="Tf")
    savefig("figs/Tf.png")

    plot(Πi[end-50:end]./G[end-50], label="Πi")
    plot!(Π[end-50:end]./G[end-50], label="Π")
    plot!(Πf[end-50:end]./G[end-50], label="Πf")
    plot!(Πb[end-50:end]./G[end-50], label="Πb")
    savefig("figs/Π.png")
end
plot_transition()

function save_mock_data()
    mock_C = deepcopy(C)
    mock_c = deepcopy(c)
    mock_I = deepcopy(I)
    mock_i = deepcopy(i)
    mock_p = deepcopy(p)
    mock_W = deepcopy(W)
    mock_G = deepcopy(G)
    mock_g = deepcopy(g)
    mock_pe = deepcopy(pe)
    mock_K = deepcopy(K)
    mock_k = deepcopy(k)
    mock_NLw = deepcopy(NLw)
    mock_NLi = deepcopy(NLi)
    mock_NLf = deepcopy(NLf)
    mock_NLb = deepcopy(NLb)
    mock_NLg = deepcopy(NLg)
    mock_M = deepcopy(M)
    mock_H = deepcopy(H)
    mock_L = deepcopy(L)
    mock_Te = deepcopy(Te)
    mock_Tf = deepcopy(Tf)
    mock_Ti = deepcopy(Ti)
    mock_u = deepcopy(u)
    mock_Π = deepcopy(Π)
    mock_NWw = deepcopy(NWw)
    mock_NWi = deepcopy(NWi)
    mock_NWf = deepcopy(NWf)
    mock_NWb = deepcopy(NWb)
    mock_NWg = deepcopy(NWg)
    return mock_C, mock_c, mock_I, mock_i, mock_p, mock_W, mock_G, mock_g, mock_pe, mock_K, mock_k, mock_NLw, mock_NLi, mock_NLf, mock_NLb, mock_NLg, mock_M, mock_H, mock_L, mock_Te, mock_Tf, mock_Ti, mock_u, mock_Π, mock_NWw, mock_NWi, mock_NWf, mock_NWb, mock_NWg
end
function plot_mock_and_control()
    plot(mock_C[end-50:end]./mock_G[end-50], label="mock")
    plot!(C[end-50:end]./G[end-50], label="control")
    savefig("figs/diff_nominal_C.png")
    plot(mock_c[end-50:end]./mock_g[end-50], label="mock")
    plot!(c[end-50:end]./g[end-50], label="control")
    savefig("figs/diff_c.png")
    plot(mock_I[end-50:end]./mock_G[end-50], label="mock")
    plot!(I[end-50:end]./G[end-50], label="control")
    savefig("figs/diff_nominal_I.png")
    plot(mock_i[end-50:end]./mock_g[end-50], label="mock")
    plot!(i[end-50:end]./g[end-50], label="control")
    savefig("figs/diff_i.png")
    plot(mock_p[end-50:end]./mock_p[end-50], label="mock")
    plot!(p[end-50:end]./p[end-50], label="control")
    savefig("figs/diff_p.png")
    plot(mock_W[end-50:end]./mock_G[end-50], label="mock")
    plot!(W[end-50:end]./G[end-50], label="control")
    savefig("figs/diff_nominal_W.png")
    plot(mock_W[end-50:end]./mock_p[end-50:end]./mock_G[end-50], label="mock")
    plot!(W[end-50:end]./p[end-50:end]./G[end-50], label="control")
    savefig("figs/diff_real_W.png")
    plot(mock_G[end-50:end]./mock_G[end-50], label="mock")
    plot!(G[end-50:end]./G[end-50], label="control")
    savefig("figs/diff_nominal_G.png")
    plot(mock_g[end-50:end]./mock_g[end-50], label="mock")
    plot!(g[end-50:end]./g[end-50], label="control")
    savefig("figs/diff_g.png")
    plot(mock_pe[end-50:end]./mock_pe[end-50], label="mock")
    plot!(pe[end-50:end]./pe[end-50], label="control")
    savefig("figs/diff_pe.png")
    plot(mock_K[end-50:end]./mock_G[end-50], label="mock")
    plot!(K[end-50:end]./G[end-50], label="control")
    savefig("figs/diff_nominal_K.png")
    plot(mock_k[end-50:end]./mock_g[end-50], label="mock")
    plot!(k[end-50:end]./g[end-50], label="control")
    savefig("figs/diff_k.png")
    plot(mock_NLw[end-50:end]./mock_G[end-50], label="mock")
    plot!(NLw[end-50:end]./G[end-50], label="control")
    savefig("figs/diff_NLw.png")
    plot(mock_NLi[end-50:end]./mock_G[end-50], label="mock")
    plot!(NLi[end-50:end]./G[end-50], label="control")
    savefig("figs/diff_NLi.png")
    plot(mock_NLf[end-50:end]./mock_G[end-50], label="mock")
    plot!(NLf[end-50:end]./G[end-50], label="control")
    savefig("figs/diff_NLf.png")
    plot(mock_NLb[end-50:end]./mock_G[end-50], label="mock")
    plot!(NLb[end-50:end]./G[end-50], label="control")
    savefig("figs/diff_NLb.png")
    plot(mock_NLg[end-50:end]./mock_G[end-50], label="mock")
    plot!(NLg[end-50:end]./G[end-50], label="control")
    savefig("figs/diff_NLg.png")
    plot(mock_M[end-50:end]./mock_G[end-50], label="mock")
    plot!(M[end-50:end]./G[end-50], label="control")
    savefig("figs/diff_M.png")
    plot(mock_H[end-50:end]./mock_G[end-50], label="mock")
    plot!(H[end-50:end]./G[end-50], label="control")
    savefig("figs/diff_H.png")
    plot(mock_L[end-50:end]./mock_G[end-50], label="mock")
    plot!(L[end-50:end]./G[end-50], label="control")
    savefig("figs/diff_L.png")
    plot(mock_Te[end-50:end]./mock_G[end-50], label="mock")
    plot!(Te[end-50:end]./G[end-50], label="control")
    savefig("figs/diff_Te.png")
    plot(mock_Tf[end-50:end]./mock_G[end-50], label="mock")
    plot!(Tf[end-50:end]./G[end-50], label="control")
    savefig("figs/diff_Tf.png")
    plot(mock_Ti[end-50:end]./mock_G[end-50], label="mock")
    plot!(Ti[end-50:end]./G[end-50], label="control")
    savefig("figs/diff_Ti.png")
    plot(mock_u[end-50:end], label="mock")
    plot!(u[end-50:end], label="control")
    savefig("figs/diff_u.png")
    plot(mock_Π[end-50:end]./mock_G[end-50], label="mock")
    plot!(Π[end-50:end]./G[end-50], label="control")
    savefig("figs/diff_Π.png")
    plot(mock_NWw[end-50:end]./mock_G[end-50], label="mock")
    plot!(NWw[end-50:end]./G[end-50], label="control")
    savefig("figs/diff_NWw.png")
    plot(mock_NWi[end-50:end]./mock_G[end-50], label="mock")
    plot!(NWi[end-50:end]./G[end-50], label="control")
    savefig("figs/diff_NWi.png")
    plot(mock_NWf[end-50:end]./mock_G[end-50], label="mock")
    plot!(NWf[end-50:end]./G[end-50], label="control")
    savefig("figs/diff_NWf.png")
    plot(mock_NWb[end-50:end]./mock_G[end-50], label="mock")
    plot!(NWb[end-50:end]./G[end-50], label="control")
    savefig("figs/diff_NWb.png")
    plot(mock_NWg[end-50:end]./mock_G[end-50], label="mock")
    plot!(NWg[end-50:end]./G[end-50], label="control")
    savefig("figs/diff_NWg.png")
end

mock_C, mock_c, mock_I, mock_i, mock_p, mock_W, mock_G, mock_g, mock_pe, mock_K, mock_k, mock_NLw, mock_NLi, mock_NLf, mock_NLb, mock_NLg, mock_M, mock_H, mock_L, mock_Te, mock_Tf, mock_Ti, mock_u, mock_Π, mock_NWw, mock_NWi, mock_NWf, mock_NWb, mock_NWg = save_mock_data()
#   初期化
G, GD, g, gD = zeros(T), zeros(T), zeros(T), zeros(T)
Wf, Wg, Wb, W, wg = zeros(T), zeros(T), zeros(T), zeros(T), zeros(T)
Cw, Ci, C, cw, ci, c, CwD, CiD = zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T)
K, k = zeros(T), zeros(T)
Mw, Mi, Mf, M = zeros(T), zeros(T), zeros(T), zeros(T)
Ei, E, Eb, ei, e, eb = zeros(T), zeros(T), zeros(T), zeros(T), fill(100.0, T), zeros(T)
Lw, Lf, L, LfD = zeros(T), zeros(T), zeros(T), zeros(T)
Hw, Hb, H = zeros(T), zeros(T), zeros(T)
NWw, NWi, NWf, NWb, NWg = zeros(T), zeros(T), zeros(T), zeros(T), zeros(T)
NLw, NLi, NLf, NLb, NLg = zeros(T), zeros(T), zeros(T), zeros(T), zeros(T)
Tef, Tei, Tew, Te = zeros(T), zeros(T), zeros(T), zeros(T)
Tii, Tiw, Ti = zeros(T), zeros(T), zeros(T)
Tew, Tei, Tef, Te = zeros(T), zeros(T), zeros(T), zeros(T)
Tff, Tfb, Tf = zeros(T), zeros(T), zeros(T)
Πi, Π, Πf, Πb = zeros(T), zeros(T), zeros(T), zeros(T)
I, i = zeros(T), zeros(T)
p, pe = zeros(T), zeros(T)
u, Δu, ue = zeros(T), zeros(T), zeros(T)
Δk, ΔK = zeros(T), zeros(T), zeros(T), zeros(T)
ΔMw, ΔMi, ΔMf, ΔM = zeros(T), zeros(T), zeros(T), zeros(T)
ΔLw, ΔLf, ΔL = zeros(T), zeros(T), zeros(T)
Δei, Δeb, ΔEi, ΔEb, ΔE = zeros(T), zeros(T), zeros(T), zeros(T), zeros(T)
ΔHw, ΔHb, ΔH = zeros(T), zeros(T), zeros(T)
Δp, Δpe = zeros(T), zeros(T), zeros(T)
Le, Ebe = zeros(T), zeros(T)
We, Tiwe, Tewe = zeros(T), zeros(T), zeros(T)
Πie, Πbe, Ee, Tiie, Teie = zeros(T), zeros(T), zeros(T), zeros(T), zeros(T)
NWie, ΔNWi = zeros(T), zeros(T)
ΔW, ΔTiw, ΔTew = zeros(T), zeros(T), zeros(T)
ΔΠi, ΔΠb, ΔTii, ΔTei = zeros(T), zeros(T), zeros(T), zeros(T)
ΔCi, ΔCw = zeros(T), zeros(T)
EiT, EbT = zeros(T), zeros(T)
K[1] = 300
Mw[1], Mi[1], Mf[1], M[1] = 400.0, 50.0, 150.0, 600.0
Ei[1], E[1], Eb[1] = 50.0, 100.0, 50.0
ei[1], e[1], eb[1] = 50.0, 100.0, 50.0
Lw[1], Lf[1], L[1] = 200.0, 0.0, 200.0
Hw[1], Hb[1], H[1] = 20.0, 580.0, 600.0
G[1], gD[1] = G0, G0
Wf[1], Wg[1], W[1], wg[1] = 1.0*G0, 0.5*G0, 1.5*G0, 0.5*G0
p[1], pe[1] = 1.0, 1.0
k[1] = K[1]/p[1]
CwD[1], CiD[1] = 1.0*G0, 0.1*G0
Cw[1], Ci[1] = CwD[1]/max(1, (CwD[1] + CiD[1])/(p[1]*ζ1*k[1])), CiD[1]/max(1, (CwD[1] + CiD[1])/(p[1]*ζ1*k[1]))
cw[1], ci[1] = Cw[1]/p[1], Ci[1]/p[1]
C[1], c[1] = Cw[1] + Ci[1], cw[1] + ci[1]
u[1] = min(1, (C[1] + Ci[1])/(p[1] * ζ1 * k[1]))
ue[1] = u[1]
Π[1] = C[1] + I[1] + G[1] - Wf[1] - Tef[1] - Tff[1] - r*Lf[1]
Πi[1] = max(0, θ1*(Π[1]-I[1])*Ei[1]/E[1] + θ2*(Mf[1]-Lf[1]))
Πb[1] = max(0, θ1*(Π[1]-I[1])*Eb[1]/E[1] + θ2*(Mf[1]-Lf[1]))
Πf[1] = Π[1] - Πi[1] - Πb[1]
Tiw[1], Tii[1], Tei[1], Tew[1] = τ1*W[1], τ1*Πi[1], γ2*(Mi[1] + Ei[1]), γ2*(Mw[1] + Hw[1])
Le[1] = L[1]
Πie[1], Πbe[1], Ee[1] = Πi[1], Πb[1], E[1]
We[1], Tiwe[1], Tewe[1] = W[1], Tiw[1], Tew[1]
Tiie[1], Teie[1] = Tii[1], Tei[1]
EiT[1], EbT[1] = Ei[1], Eb[1]
NLw[1] = -Cw[1] + W[1] - Tiw[1] - Tew[1] - r*L[1]
ΔMw[1] = NLw[1] # ΔLw[1] = ΔHw[1] = 0ということにする
NLi[1] = Πi[1] - Ci[1] - Tii[1] - Tei[1]  # i[1] = 0を仮定
ΔMi[1] = NLi[1]  # Δei[1] = 0ということにする
NLf[1] = Πf[1]  # i[1] = 0を仮定
ΔMf[1] = NLf[1] # ΔLf[1] = 0を仮定
Ti[1] = Tiw[1] + Tii[1]
Te[1] = Tew[1] + Tei[1] + Tef[1]
ΔM[1] = ΔMw[1] + ΔMi[1] + ΔMf[1]
ΔL[1] = ΔLw[1] + ΔLf[1]
NLb[1] = -Tfb[1] + Πb[1] + r*L[1]
ΔHb[1] = NLb[1] + ΔM[1] - ΔL[1] - pe[1]*Δeb[1]
ΔH[1] = ΔHw[1] + ΔHb[1]
NLg[1] = -G[1] - Wg[1] + Ti[1] + Te[1] + Tf[1]
NWw[1] = Mw[1] - Lw[1] + Hw[1]
NWi[1] = Mi[1] + Ei[1]
NWf[1] = K[1] + Mf[1] - E[1] - Lf[1]
NWb[1] = -M[1] + Eb[1] + L[1] + Hb[1]
NWg[1] = -H[1]
println("------initial------")
println("-G-Wg+Ti+Te+Tf+ΔH=", -G[1]-Wg[1]+Ti[1]+Te[1]+Tf[1]+ΔH[1])
println("NLw+NLi+NLf+NLb+NLg=", NLw[1]+NLi[1]+NLf[1]+NLb[1]+NLg[1])
println("-M+Mw+Mi+Mf=", -M[1]+Mw[1]+Mi[1]+Mf[1])
println("-E+Ei+Eb=", -E[1]+Ei[1]+Eb[1])
println("-L+Lw+Lf=", -L[1]+Lw[1]+Lf[1])
println("-H+Hw+Hf=", -H[1]+Hw[1]+Hb[1])
println("NWw+NWi+NWf+NWb+NWg-K=", NWw[1]+NWi[1]+NWf[1]+NWb[1]+NWg[1]-K[1])
#   パラメータ変更

#   controlのシミュレーションと、mockとcontrolの比較のプロット
run()
plot_mock_and_control()
