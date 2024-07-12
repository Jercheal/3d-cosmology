η = [[1, 0, 0], [0, -1, 0], [0, 0, -1]]

function Vol3SL(a0::Float64,a1::Float64,b::Float64)

    res = (a0^2 + a0*a1 + a1^2) * sqrt(complex((a0-a1)^2/2 - b^2))/3

    real(res)

end

function Vol3TL(a0::Float64,a1::Float64,b::Float64)

    res = (a0^2 + a0*a1 + a1^2) * sqrt(complex((a0-a1)^2/2 + b^2))/3

    real(res)

end

function generate_null(v::Array{Float64, 1})

    res = [cosh(asinh(-v[1])), cos(atan(v[2],v[3])) + v[1]*sin(atan(v[2],v[3])), -sin(atan(v[2],v[3]))+v[1]*cos(atan(v[2],v[3]))]

    real(res)
    
end

function S_Regge_I(a0::Float64,a1::Float64,b::Float64)

    res = -4*abs(a0-a1)*log(abs(a0-a1) + sqrt(complex2(2*(a0-a1)^2 - 4*b^2))) 
    + 4*abs(b)*log((a0-a1)^2 + sqrt(complex(8*b^2*((a0-a1)^2 - 2*b^2)))) 
    + (4*abs(b) - 2*abs(a0 - a1))*log(4*b^2 - (a0-a1)^2)

    real(res)

end

function S_Regge_null(a0::Float64,a1::Float64,b::Float64)

    res = 4*abs(a0 - a1)*log(abs(a0 - a1)) - 4*abs(b)*log((a0 - a1)^2/2)

    real(res)

end

function S_Regge_II(a0::Float64,a1::Float64,b::Float64)

    res = -4*(a0 - a1)*asinh((a0 - a1)/sqrt(complex((a0 - a1)^2 - 4*b^2))) + 4*abs(b)*acosh((a0 - a1)^2/((a0 - a1)^2 - 4*b^2))

    real(res)

end

function S_Regge_III(a0::Float64,a1::Float64,b::Float64)

    res = -4*(a0 - a1)*asinh((a0 - a1)/sqrt(complex((a0 - a1)^2 - 4*b^2))) + 4*abs(b)*(pi/2 - acos((a0 - a1)^2/(4*b^2 + (a0 - a1)^2)))

    real(res)

end

function S_phi_SL(a0::Float64,a1::Float64,b::Float64,ϕ0::Float64,ϕ1::Float64,m::Float64)

    res = (a0 + a1)^2/(8*sqrt(complex((a0 - a1)^2/2 - b^2)))*(ϕ0 - ϕ1)^2 - m^2/2*Vol3SL(a0,a1,b)*(ϕ0^2 + ϕ1^2)/2

    real(res)

end

function S_phi_TL(a0::Float64,a1::Float64,b::Float64,ϕ0::Float64,ϕ1::Float64,m::Float64)

    res = (a0 + a1)^2/(8*sqrt(complex((a0 - a1)^2/2 - b^2)))*(ϕ0 - ϕ1)^2 + m^2/2*Vol3TL(a0,a1,b)*(ϕ0^2 + ϕ1^2)/2

    real(res)

end

function ϑa(a0::Float64,a1::Float64,b::Float64)

    res = sqrt(complex((a0 - a1 + sqrt(complex(2(a0 - a1)^2 - 4*b^2)))/(a0 - a1 - sqrt(complex(2(a0 - a1)^2 - 4*b^2)))))

    real(res)

end

function ϑb(a0::Float64,a1::Float64,b::Float64)

    res = sqrt(complex(((a0 - a1)^2 - 2*abs(b)*sqrt(complex(2(a0 - a1)^2 - 4*b^2)))/((a0 - a1)^2 + 2*abs(b)*sqrt(complex(2(a0 - a1)^2 - 4*b^2)))))

    real(res)

end

n13 = float([0, 0, 1])
n14 = float([0, -1, 0])
n15 = float([0, 0, -1])
n16 = float([0, 1, 0])
n23 = float([0, 0, -1])
n24 = float([0, 1, 0])
n25 = float([0, 0, 1])
n26 = float([0, -1, 0])

function n13_refl(a0::Float64, a1::Float64, b::Float64)

    res = [-(2*(a0 - a1)*sqrt(complex(2(a0 - a1)^2 - 4*b^2)))/((a0 - a1)^2 - 4*b^2), 0, (-3*(a0 - a1)^2 + 4*b^2)/((a0 - a1)^2 - 4*b^2)]

    real(res)
    
end

function n15_refl(a0::Float64, a1::Float64, b::Float64)

    res = [(2*(a0 - a1)*sqrt(complex(2(a0 - a1)^2 - 4*b^2)))/((a0 - a1)^2 - 4*b^2), 0, (3*(a0 - a1)^2 - 4*b^2)/((a0 - a1)^2 - 4*b^2)]

    real(res)
    
end

function n23_refl(a0::Float64, a1::Float64, b::Float64)

    res = [(2*(a0 - a1)*sqrt(complex(2(a0 - a1)^2 - 4*b^2)))/((a0 - a1)^2 - 4*b^2), 0, (3*(a0 - a1)^2 - 4*b^2)/((a0 - a1)^2 - 4*b^2)]

    real(res)
    
end

function n25_refl(a0::Float64, a1::Float64, b::Float64)

    res = [-(2*(a0 - a1)*sqrt(complex(2(a0 - a1)^2 - 4*b^2)))/((a0 - a1)^2 - 4*b^2), 0, (-3*(a0 - a1)^2 + 4*b^2)/((a0 - a1)^2 - 4*b^2)]

    real(res)
    
end

function n34_refl(a0::Float64, a1::Float64, b::Float64)

    res = [((5*(a0 - a1)^2 - 4*b^2)*sqrt(complex((a0 - a1)^2 -2*b^2)))/(sqrt(2)*(-(a0 - a1)^2*b + 4*b^3)), (a0 - a1)/(2*b), (7*(a0 - a1)^3 + 12*(-a0 + a1)*b^2)/(-2*(a0 - a1)^2*b + 8*b^3)]

    real(res)
    
end

function n36_refl(a0::Float64, a1::Float64, b::Float64)

    res = [(sqrt(complex((a0 - a1)^2/2 - b^2)))/b, -(a0 - a1)/(2*b), (a0 - a1)/(2*b)]

    real(res)
    
end

function n45_refl(a0::Float64, a1::Float64, b::Float64)

    res = [((5*(a0 - a1)^2 - 4*b^2)*sqrt(complex((a0 - a1)^2 -2*b^2)))/(sqrt(2)*(-(a0 - a1)^2*b + 4*b^3)), -(a0 - a1)/(2*b), (7*(a0 - a1)^3 + 12*(-a0 + a1)*b^2)/(-2*(a0 - a1)^2*b + 8*b^3)]
    
    real(res)

end

function n56_refl(a0::Float64, a1::Float64, b::Float64)

    res = [-(sqrt(complex((a0 - a1)^2/2 - b^2)))/b, -(a0 - a1)/(2*b), -(a0 - a1)/(2*b)]

    real(res)
    
end

function edgevectorsets_refl(a0::Float64, a1::Float64, b::Float64)

    res = [ [0, 0, n13_refl(a0, a1, b), n14, n15_refl(a0, a1, b), n16],
            [0, 0, n23_refl(a0, a1, b), n24, n25_refl(a0, a1, b), n26],
            [n13_refl(a0, a1, b), n23, 0, n34_refl(a0, a1, b), 0, n36_refl(a0, a1, b)],
            [n14, n24, n34_refl(a0, a1, b), 0, n45_refl(a0, a1, b), 0],
            [n15_refl(a0, a1, b), n25_refl(a0, a1, b), 0, n45_refl(a0, a1,b), 0, n56_refl(a0, a1, b)],
            [n16, n26, n36_refl(a0, a1, b), 0, n56_refl(a0, a1, b), 0]]

            real(res)
    
end

function nullvectorsets_refl(a0::Float64, a1::Float64, b::Float64)

    res = [ [0, 0, generate_null(n13_refl(a0, a1, b)), generate_null(n14), generate_null(n15_refl(a0, a1, b)), generate_null(n16)],
            [0, 0, generate_null(n23_refl(a0, a1, b)), generate_null(n24), generate_null(n25_refl(a0, a1, b)), generate_null(n26)],
            [generate_null(n13_refl(a0, a1, b)), generate_null(n23), 0, generate_null(n34_refl(a0, a1, b)), 0, generate_null(n36_refl(a0, a1, b))],
            [generate_null(n14), generate_null(n24), generate_null(n34_refl(a0, a1, b)), 0, generate_null(n45_refl(a0, a1, b)), 0],
            [generate_null(n15_refl(a0, a1, b)), generate_null(n25_refl(a0, a1, b)), 0, generate_null(n45_refl(a0, a1,b)), 0, generate_null(n56_refl(a0, a1, b))],
            [generate_null(n16), generate_null(n26), generate_null(n36_refl(a0, a1, b)), 0, generate_null(n56_refl(a0, a1, b)), 0]]

    real(res)
    
end

indexsets = [[3, 4, 5, 6], [3, 4, 5, 6], [1, 2, 4, 6], [1, 2, 3, 5], [1, 2, 4, 6], [1, 2, 3, 5]]

function lengthsetsSL(a0::Float64, a1::Float64, b::Float64)
    
    res = [[0, 0, a0, a0, a0, a0], [0, 0, a1, a1, a1, a1], [a0, a1, 0, b, 0, b], [a0, a1, b, 0, b, 0], [a0, a1, 0, b, 0, b], [a0, a1, b, 0, b, 0]]

    real(res)

end

ϵsets = [float([0, 0, 1, 1, 1, 1]), float([0, 0, 1, 1, 1, 1]), float([-1, -1, 0, 1, 0, 1]), float([-1, -1, -1, 0, 1, 0]), float([-1, -1, 0, -1, 0, 1]), float([-1, -1, -1, 0, -1, 0])]

function ϑsets(a0::Float64, a1::Float64, b::Float64)

    res = [ [0, 0, ϑa(a0, a1, b)^(-2), ϑa(a0, a1, b)^(-2), ϑa(a0, a1, b)^(-2), ϑa(a0, a1, b)^(-2)],
            [0, 0, ϑa(a0, a1, b)^2, ϑa(a0, a1, b)^2, ϑa(a0, a1, b)^2, ϑa(a0, a1, b)^2],
            [ϑa(a0, a1, b)^(-2), ϑa(a0, a1, b)^2, 0, ϑb(a0, a1, b)^(-2), 0 ,ϑb(a0, a1, b)^(-2)],
            [ϑa(a0, a1, b)^(-2), ϑa(a0, a1, b)^2, ϑb(a0, a1, b)^(-2), 0 ,ϑb(a0, a1, b)^(-2), 0],
            [ϑa(a0, a1, b)^(-2), ϑa(a0, a1, b)^2, 0, ϑb(a0, a1, b)^(-2), 0 ,ϑb(a0, a1, b)^(-2)],
            [ϑa(a0, a1, b)^(-2), ϑa(a0, a1, b)^2, ϑb(a0, a1, b)^(-2), 0 ,ϑb(a0, a1, b)^(-2), 0]]

    real(res)

end

function HessianSLϑ(a0::Float64, a1::Float64, b0::Float64)

    H = zeros(ComplexF64,18,18)

    for a in 1:6
        for i in 1:3
            for j in 1:3

                A_block = 3(a-1) + i
                C_block = 3(a-1) + j

                H[A_block, C_block] =   -im/2 * sum([ lengthsetsSL(a0, a1, b0)[a][c] * (η[i][j] + edgevectorsets_refl(a0, a1, b0)[a][c][i]*edgevectorsets_refl(a0, a1, b0)[a][c][j] - 
                im*ϑsets(a0, a1, b0)[a][c]*nullvectorsets_refl(a0, a1, b0)[a][c][i]*nullvectorsets_refl(a0, a1, b0)[a][c][j]) for c in indexsets[a]])

                for b in indexsets[a]

                    B_block = 3(b-1) + j             
                    
                    H[A_block, B_block] = im/2 * lengthsetsSL(a0, a1, b0)[a][b] * 
                    (η[i][j] + edgevectorsets_refl(a0, a1, b0)[a][b][i]*edgevectorsets_refl(a0, a1, b0)[a][b][j] - 
                    im*ϑsets(a0, a1, b0)[a][b]*
                    nullvectorsets_refl(a0, a1, b0)[a][b][i]*nullvectorsets_refl(a0, a1, b0)[a][b][j] + 
                    ϵsets[a][b]*sum([levicivita([i,j,k])*η[k][l]*edgevectorsets_refl(a0, a1, b0)[a][b][l] for k in 1:3, l in 1:3]))

                end
            end
        end
    end  
    
    H
    
end

function DetSLId(a0::Float64, a1::Float64, b::Float64)

    res = (1/b^4)*((1/32 + im/32)*a0^3*a1^3*sqrt(complex(2*(a0 - a1)^2 - 4*b^2))*((4 + 7*im)*a0^12 - (44 + 56*im)*a0^11*a1 + (208 + 214*im)*a0^10*a1^2 - (572 + 536*im)*a0^9*a1^3 + (1052 + 1001*im)*a0^8*a1^4 - (1432 + 1456*im)*a0^7*a1^5 + (1568 + 1652*im)*a0^6*a1^6 - (1432 + 1456*im)*a0^5*a1^7 + (1052 + 1001*im)*a0^4*a1^8 - (572 + 536*im)*a0^3*a1^9 + (208 + 214*im)*a0^2*a1^10 - (44 + 56*im)*a0*a1^11 + (4 + 7*im)*a1^12 + (16 + 24*im)*a0^11*b - (128 + 148*im)*a0^10*a1*b + (448 + 380*im)*a0^9*a1^2*b - (880 + 520*im)*a0^8*a1^3*b + (992 + 400*im)*a0^7*a1^4*b - (448 + 136*im)*a0^6*a1^5*b - (448 + 136*im)*a0^5*a1^6*b + (992 + 400*im)*a0^4*a1^7*b - (880 + 520*im)*a0^3*a1^8*b + (448 + 380*im)*a0^2*a1^9*b - (128 + 148*im)*a0*a1^10*b + (16 + 24*im)*a1^11*b + (18 - 16*im)*a0^10*b^2 + (10 + 126*im)*a0^9*a1*b^2 - (458 + 492*im)*a0^8*a1^2*b^2 + (1648 + 1232*im)*a0^7*a1^3*b^2 - (3080 + 2116*im)*a0^6*a1^4*b^2 + (3724 + 2532*im)*a0^5*a1^5*b^2 - (3080 + 2116*im)*a0^4*a1^6*b^2 + (1648 + 1232*im)*a0^3*a1^7*b^2 - (458 + 492*im)*a0^2*a1^8*b^2 + (10 + 126*im)*a0*a1^9*b^2 + (18 - 16*im)*a1^10*b^2 + (16 - 72*im)*a0^9*b^3 + (44 + 232*im)*a0^8*a1*b^3 - (444 + 176*im)*a0^7*a1^2*b^3 + (908 - 144*im)*a0^6*a1^3*b^3 - (524 - 160*im)*a0^5*a1^4*b^3 - (524 - 160*im)*a0^4*a1^5*b^3 + (908 - 144*im)*a0^3*a1^6*b^3 - (444 + 176*im)*a0^2*a1^7*b^3 + (44 + 232*im)*a0*a1^8*b^3 + (16 - 72*im)*a1^9*b^3 - 8*a0^8*b^4 - (308 + 172*im)*a0^7*a1*b^4 + (1272 + 744*im)*a0^6*a1^2*b^4 - (2188 + 1428*im)*a0^5*a1^3*b^4 + (2464 + 1712*im)*a0^4*a1^4*b^4 - (2188 + 1428*im)*a0^3*a1^5*b^4 + (1272 + 744*im)*a0^2*a1^6*b^4 - (308 + 172*im)*a0*a1^7*b^4 - 8*a1^8*b^4 - (104 - 8*im)*a0^7*b^5 + (96 + 256*im)*a0^6*a1*b^5 + (272 - 592*im)*a0^5*a1^2*b^5 - (264 - 328*im)*a0^4*a1^3*b^5 - (264 - 328*im)*a0^3*a1^4*b^5 + (272 - 592*im)*a0^2*a1^5*b^5 + (96 + 256*im)*a0*a1^6*b^5 - (104 - 8*im)*a1^7*b^5 - (32 - 48*im)*a0^6*b^6 + (456 + 216*im)*a0^5*a1*b^6 - (512 + 528*im)*a0^4*a1^2*b^6 + (176 + 528*im)*a0^3*a1^3*b^6 - (512 + 528*im)*a0^2*a1^4*b^6 + (456 + 216*im)*a0*a1^5*b^6 - (32 - 48*im)*a1^6*b^6 + (96 + 96*im)*a0^5*b^7 + (64 - 256*im)*a0^4*a1*b^7 + (32 - 32*im)*a0^3*a1^2*b^7 + (32 - 32*im)*a0^2*a1^3*b^7 + (64 - 256*im)*a0*a1^4*b^7 + (96 + 96*im)*a1^5*b^7 - 64*im*a0^4*b^8 - (32 + 96*im)*a0^3*a1*b^8 - (64 + 64*im)*a0^2*a1^2*b^8 - (32 + 96*im)*a0*a1^3*b^8 - 64*im*a1^4*b^8)) + (1/b^4)*((1/32 + im/32)*a0^3*a1^3*((8 + 7*im)*a0^13 - (72 + 57*im)*a0^12*a1 + (288 + 198*im)*a0^11*a1^2 - (672 + 378*im)*a0^10*a1^3 + (1000 + 425*im)*a0^9*a1^4 - (936 + 279*im)*a0^8*a1^5 + (384 + 84*im)*a0^7*a1^6 + (384 + 84*im)*a0^6*a1^7 - (936 + 279*im)*a0^5*a1^8 + (1000 + 425*im)*a0^4*a1^9 - (672 + 378*im)*a0^3*a1^10 + (288 + 198*im)*a0^2*a1^11 - (72 + 57*im)*a0*a1^12 + (8 + 7*im)*a1^13 + (12 + 28*im)*a0^12*b - (104 + 252*im)*a0^11*a1*b + (360 + 1056*im)*a0^10*a1^2*b - (584 + 2764*im)*a0^9*a1^3*b + (244 + 5124*im)*a0^8*a1^4*b + (688 - 7224*im)*a0^7*a1^5*b - (1232 - 8064*im)*a0^6*a1^6*b + (688 - 7224*im)*a0^5*a1^7*b + (244 + 5124*im)*a0^4*a1^8*b - (584 + 2764*im)*a0^3*a1^9*b + (360 + 1056*im)*a0^2*a1^10*b - (104 + 252*im)*a0*a1^11*b + (12 + 28*im)*a1^12*b + (2 + 2*im)*a0^11*b^2 - (68 + 52*im)*a0^10*a1*b^2 + (384 + 344*im)*a0^9*a1^2*b^2 - (970 + 1010*im)*a0^8*a1^3*b^2 + (1256 + 1456*im)*a0^7*a1^4*b^2 - (604 + 740*im)*a0^6*a1^5*b^2 - (604 + 740*im)*a0^5*a1^6*b^2 + (1256 + 1456*im)*a0^4*a1^7*b^2 - (970 + 1010*im)*a0^3*a1^8*b^2 + (384 + 344*im)*a0^2*a1^9*b^2 - (68 + 52*im)*a0*a1^10*b^2 + (2 + 2*im)*a1^11*b^2 + (60 - 108*im)*a0^10*b^3 - (332 - 944*im)*a0^9*a1*b^3 + (908 - 3548*im)*a0^8*a1^2*b^3 - (1808 - 7808*im)*a0^7*a1^3*b^3 + (2872 - 11704*im)*a0^6*a1^4*b^3 - (3400 - 13216*im)*a0^5*a1^5*b^3 + (2872 - 11704*im)*a0^4*a1^6*b^3 - (1808 - 7808*im)*a0^3*a1^7*b^3 + (908 - 3548*im)*a0^2*a1^8*b^3 - (332 - 944*im)*a0*a1^9*b^3 + (60 - 108*im)*a1^10*b^3 - (12 + 56*im)*a0^9*b^4 + (504 + 428*im)*a0^8*a1*b^4 - (2028 + 1420*im)*a0^7*a1^2*b^4 + (3180 + 2252*im)*a0^6*a1^3*b^4 - (1644 + 1204*im)*a0^5*a1^4*b^4 - (1644 + 1204*im)*a0^4*a1^5*b^4 + (3180 + 2252*im)*a0^3*a1^6*b^4 - (2028 + 1420*im)*a0^2*a1^7*b^4 + (504 + 428*im)*a0*a1^8*b^4 - (12 + 56*im)*a1^9*b^4 - (200 - 88*im)*a0^8*b^5 + (1024 - 1272*im)*a0^7*a1*b^5 - (2208 - 4016*im)*a0^6*a1^2*b^5 + (2816 - 6024*im)*a0^5*a1^3*b^5 - (2864 - 6384*im)*a0^4*a1^4*b^5 + (2816 - 6024*im)*a0^3*a1^5*b^5 - (2208 - 4016*im)*a0^2*a1^6*b^5 + (1024 - 1272*im)*a0*a1^7*b^5 - (200 - 88*im)*a1^8*b^5 + (16 - 48*im)*a0^7*b^6 - (944 + 208*im)*a0^6*a1*b^6 + (2064 + 816*im)*a0^5*a1^2*b^6 - (1136 + 560*im)*a0^4*a1^3*b^6 - (1136 + 560*im)*a0^3*a1^4*b^6 + (2064 + 816*im)*a0^2*a1^5*b^6 - (944 + 208*im)*a0*a1^6*b^6 + (16 - 48*im)*a1^7*b^6 + 32*im*a0^6*b^7 - (592 - 848*im)*a0^5*a1*b^7 + (960 - 1120*im)*a0^4*a1^2*b^7 - (736 - 480*im)*a0^3*a1^3*b^7 + (960 - 1120*im)*a0^2*a1^4*b^7 - (592 - 848*im)*a0*a1^5*b^7 + 32*im*a1^6*b^7 - (128 - 192*im)*a0^5*b^8 + (576 + 128*im)*a0^4*a1*b^8 - (64 - 64*im)*a0^3*a1^2*b^8 - (64 - 64*im)*a0^2*a1^3*b^8 + (576 + 128*im)*a0*a1^4*b^8 - (128 - 192*im)*a1^5*b^8 + 128*a0^4*b^9 + (192 - 64*im)*a0^3*a1*b^9 + (128 - 128*im)*a0^2*a1^2*b^9 + (192 - 64*im)*a0*a1^3*b^9 + 128*a1^4*b^9))
    
    res

end