#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# This julia file contains the SYMBOLIC definitions of the geometric quantities as well as the Hessian determinant that enter the amplitude of the effective cosmological model #
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

#----------------------------------------------

function generate_null_sym(v::Array{Num, 1})

    res = [cosh(asinh(-v[1])), cos(atan(v[2],v[3])) + v[1]*sin(atan(v[2],v[3])), -sin(atan(v[2],v[3]))+v[1]*cos(atan(v[2],v[3]))]

    
end

function ϑa_sym(a0::Num, a1::Num, b::Num)

    res = sqrt(((a0 - a1 + sqrt(2(a0 - a1)^2 - 4*b^2))/((a0 - a1 - sqrt(2(a0 - a1)^2 - 4*b^2)))))

end

function ϑb_sym(a0::Num, a1::Num, b::Num)

    res = sqrt(((a0 - a1)^2 - 2*abs(b)*sqrt(2(a0 - a1)^2 - 4*b^2))/((a0 - a1)^2 + 2*abs(b)*sqrt(2(a0 - a1)^2 - 4*b^2)))

end

function n13_refl_sym(a0::Num, a1::Num, b::Num)

    res = [-(2*(a0 - a1)*sqrt(2(a0 - a1)^2 - 4*b^2))/((a0 - a1)^2 - 4*b^2), 0, (-3*(a0 - a1)^2 + 4*b^2)/((a0 - a1)^2 - 4*b^2)]
    
end

function n15_refl_sym(a0::Num, a1::Num, b::Num)

    res = [(2*(a0 - a1)*sqrt(2(a0 - a1)^2 - 4*b^2))/((a0 - a1)^2 - 4*b^2), 0, (3*(a0 - a1)^2 - 4*b^2)/((a0 - a1)^2 - 4*b^2)]
    
end

function n23_refl_sym(a0::Num, a1::Num, b::Num)

    res = [(2*(a0 - a1)*sqrt(2(a0 - a1)^2 - 4*b^2))/((a0 - a1)^2 - 4*b^2), 0, (3*(a0 - a1)^2 - 4*b^2)/((a0 - a1)^2 - 4*b^2)]
    
end

function n25_refl_sym(a0::Num, a1::Num, b::Num)

    res = [-(2*(a0 - a1)*sqrt(2(a0 - a1)^2 - 4*b^2))/((a0 - a1)^2 - 4*b^2), 0, (-3*(a0 - a1)^2 + 4*b^2)/((a0 - a1)^2 - 4*b^2)]
    
end

function n34_refl_sym(a0::Num, a1::Num, b::Num)

    res = [((5*(a0 - a1)^2 - 4*b^2)*sqrt((a0 - a1)^2 -2*b^2))/(sqrt(2)*(-(a0 - a1)^2*b + 4*b^3)), (a0 - a1)/(2*b), (7*(a0 - a1)^3 + 12*(-a0 + a1)*b^2)/(-2*(a0 - a1)^2*b + 8*b^3)]
    
end

function n36_refl_sym(a0::Num, a1::Num, b::Num)

    res = [(sqrt((a0 - a1)^2/2 - b^2))/b, -(a0 - a1)/(2*b), (a0 - a1)/(2*b)]
    
end

function n45_refl_sym(a0::Num, a1::Num, b::Num)

    res = [((5*(a0 - a1)^2 - 4*b^2)*sqrt((a0 - a1)^2 -2*b^2))/(sqrt(2)*(-(a0 - a1)^2*b + 4*b^3)), -(a0 - a1)/(2*b), (7*(a0 - a1)^3 + 12*(-a0 + a1)*b^2)/(-2*(a0 - a1)^2*b + 8*b^3)]

end

function n56_refl_sym(a0::Num, a1::Num, b::Num)

    res = [-(sqrt((a0 - a1)^2/2 - b^2))/b, -(a0 - a1)/(2*b), -(a0 - a1)/(2*b)]
    
end

function edgevectorsets_refl_sym(a0::Num, a1::Num, b::Num)

    res = [ [0, 0, n13_refl_sym(a0, a1, b), n14, n15_refl_sym(a0, a1, b), n16],
            [0, 0, n23_refl_sym(a0, a1, b), n24, n25_refl_sym(a0, a1, b), n26],
            [n13_refl_sym(a0, a1, b), n23, 0, n34_refl_sym(a0, a1, b), 0, n36_refl_sym(a0, a1, b)],
            [n14, n24, n34_refl_sym(a0, a1, b), 0, n45_refl_sym(a0, a1, b), 0],
            [n15_refl_sym(a0, a1, b), n25_refl_sym(a0, a1, b), 0, n45_refl_sym(a0, a1,b), 0, n56_refl_sym(a0, a1, b)],
            [n16, n26, n36_refl_sym(a0, a1, b), 0, n56_refl_sym(a0, a1, b), 0]]
    
end

function nullvectorsets_refl_sym(a0::Num, a1::Num, b::Num)

    res = [ [0, 0, generate_null_sym(n13_refl_sym(a0, a1, b)), generate_null(n14), generate_null_sym(n15_refl_sym(a0, a1, b)), generate_null(n16)],
            [0, 0, generate_null_sym(n23_refl_sym(a0, a1, b)), generate_null(n24), generate_null_sym(n25_refl_sym(a0, a1, b)), generate_null(n26)],
            [generate_null_sym(n13_refl_sym(a0, a1, b)), generate_null(n23), 0, generate_null_sym(n34_refl_sym(a0, a1, b)), 0, generate_null_sym(n36_refl_sym(a0, a1, b))],
            [generate_null(n14), generate_null(n24), generate_null_sym(n34_refl_sym(a0, a1, b)), 0, generate_null_sym(n45_refl_sym(a0, a1, b)), 0],
            [generate_null_sym(n15_refl_sym(a0, a1, b)), generate_null_sym(n25_refl_sym(a0, a1, b)), 0, generate_null_sym(n45_refl_sym(a0, a1,b)), 0, generate_null_sym(n56_refl_sym(a0, a1, b))],
            [generate_null(n16), generate_null(n26), generate_null_sym(n36_refl_sym(a0, a1, b)), 0, generate_null_sym(n56_refl_sym(a0, a1, b)), 0]]
    
end

function lengthsetsSL_sym(a0::Num, a1::Num, b::Num)
    
    res = [[0, 0, a0, a0, a0, a0], [0, 0, a1, a1, a1, a1], [a0, a1, 0, b, 0, b], [a0, a1, b, 0, b, 0], [a0, a1, 0, b, 0, b], [a0, a1, b, 0, b, 0]]

end

function ϑsets_sym(a0::Num, a1::Num, b::Num)

    res = [ [0, 0, ϑa_sym(a0, a1, b)^(-2), ϑa_sym(a0, a1, b)^(-2), ϑa_sym(a0, a1, b)^(-2), ϑa_sym(a0, a1, b)^(-2)],
            [0, 0, ϑa_sym(a0, a1, b)^2, ϑa_sym(a0, a1, b)^2, ϑa_sym(a0, a1, b)^2, ϑa_sym(a0, a1, b)^2],
            [ϑa_sym(a0, a1, b)^(-2), ϑa_sym(a0, a1, b)^2, 0, ϑb_sym(a0, a1, b)^(-2), 0 ,ϑb_sym(a0, a1, b)^(-2)],
            [ϑa_sym(a0, a1, b)^(-2), ϑa_sym(a0, a1, b)^2, ϑb_sym(a0, a1, b)^(-2), 0 ,ϑb_sym(a0, a1, b)^(-2), 0],
            [ϑa_sym(a0, a1, b)^(-2), ϑa_sym(a0, a1, b)^2, 0, ϑb_sym(a0, a1, b)^(-2), 0 ,ϑb_sym(a0, a1, b)^(-2)],
            [ϑa_sym(a0, a1, b)^(-2), ϑa_sym(a0, a1, b)^2, ϑb_sym(a0, a1, b)^(-2), 0 ,ϑb_sym(a0, a1, b)^(-2), 0] ]
            
end

function HessianSLϑ_sym(a0::Num, a1::Num, b0::Num)

    @variables M[1:18, 1:18]::Complex{Real}
    H = [M[i,j] for i in 1:18, j in 1:18]

    for a in 1:6
        for i in 1:3
            for j in 1:3

                A_block = 3(a-1) + i
                
                for b in 1:6

                    B_block = 3(b-1) + j
                    
                    if b in indexsets[a]            
                    
                        H[A_block, B_block] = im/2 * lengthsetsSL_sym(a0, a1, b0)[a][b] * 
                        (η[i][j] + edgevectorsets_refl_sym(a0, a1, b0)[a][b][i]*edgevectorsets_refl_sym(a0, a1, b0)[a][b][j] - 
                        im*ϑsets_sym(a0, a1, b0)[a][b]*
                        nullvectorsets_refl_sym(a0, a1, b0)[a][b][i]*nullvectorsets_refl_sym(a0, a1, b0)[a][b][j] + 
                        ϵsets[a][b]*sum([levicivita([i,j,k])*η[k][l]*edgevectorsets_refl_sym(a0, a1, b0)[a][b][l] for k in 1:3, l in 1:3]))

                    else

                        H[A_block,B_block] = 0

                    end
                end

                C_block = 3(a-1) + j

                H[A_block, C_block] =   -im/2 * sum([ lengthsetsSL_sym(a0, a1, b0)[a][c] * (η[i][j] + edgevectorsets_refl_sym(a0, a1, b0)[a][c][i]*edgevectorsets_refl_sym(a0, a1, b0)[a][c][j] - 
                im*ϑsets_sym(a0, a1, b0)[a][c]*nullvectorsets_refl_sym(a0, a1, b0)[a][c][i]*nullvectorsets_refl_sym(a0, a1, b0)[a][c][j]) for c in indexsets[a]])

            end
        end
    end  
    
    H
    
end