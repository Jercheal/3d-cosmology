function w(a0::Float64, a1::Float64, b::Float64)

    res = (a0 + a1)^2 / (8 * sqrt(complex((a0 - a1)^2/2 + b^2)))

    real(res)

end

function M(a0::Float64, a1::Float64, b::Float64, m::Float64)

    res = m^2/4 * Vol3TL(a0, a1, b)

end

function ϕintegral(a0::Float64, a1::Float64, a2::Float64, b0::Float64, b1::Float64, ϕ0::Float64, ϕ2::Float64, m::Float64)

    res = exp(im*((ϕ0 - ϕ2)^2*((w(a0, a1, b0)*w(a1, a2, b1))/(w(a0, a1, b0) + w(a1, a2, b1) - M(a0, a1, b0, m) - M(a1, a2, b1, m))) 
    - ϕ0^2*(M(a0, a1, b0, m) + (w(a0, a1, b0)*(M(a0, a1, b0, m) + M(a1, a2, b1, m)))/(w(a0, a1, b0) + w(a1, a2, b1) - M(a0, a1, b0, m) - M(a1, a2, b1, m))) 
    - ϕ2^2*(M(a1, a2, b1, m) + (w(a1, a2, b1)*(M(a0, a1, b0, m) + M(a1, a2, b1, m)))/(w(a0, a1, b0) + w(a1, a2, b1) - M(a0, a1, b0, m) - M(a1, a2, b1, m)))))

    res *= sqrt(complex((im * pi)/(w(a0, a1, b0) + w(a1, a2, b1) - M(a0, a1, b0, m) - M(a1, a2, b1, m))))

end

function Ampl_1slice(a0::Float64, a1::Float64, a2::Float64, b0::Float64, b1::Float64, ϕ0::Float64, ϕ2::Float64, m::Float64)

    res = Ampl_vertex_III(a0, a1, b0) * Ampl_vertex_III(a1, a2, b1) * ϕintegral(a0, a1, a2, b0, b1, ϕ0, ϕ2, m)
    
end
