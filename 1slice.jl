#--------------------------------------------------------------------------------------------------------------#
# This julia file contains the definition of the 3d effective cosmological spin foam model with one bulk slice #
#--------------------------------------------------------------------------------------------------------------#

#----------------------------------------------

### Definition of weights entering the scalar field action

function w(a0::Float64, a1::Float64, b::Float64)

    res = (a0 + a1)^2 / (8 * sqrt(complex((a0 - a1)^2/2 + b^2)))

    real(res)

end

function w_big(a0::BigFloat, a1::BigFloat, b::BigFloat)

    res = (a0 + a1)^2 / (8 * sqrt(complex((a0 - a1)^2/2 + b^2)))

    real(res)

end

function M(a0::Float64, a1::Float64, b::Float64, m::Float64)

    res = m^2/4 * Vol3TL(a0, a1, b)

end

function M_big(a0::BigFloat, a1::BigFloat, b::BigFloat, m::BigFloat)

    res = m^2/4 * Vol3TL_big(a0, a1, b)

end

#----------------------------------------------

### Scalar field solutions for single bulk slice

function ϕ1_sol(a0::Float64, a1::Float64, a2::Float64, b0::Float64, b1::Float64, ϕ0::Float64, ϕ2::Float64, m::Float64)

    res = (ϕ0 * w(a0, a1, b0) + ϕ2 * w(a1, a2, b1))/(w(a0, a1, b0) + w(a1, a2, b1) - M(a0, a1, b0, m) + M(a1, a2, b1, m))

end

#----------------------------------------------

### Scalar field amplitude after integrating out the bulk scalar field 

function S_phi_eff(a0::Float64, a1::Float64, a2::Float64, b0::Float64, b1::Float64, ϕ0::Float64, ϕ2::Float64, m::Float64)

    (ϕ0 - ϕ2)^2*((w(a0, a1, b0)*w(a1, a2, b1))/(w(a0, a1, b0) + w(a1, a2, b1) - M(a0, a1, b0, m) - M(a1, a2, b1, m))) - ϕ0^2*(M(a0, a1, b0, m) + (w(a0, a1, b0)*(M(a0, a1, b0, m) + M(a1, a2, b1, m)))/(w(a0, a1, b0) + w(a1, a2, b1) - M(a0, a1, b0, m) - M(a1, a2, b1, m))) - ϕ2^2*(M(a1, a2, b1, m) + (w(a1, a2, b1)*(M(a0, a1, b0, m) + M(a1, a2, b1, m)))/(w(a0, a1, b0) + w(a1, a2, b1) - M(a0, a1, b0, m) - M(a1, a2, b1, m)))
    
end

function S_phi_eff_big(a0::BigFloat, a1::BigFloat, a2::BigFloat, b0::BigFloat, b1::BigFloat, ϕ0::BigFloat, ϕ2::BigFloat, m::BigFloat)

    (ϕ0 - ϕ2)^2*((w_big(a0, a1, b0)*w_big(a1, a2, b1))/(w_big(a0, a1, b0) + w_big(a1, a2, b1) - M_big(a0, a1, b0, m) - M_big(a1, a2, b1, m))) - ϕ0^2*(M_big(a0, a1, b0, m) + (w_big(a0, a1, b0)*(M_big(a0, a1, b0, m) + M_big(a1, a2, b1, m)))/(w_big(a0, a1, b0) + w_big(a1, a2, b1) - M_big(a0, a1, b0, m) - M_big(a1, a2, b1, m))) - ϕ2^2*(M_big(a1, a2, b1, m) + (w_big(a1, a2, b1)*(M_big(a0, a1, b0, m) + M_big(a1, a2, b1, m)))/(w_big(a0, a1, b0) + w_big(a1, a2, b1) - M_big(a0, a1, b0, m) - M_big(a1, a2, b1, m)))

end

function ϕintegral(a0::Float64, a1::Float64, a2::Float64, b0::Float64, b1::Float64, ϕ0::Float64, ϕ2::Float64, m::Float64)

    res = exp(im*((ϕ0 - ϕ2)^2*((w(a0, a1, b0)*w(a1, a2, b1))/(w(a0, a1, b0) + w(a1, a2, b1) - M(a0, a1, b0, m) - M(a1, a2, b1, m))) - ϕ0^2*(M(a0, a1, b0, m) + (w(a0, a1, b0)*(M(a0, a1, b0, m) + M(a1, a2, b1, m)))/(w(a0, a1, b0) + w(a1, a2, b1) - M(a0, a1, b0, m) - M(a1, a2, b1, m))) - ϕ2^2*(M(a1, a2, b1, m) + (w(a1, a2, b1)*(M(a0, a1, b0, m) + M(a1, a2, b1, m)))/(w(a0, a1, b0) + w(a1, a2, b1) - M(a0, a1, b0, m) - M(a1, a2, b1, m)))))

    res *= sqrt(complex((im * pi)/(w(a0, a1, b0) + w(a1, a2, b1) - M(a0, a1, b0, m) - M(a1, a2, b1, m))))

end

function ϕintegral_big(a0::BigFloat, a1::BigFloat, a2::BigFloat, b0::BigFloat, b1::BigFloat, ϕ0::BigFloat, ϕ2::BigFloat, m::BigFloat)

    res = exp(im*((ϕ0 - ϕ2)^2*((w_big(a0, a1, b0)*w_big(a1, a2, b1))/(w_big(a0, a1, b0) + w_big(a1, a2, b1) - M_big(a0, a1, b0, m) - M_big(a1, a2, b1, m))) - ϕ0^2*(M_big(a0, a1, b0, m) + (w_big(a0, a1, b0)*(M_big(a0, a1, b0, m) + M_big(a1, a2, b1, m)))/(w_big(a0, a1, b0) + w_big(a1, a2, b1) - M_big(a0, a1, b0, m) - M_big(a1, a2, b1, m))) - ϕ2^2*(M_big(a1, a2, b1, m) + (w_big(a1, a2, b1)*(M_big(a0, a1, b0, m) + M_big(a1, a2, b1, m)))/(w_big(a0, a1, b0) + w_big(a1, a2, b1) - M_big(a0, a1, b0, m) - M_big(a1, a2, b1, m)))))

    res *= sqrt(complex((im * pi)/(w_big(a0, a1, b0) + w_big(a1, a2, b1) - M_big(a0, a1, b0, m) - M_big(a1, a2, b1, m))))

end

#----------------------------------------------

### Full amplitude in Sector III before integrating out the scalar field

function Ampl_1slice_phi1sol(a0::Float64, a1::Float64, a2::Float64, b0::Float64, b1::Float64, ϕ0::Float64, ϕ1sol::Float64, ϕ2::Float64, m::Float64)

    res = Ampl_face_SL(a1) * Ampl_vertex_III(a0, a1, b0) * Ampl_vertex_III(a1, a2, b1) * exp(im * (S_phi_TL(a0, a1, b0, ϕ0, ϕ1sol, m) + S_phi_TL(a1, a2, b1, ϕ1sol, ϕ2, m)))
    
end

#----------------------------------------------

### Full amplitude in Sector III after integrating out the scalar field

function Ampl_1slice(a0::Float64, a1::Float64, a2::Float64, b0::Float64, b1::Float64, ϕ0::Float64, ϕ2::Float64, m::Float64)

    res = Ampl_face_SL(a1) * Ampl_vertex_III(a0, a1, b0) * Ampl_vertex_III(a1, a2, b1) * ϕintegral(a0, a1, a2, b0, b1, ϕ0, ϕ2, m)
    
end

#----------------------------------------------

### Full ESF amplitude in Sector III after integrating out the scalar field

function Ampl_1slice_ESF(a0::Float64, a1::Float64, a2::Float64, b0::Float64, b1::Float64, ϕ0::Float64, ϕ2::Float64, m::Float64)

    res = Ampl_face_SL(a1) * Ampl_vertex_III_ESF(a0, a1, b0) * Ampl_vertex_III_ESF(a1, a2, b1) * ϕintegral(a0, a1, a2, b0, b1, ϕ0, ϕ2, m)
    
end

function Ampl_1slice_ESF_big(a0::BigFloat, a1::BigFloat, a2::BigFloat, b0::BigFloat, b1::BigFloat, ϕ0::BigFloat, ϕ2::BigFloat, m::BigFloat)

    res = Ampl_vertex_III_ESF_big(a0, a1, b0) * Ampl_vertex_III_ESF_big(a1, a2, b1) * ϕintegral_big(a0, a1, a2, b0, b1, ϕ0, ϕ2, m)
    
end

#----------------------------------------------

### Toy amplitude with modified ESF-measure and no phi-measure

function Ampl_1slice_toy(a0::Float64, a1::Float64, a2::Float64, b0::Float64, b1::Float64, ϕ0::Float64, ϕ2::Float64, m::Float64)

    res = Ampl_face_SL(a1) * Ampl_face_TL(b0) * μcont_TL(a0, a1, b0)^2 * exp(im * S_Regge_III(a0, a1, b0)) * Ampl_face_TL(b1) * μcont_TL(a1, a2, b1)^2 * exp(im * S_Regge_III(a1, a2, b1)) * exp(im * S_phi_eff(a0, a1, a2, b0, b1, ϕ0, ϕ2, m))
    
end

#----------------------------------------------

### Toy amplitude with modified ESF-measure and phi-measure

function Ampl_1slice_toy_phi_measure(a0::Float64, a1::Float64, a2::Float64, b0::Float64, b1::Float64, ϕ0::Float64, ϕ2::Float64, m::Float64)

    res = Ampl_face_SL(a1) * Ampl_face_TL(b0) * μcont_TL(a0, a1, b0)^2 * exp(im * S_Regge_III(a0, a1, b0)) * Ampl_face_TL(b1) * μcont_TL(a1, a2, b1)^2 * exp(im * S_Regge_III(a1, a2, b1)) *  ϕintegral(a0, a1, a2, b0, b1, ϕ0, ϕ2, m)
    
end

#----------------------------------------------

### Toy measure and amplitude

function μ_toy_1slice(a0::Float64, a1::Float64, a2::Float64, b0::Float64, b1::Float64)

    (b0*a1)/((a0-a1)^2/2+b0^2)^2 * (b1*a1)/((a1-a2)^2/2+b1^2)^2

end

function Ampl_1slice_toy_working(a0::Float64, a1::Float64, a2::Float64, b0::Float64, b1::Float64, ϕ0::Float64, ϕ2::Float64, m::Float64)

    res =  μ_toy_1slice(a0, a1, a2, b0, b1) * exp(im * S_Regge_III(a0, a1, b0)) * exp(im * S_Regge_III(a1, a2, b1)) * ϕintegral(a0, a1, a2, b0, b1, ϕ0, ϕ2, m)
    
end