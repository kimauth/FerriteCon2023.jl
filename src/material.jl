struct Elasticity
    G::Float64
    K::Float64
end

function material_routine(material::Elasticity, ε::SymmetricTensor{2})
    (; G, K) = material
    stress(ε) = 2G * dev(ε) + K * tr(ε) * one(ε)
    ∂σ∂ε, σ = gradient(stress, ε, :all)
    return σ, ∂σ∂ε
end
