# compressible linear elasticity
function assemble_cell!(ke, cellvalues::CellValues, material, xe)
    fill!(ke, 0.0)

    reinit!(cellvalues, xe)

    n_basefuncs = getnbasefunctions(cellvalues)
    for q_point in 1:getnquadpoints(cellvalues)
        ## For each integration point, compute strain, stress and material stiffness
        _ε = shape_symmetric_gradient(cellvalues, q_point, 1)
        _, ∂σ∂ε = material_routine(material, _ε) # _ε needed for dimension of stiffness

        dΩ = getdetJdV(cellvalues, q_point)
        for i in 1:n_basefuncs
            ∇Nᵢ = shape_gradient(cellvalues, q_point, i)
            for j in 1:n_basefuncs
                ∇ˢʸᵐNⱼ = shape_symmetric_gradient(cellvalues, q_point, j)
                ke[i, j] += (∂σ∂ε ⊡ ∇ˢʸᵐNⱼ) ⊡ ∇Nᵢ * dΩ
            end
        end
    end
end

# incompressible elasticity
function assemble_cell!(ke, values::Tuple{<:CellValues, <:CellValues}, material, xe)
    fill!(ke, 0.0)

    (cellvalues_u, cellvalues_p) = values
    (; G, K) = material

    n_basefuncs_u = getnbasefunctions(cellvalues_u)
    n_basefuncs_p = getnbasefunctions(cellvalues_p)

    u▄, p▄ = 1, 2

    reinit!(cellvalues_u, xe)
    reinit!(cellvalues_p, xe)

    for q_point in 1:getnquadpoints(cellvalues_u)
        dΩ = getdetJdV(cellvalues_u, q_point)
        for i in 1:n_basefuncs_u
            δε = shape_symmetric_gradient(cellvalues_u, q_point, i)
            divδu = shape_divergence(cellvalues_u, q_point, i)
            δu = shape_value(cellvalues_u, q_point, i)
            for j in 1:n_basefuncs_u
                εdev = dev(shape_symmetric_gradient(cellvalues_u, q_point, j))
                ke[BlockIndex((u▄, u▄), (i, j))] += 2G * ɛdev ⊡ δɛ * dΩ
            end
            for j in 1:n_basefuncs_p
                p = shape_value(cellvalues_p, q_point, j)
                ke[BlockIndex((u▄, p▄), (i, j))] += -p * divδu * dΩ
            end
        end

        for i in 1:n_basefuncs_p
            δp = shape_value(cellvalues_p, q_point, i)
            for j in 1:n_basefuncs_u
                divδu = shape_divergence(cellvalues_u, q_point, j)
                ke[BlockIndex((p▄, u▄), (i, j))] += -δp * divδu * dΩ
            end
            for j in 1:n_basefuncs_p
                p = shape_value(cellvalues_p, q_point, j)
                ke[BlockIndex((p▄, p▄), (i, j))] += - p * δp / K * dΩ
            end

        end
    end
end
