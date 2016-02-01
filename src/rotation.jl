

function rotmat(u,θ)
    ux,uy,uz = u
    cos(θ)*eye(3) + sin(θ)*[[0 -uz uy];[uz 0 -ux];[-uy ux 0]] + (1-cos(θ))*
        [[ux^2 ux*uy ux*uz];[ux*uy uy^2 uy*uz];[ux*uz uy*uz uz^2]]
end

function randrotmat()
    u = rand(3)
    u ./= norm(u)
    θ = 2π*rand()
    rotmat(u,θ)
end

function rotate!(Ω::AngularGrid, R::Matrix)
    for i in 1:length(Ω)
        Ω[i] = R * Ω[i]
    end
end





