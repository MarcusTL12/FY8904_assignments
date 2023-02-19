using LinearAlgebra
using ForwardDiff
using ReverseDiff

function compute_hamiltonian_anisotropy(S, dz)
    H = 0.0

    @inbounds @simd for i in axes(S, 1)
        H += S[i, 3]^2
    end

    -dz * H
end

function compute_hamiltonian_magnetic_field(S, B)
    H = 0.0

    @inbounds @simd for i in axes(S, 1)
        H += B[1] * S[i, 1] + B[2] * S[i, 2] + B[3] * S[i, 3]
    end

    -H
end

function Sdot(S, x1, y1, z1, x2, y2, z2)
    @inbounds @fastmath begin
        S[x1, y1, z1, 1] * S[x2, y2, z2, 1] +
        S[x1, y1, z1, 2] * S[x2, y2, z2, 2] +
        S[x1, y1, z1, 3] * S[x2, y2, z2, 3]
    end
end

function compute_hamiltonian(S, J, dz, B)
    nx, ny, nz, _ = size(S)
    n = nx * ny * nz

    S_linear = reshape(S, n, 3)

    H = 0.0

    @simd for z in 1:nz-1
        @simd for y in 1:ny-1
            @simd for x in 1:nx-1
                H += @inline Sdot(S, x, y, z, x + 1, y, z)
                H += @inline Sdot(S, x, y, z, x, y + 1, z)
                H += @inline Sdot(S, x, y, z, x, y, z + 1)
            end

            H += @inline Sdot(S, nx, y, z, 1, y, z)
            H += @inline Sdot(S, nx, y, z, nx, y + 1, z)
            H += @inline Sdot(S, nx, y, z, nx, y, z + 1)
        end

        @simd for x in 1:nx-1
            H += @inline Sdot(S, x, ny, z, x + 1, ny, z)
            H += @inline Sdot(S, x, ny, z, x, 1, z)
            H += @inline Sdot(S, x, ny, z, x, ny, z + 1)
        end

        H += @inline Sdot(S, nx, ny, z, 1, ny, z)
        H += @inline Sdot(S, nx, ny, z, nx, 1, z)
        H += @inline Sdot(S, nx, ny, z, nx, ny, z + 1)
    end

    @simd for y in 1:ny-1
        @simd for x in 1:nx-1
            H += @inline Sdot(S, x, y, nz, x + 1, y, nz)
            H += @inline Sdot(S, x, y, nz, x, y + 1, nz)
            H += @inline Sdot(S, x, y, nz, x, y, 1)
        end

        H += @inline Sdot(S, nx, y, nz, 1, y, nz)
        H += @inline Sdot(S, nx, y, nz, nx, y + 1, nz)
        H += @inline Sdot(S, nx, y, nz, nx, y, 1)
    end

    @simd for x in 1:nx-1
        H += @inline Sdot(S, x, ny, nz, x + 1, ny, nz)
        H += @inline Sdot(S, x, ny, nz, x, 1, nz)
        H += @inline Sdot(S, x, ny, nz, x, ny, 1)
    end

    H += @inline Sdot(S, nx, ny, nz, 1, ny, nz)
    H += @inline Sdot(S, nx, ny, nz, nx, 1, nz)
    H += @inline Sdot(S, nx, ny, nz, nx, ny, 1)

    -J * H +
    @inline compute_hamiltonian_anisotropy(S_linear, dz) +
            @inline compute_hamiltonian_magnetic_field(S_linear, B)
end

function make_H_func(J, dz, B)
    S -> compute_hamiltonian(S, J, dz, B)
end

function add_grad_contrib!(∇H, S, J, x1, y1, z1, x2, y2, z2)
    @inbounds @fastmath begin
        ∇H[x1, y1, z1, 1] -= J * S[x2, y2, z2, 1]
        ∇H[x1, y1, z1, 2] -= J * S[x2, y2, z2, 2]
        ∇H[x1, y1, z1, 3] -= J * S[x2, y2, z2, 3]

        ∇H[x2, y2, z2, 1] -= J * S[x1, y1, z1, 1]
        ∇H[x2, y2, z2, 2] -= J * S[x1, y1, z1, 2]
        ∇H[x2, y2, z2, 3] -= J * S[x1, y1, z1, 3]
    end
end

function compute_∇H!(∇H, S, J, dz, B)
    nx, ny, nz, _ = size(S)
    # Magnetic field
    @inbounds @fastmath begin
        fill!((@view ∇H[:, :, :, 1]), -B[1])
        fill!((@view ∇H[:, :, :, 2]), -B[2])
        fill!((@view ∇H[:, :, :, 3]), -B[3])

        # Spin coupling term
        @simd for z in 1:nz-1
            @simd for y in 1:ny-1
                @simd for x in 1:nx-1
                    @inline add_grad_contrib!(∇H, S, J, x, y, z, x + 1, y, z)
                    @inline add_grad_contrib!(∇H, S, J, x, y, z, x, y + 1, z)
                    @inline add_grad_contrib!(∇H, S, J, x, y, z, x, y, z + 1)
                end

                @inline add_grad_contrib!(∇H, S, J, nx, y, z, 1, y, z)
                @inline add_grad_contrib!(∇H, S, J, nx, y, z, nx, y + 1, z)
                @inline add_grad_contrib!(∇H, S, J, nx, y, z, nx, y, z + 1)
            end

            @simd for x in 1:nx-1
                @inline add_grad_contrib!(∇H, S, J, x, ny, z, x + 1, ny, z)
                @inline add_grad_contrib!(∇H, S, J, x, ny, z, x, 1, z)
                @inline add_grad_contrib!(∇H, S, J, x, ny, z, x, ny, z + 1)
            end

            @inline add_grad_contrib!(∇H, S, J, nx, ny, z, 1, ny, z)
            @inline add_grad_contrib!(∇H, S, J, nx, ny, z, nx, 1, z)
            @inline add_grad_contrib!(∇H, S, J, nx, ny, z, nx, ny, z + 1)
        end

        @simd for y in 1:ny-1
            @simd for x in 1:nx-1
                @inline add_grad_contrib!(∇H, S, J, x, y, nz, x + 1, y, nz)
                @inline add_grad_contrib!(∇H, S, J, x, y, nz, x, y + 1, nz)
                @inline add_grad_contrib!(∇H, S, J, x, y, nz, x, y, 1)
            end

            @inline add_grad_contrib!(∇H, S, J, nx, y, nz, 1, y, nz)
            @inline add_grad_contrib!(∇H, S, J, nx, y, nz, nx, y + 1, nz)
            @inline add_grad_contrib!(∇H, S, J, nx, y, nz, nx, y, 1)
        end

        @simd for x in 1:nx-1
            @inline add_grad_contrib!(∇H, S, J, x, ny, nz, x + 1, ny, nz)
            @inline add_grad_contrib!(∇H, S, J, x, ny, nz, x, 1, nz)
            @inline add_grad_contrib!(∇H, S, J, x, ny, nz, x, ny, 1)
        end

        @inline add_grad_contrib!(∇H, S, J, nx, ny, nz, 1, ny, nz)
        @inline add_grad_contrib!(∇H, S, J, nx, ny, nz, nx, 1, nz)
        @inline add_grad_contrib!(∇H, S, J, nx, ny, nz, nx, ny, 1)

        # Anisotropic term
        for z in 1:nz, y in 1:ny, x in 1:nx
            ∇H[x, y, z, 3] -= 2dz * S[x, y, z, 3]
        end
    end
end
