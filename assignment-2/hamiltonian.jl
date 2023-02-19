using LinearAlgebra
using ForwardDiff
using SIMD

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

function add_spin_coupling_term!(∇H, S, J, nx, ny, nz)
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
end

function compute_∇H!(∇H, S, J, dz, B)
    nx, ny, nz, _ = size(S)
    # Magnetic field
    @inbounds @fastmath begin
        fill!((@view ∇H[:, :, :, 1]), -B[1])
        fill!((@view ∇H[:, :, :, 2]), -B[2])
        fill!((@view ∇H[:, :, :, 3]), -B[3])

        # Spin coupling term
        @inline add_spin_coupling_term_simd!(∇H, S, J, nx, ny, nz)

        # Anisotropic term
        for z in 1:nz, y in 1:ny, x in 1:nx
            ∇H[x, y, z, 3] -= 2dz * S[x, y, z, 3]
        end
    end
end

function range_chunks(r, n)
    rv = @inbounds r[1:n:end-(n-1)]
    rr = @inbounds r[length(rv)*n+1:end]
    rv, rr
end

function add_spin_coupling_term_simd!(∇H, S, J, nx, ny, nz)
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

    N = 8
    lane = VecRange{N}(0)
    xv, xr = range_chunks(1:nx-1, N)

    for z in 1:nz-1
        for y in 1:ny-1
            for x in xv
                xl = lane + x
                @inline add_grad_contrib!(∇H, S, J, xl, y, z, xl + 1, y, z)
                @inline add_grad_contrib!(∇H, S, J, xl, y, z, xl, y + 1, z)
                @inline add_grad_contrib!(∇H, S, J, xl, y, z, xl, y, z + 1)
            end

            for x in xr
                @inline add_grad_contrib!(∇H, S, J, x, y, z, x + 1, y, z)
                @inline add_grad_contrib!(∇H, S, J, x, y, z, x, y + 1, z)
                @inline add_grad_contrib!(∇H, S, J, x, y, z, x, y, z + 1)
            end

            @inline add_grad_contrib!(∇H, S, J, nx, y, z, 1, y, z)
            @inline add_grad_contrib!(∇H, S, J, nx, y, z, nx, y + 1, z)
            @inline add_grad_contrib!(∇H, S, J, nx, y, z, nx, y, z + 1)
        end

        for x in xv
            xl = lane + x
            @inline add_grad_contrib!(∇H, S, J, xl, ny, z, xl + 1, ny, z)
            @inline add_grad_contrib!(∇H, S, J, xl, ny, z, xl, 1, z)
            @inline add_grad_contrib!(∇H, S, J, xl, ny, z, xl, ny, z + 1)
        end

        for x in xr
            @inline add_grad_contrib!(∇H, S, J, x, ny, z, x + 1, ny, z)
            @inline add_grad_contrib!(∇H, S, J, x, ny, z, x, 1, z)
            @inline add_grad_contrib!(∇H, S, J, x, ny, z, x, ny, z + 1)
        end

        @inline add_grad_contrib!(∇H, S, J, nx, ny, z, 1, ny, z)
        @inline add_grad_contrib!(∇H, S, J, nx, ny, z, nx, 1, z)
        @inline add_grad_contrib!(∇H, S, J, nx, ny, z, nx, ny, z + 1)
    end

    for y in 1:ny-1
        for x in xv
            xl = lane + x
            @inline add_grad_contrib!(∇H, S, J, xl, y, nz, xl + 1, y, nz)
            @inline add_grad_contrib!(∇H, S, J, xl, y, nz, xl, y + 1, nz)
            @inline add_grad_contrib!(∇H, S, J, xl, y, nz, xl, y, 1)
        end

        for x in xr
            @inline add_grad_contrib!(∇H, S, J, x, y, nz, x + 1, y, nz)
            @inline add_grad_contrib!(∇H, S, J, x, y, nz, x, y + 1, nz)
            @inline add_grad_contrib!(∇H, S, J, x, y, nz, x, y, 1)
        end

        @inline add_grad_contrib!(∇H, S, J, nx, y, nz, 1, y, nz)
        @inline add_grad_contrib!(∇H, S, J, nx, y, nz, nx, y + 1, nz)
        @inline add_grad_contrib!(∇H, S, J, nx, y, nz, nx, y, 1)
    end

    for x in xv
        xl = lane + x
        @inline add_grad_contrib!(∇H, S, J, xl, ny, nz, xl + 1, ny, nz)
        @inline add_grad_contrib!(∇H, S, J, xl, ny, nz, xl, 1, nz)
        @inline add_grad_contrib!(∇H, S, J, xl, ny, nz, xl, ny, 1)
    end

    for x in xr
        @inline add_grad_contrib!(∇H, S, J, x, ny, nz, x + 1, ny, nz)
        @inline add_grad_contrib!(∇H, S, J, x, ny, nz, x, 1, nz)
        @inline add_grad_contrib!(∇H, S, J, x, ny, nz, x, ny, 1)
    end

    @inline add_grad_contrib!(∇H, S, J, nx, ny, nz, 1, ny, nz)
    @inline add_grad_contrib!(∇H, S, J, nx, ny, nz, nx, 1, nz)
    @inline add_grad_contrib!(∇H, S, J, nx, ny, nz, nx, ny, 1)
end

function add_thermal_noise!(∇H, c)
    for i in eachindex(∇H)
        ∇H[i] += c * randn()
    end
end

function cross_spin!(∇H, S)
    nx, ny, nz, _ = size(S)

    @inbounds @fastmath begin
        @simd for z in 1:nz
            @simd for y in 1:ny
                @simd for x in 1:nx
                    Sx = S[x, y, z, 1]
                    Sy = S[x, y, z, 2]
                    Sz = S[x, y, z, 3]

                    Hx = ∇H[x, y, z, 1]
                    Hy = ∇H[x, y, z, 2]
                    Hz = ∇H[x, y, z, 3]

                    ∇H[x, y, z, 1] = Sy * Hz - Sz * Hy
                    ∇H[x, y, z, 2] = Sz * Hx - Sx * Hz
                    ∇H[x, y, z, 3] = Sx * Hy - Sy * Hx
                end
            end
        end
    end
end

function cross_spin_simd!(∇H, S)
    nx, ny, nz, _ = size(S)
    n = nx * ny * nz

    Hxs = @view ∇H[:, :, :, 1]
    Hys = @view ∇H[:, :, :, 2]
    Hzs = @view ∇H[:, :, :, 3]

    Sxs = @view S[:, :, :, 1]
    Sys = @view S[:, :, :, 2]
    Szs = @view S[:, :, :, 3]

    N = 8
    lane = VecRange{N}(0)
    iv, ir = range_chunks(1:n, N)

    @inbounds @fastmath begin
        for i in iv
            il = lane + i
            Sx = Sxs[il]
            Sy = Sys[il]
            Sz = Szs[il]

            Hx = Hxs[il]
            Hy = Hys[il]
            Hz = Hzs[il]

            Hxs[il] = Sy * Hz - Sz * Hy
            Hys[il] = Sz * Hx - Sx * Hz
            Hzs[il] = Sx * Hy - Sy * Hx
        end

        for i in ir
            Sx = Sxs[i]
            Sy = Sys[i]
            Sz = Szs[i]

            Hx = Hxs[i]
            Hy = Hys[i]
            Hz = Hzs[i]

            Hxs[i] = Sy * Hz - Sz * Hy
            Hys[i] = Sz * Hx - Sx * Hz
            Hzs[i] = Sx * Hy - Sy * Hx
        end
    end
end
