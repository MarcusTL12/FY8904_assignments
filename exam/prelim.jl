using Unitful

function prelim()
    n_tert = 4.0^299

    @show n_tert

    total_time = uconvert(u"yr", n_tert * 1u"ps")

    @show total_time

    @show maxintfloat(Float64) maxintfloat(Float32) maxintfloat(Float16)

    @show eps(Float64) eps(Float32) eps(Float16)

    ;
end
