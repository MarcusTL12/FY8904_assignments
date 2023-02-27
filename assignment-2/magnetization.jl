using Statistics

function space_averaged_magnetization(S)
    mean(@view S[:, :, :, 3])
end
