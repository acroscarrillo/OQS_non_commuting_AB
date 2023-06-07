using SparseArrays

a = 2

const sigma_x = [[0.0 1.0]; [1.0 0.0]]
const sigma_y = [[0.0 -1.0 * im]; [1.0 * im 0.0]]
const sigma_z = [[1.0 0.0]; [0.0 -1.0]]
const sigma_p = 0.5 * (sigma_x + 1 * im * sigma_y)
const sigma_m = sigma_p'
const sigma_0 = [[1.0 0.0]; [0.0 1.0]]

⨷(a, b) = kron(a, b)

a=2

function pauli_string(L, j)
    if j > L
        throw(ArgumentError("Requested pauli string at site j with j>L."))
    end

    if j > 1
        string = sparse(-sigma_z)
        for _ in 2:(j - 1)
            string = string ⨷ (-sigma_z)
        end
        for _ in j:L
            string = string ⨷ sigma_0
        end
        return string
    else
        return sparse(I(2^L))
    end
end

function create(L, j)
    I_left = sparse(I(2^(j - 1)))
    I_right = sparse(I(2^(L - j)))
    return pauli_string(L, j) * (I_left ⨷ sparse(sigma_p) ⨷ I_right)
end

function destroy(L, j)
    return create(L, j)'
end

function pwr_law_mat(L, p, real=true)
    if real
        temp = zeros(Float64, (L, L))
        for col in 1:L
            for row in 1:L
                temp[row, col] = randn() / (abs(row - col) + 1)^p
            end
        end
        temp_sym = temp * (temp)'
        max_lamb = eigs(temp_sym; nev=1, ritzvec=false, which=:LR)[1][1]
        return temp_sym / max_lamb
    else
        temp = zeros(ComplexF64, (L, L))
        for col in 1:L
            for row in 1:L
                temp[row, col] = (randn() + im * randn()) / (abs(row - col) + 1)^p
            end
        end
        temp_sym = temp * (temp)'
        max_lamb = eigs(temp_sym; nev=1, ritzvec=false, which=:LR)[1][1]
        return temp_sym / max_lamb
    end
end