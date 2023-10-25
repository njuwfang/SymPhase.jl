@enum GateType begin
    MX
    MY
    M
    RX
    RY
    R
    CX
    H
    DEPOLARIZE1
    DEPOLARIZE2
    X_ERROR
    Y_ERROR
    Z_ERROR
    I
    X
    Y
    Z
    S
end

GateDict = Dict(
    "MX" => MX,
    "MY" => MY,
    "M" => M,
    "MZ" => M,
    "RX" => RX,
    "RY" => RY,
    "R" => R,
    "RZ" => R,
    "CX" => CX,
    "H" => H,
    "DEPOLARIZE1" => DEPOLARIZE1,
    "DEPOLARIZE2" => DEPOLARIZE2,
    "X_ERROR" => X_ERROR,
    "Y_ERROR" => Y_ERROR,
    "Z_ERROR" => Z_ERROR,
    "I" => I,
    "X" => X,
    "Y" => Y,
    "Z" => Z,
    "S" => S,
)

const _SYMGATE_ = [MX, MY, M, DEPOLARIZE1, DEPOLARIZE2, X_ERROR, Y_ERROR, Z_ERROR]
const _SYMGATE_M_ = [MX, MY, M]
const _SYMGATE_E_ = [X_ERROR, Y_ERROR, Z_ERROR]

mutable struct _Vector{T}
    len::Int
    values::Vector{T}
end

_Vector{T}() where T = _Vector{T}(0, Vector{T}(undef, 128))

function Base.push!(v::_Vector{T}, vv) where T
    vv = convert(T, vv)

    if v.len == length(v.values)
        nv = Vector{T}(undef, v.len<<1)
        @inbounds nv[1:v.len] = v.values
        v.values = nv
    end

    v.len += 1
    @inbounds v.values[v.len] = vv
    v
end

mutable struct Circuit
    nqubits::Int
    nsymbols::Int
    nmeasurements::Int
    nops::Int
    ops::Vector{GateType}
    args::Vector{Float32}
    target_inds::Vector{Int64}
    targets::Vector{Int16}
end

Circuit() = Circuit(0, 0, 0, 0, Vector{GateType}(undef, 128), Vector{Float32}(undef, 128), Vector{Int64}(undef, 128), Vector{Int16}(undef, 128))

@eval function add_symbol!(circuit::Circuit)
    op = circuit.ops[circuit.nops]
    if $(Meta.parse(join(["op == $(op)" for op in _SYMGATE_M_], " || ")))
        circuit.nsymbols += 1
        circuit.nmeasurements += 1
    elseif $(Meta.parse(join(["op == $(op)" for op in _SYMGATE_E_], " || ")))
        circuit.nsymbols += 1
    elseif op == DEPOLARIZE1
        circuit.nsymbols += 2
    elseif op == DEPOLARIZE2
        circuit.nsymbols += 4
    end
end

function _double_(v)
    len = size(v, 1)
    nv = Vector{eltype(v)}(undef, len<<1)
    nv[1:len] = v
    nv
end