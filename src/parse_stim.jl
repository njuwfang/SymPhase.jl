_takewhile_(f) = (s, i=0) -> begin
    isempty(s) && return i
    i += 1
    while i <= lastindex(s)
        c = s[i]
        if !f(c)
            break
        end
        i += 1
    end
    i-1
end

take_ws = _takewhile_(c -> c == ' ' || c == '\t')
@eval @inline _isdigit(c) = $(Meta.parse(join(["c == '$(j)'" for j in 0:9], " || ")))
@eval @inline _isletter(c) = $(Meta.parse(join([["c == '$('A'+j)'" for j in 0:25];["c == '$('a'+j)'" for j in 0:25]], " || ")))

take_digits = _takewhile_(isdigit)
take_aa = _takewhile_(c -> c != ' ' && c != '\t')

#=function take_ws(s::String, i=0)
    isempty(s) && return i
    i += 1
    while i <= lastindex(s)
        c = s[i]
        if !(c == ' ' || c == '\t')
            break
        end
        i += 1
    end
    i - 1
end=#

function take_id(s, i=0)
    isempty(s) && return i
    i += 1
    if i <= lastindex(s)
        c = s[i]
        if !_isletter(c)
            return i-1
        end
    end
    while i <= lastindex(s)
        c = s[i]
        if !(_isletter(c) || isdigit(c) || c == '_')
            break
        end
        i += 1
    end
    i-1
end

@inline function take_arg(s, i=0)
    isempty(s) && return i
    i += 1
    if i <= lastindex(s)
        c = s[i]
        if c != '('
            return i-1
        end
    end
    while i <= lastindex(s)
        c = s[i]
        if c == ')'
            return i
        end
        i += 1
    end
    i-1
end

@eval Base.@propagate_inbounds @inline function SetGateType!(ops, i, k::AbstractString)
    $(
        [quote if k == $(key) return ops[i] = $(value) end end for (key,value) in GateDict]...
    )
end

function parse_targets!(line::String, i_start, circuit)
    n = length(line)
    flag = true
    i = i_start
    jj = 0
    target_ind = circuit.nops == 1 ? 0 : circuit.target_inds[circuit.nops-1]
    new_target = 0
    for j in i_start:n
        c = line[j]
        if c == '#'
            jj = j
            break
        end
        if (c == ' ' || c == '\t')
            if !flag
                target_ind += 1
                if target_ind > size(circuit.targets, 1)
                    circuit.targets = _double_(circuit.targets)
                end
                @inbounds circuit.targets[target_ind] = new_target = parse(UInt16, SubString(line, i, j-1)) + 1
                circuit.nqubits = max(circuit.nqubits, new_target)
                add_symbol!(circuit)
                flag = true
            end
        else
            if flag
                i = j
                flag = false
            end
        end
    end
    if !flag
        target_ind += 1
        if target_ind > size(circuit.targets, 1)
            circuit.targets = _double_(circuit.targets)
        end
        if jj != 0
            @inbounds circuit.targets[target_ind] = new_target = parse(UInt16, SubString(line, i, jj-1))+1
        else
            @inbounds circuit.targets[target_ind] = new_target = parse(UInt16, SubString(line, i, n))+1
        end
        circuit.nqubits = max(circuit.nqubits, new_target)
        add_symbol!(circuit)
    end
    circuit.target_inds[circuit.nops] = target_ind
    nothing
end

function parse_stim(stim_file_name::String)
    circuit = Circuit()
    i_start = 0
    i_end = 0
    for line in eachline(stim_file_name)
        i_start = take_ws(line)
        i_end = take_id(line, i_start)
        if i_start != i_end
            if circuit.nops == size(circuit.ops, 1)
                circuit.ops = _double_(circuit.ops)
                circuit.args = _double_(circuit.args)
                circuit.target_inds = _double_(circuit.target_inds)
            end

            circuit.nops += 1
            @inbounds SetGateType!(circuit.ops, circuit.nops, SubString(line, i_start+1, i_end))
                
            i_start = take_ws(line, i_end)
            i_end = take_arg(line, i_start)
            if i_end-i_start >= 3
                circuit.args[circuit.nops] = parse(Float32, SubString(line, i_start+2, i_end-1))
            end

            parse_targets!(line, i_end+1, circuit)
        end
    end
    #circuit.nqubits += 1
    circuit
end