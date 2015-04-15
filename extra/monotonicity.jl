#!/usr/bin/env julia

# Compute probabilities that feature-level splicing changes are occuring
# monotonically across time points.
#
# As input it it takes a YAML formatted file listing timepoints. Each timepoint
# is a just a list of condition names, which should match the condition names in
# the experiment file fed to isolator. E.g.
#
#     - [h7-day20, empty-vector]
#     - [adult, h7-1yr, let-7g-oe, fetal-ventricle, fetal-atrium, IMR90]
#
# Any condition that was present in the experiment but not present here will be
# ignored when considering monotonicity.
#
# This will print a table containing any feature with > 0.5 probability of
# changing across timepoints in a specific pattern.


using DataStructures, DataFrames, YAML

if length(ARGS) != 1
    println(STDERR,
    """
        Usage: isolator summarize condition-feature-splicing-samples isolator-output.h5 | ./monotonicity.jl timepoints.yml
    """)
    exit(1)
end


function read_permutations(timepoints_filename)
    equiv_classes = YAML.load_file(timepoints_filename)

    header = split(strip(readline(STDIN)), '\t')
    condition_names = [replace(name, r"_splice_rate_samples$", "") for name in header[7:end]]

    equiv = Dict{Int, Int}()
    for (i, equiv_class) in enumerate(equiv_classes)
        for cond in equiv_class
            equiv[findfirst(condition_names, cond)] = i
        end
    end

    k = maximum(values(equiv))
    permutations = Vector{Uint64}[]

    for i in 1:length(condition_names)
        if !haskey(equiv, i)
            equiv[i] = 0
        end
    end

    features = Any[]
    for line in eachline(STDIN)
        entry = split(strip(line), '\t')
        push!(features, entry[1:6])
        rows = [split(row, ',') for row in entry[7:end]]

        M = [parse(Float32, rows[i][j])
             for i in 1:length(rows), j in 1:length(rows[1])]

        permutations_row = Array(Uint64, length(rows[1]))

        for j in 1:size(M, 2)
            ps = sortperm(M[:,j])
            for i in 1:length(ps)
                ps[i] = equiv[ps[i]]
            end

            u = Uint64(0)
            for p in ps
                if p == 0
                    continue
                end
                u *= k
                u += p
            end

            permutations_row[j] = u
        end
        push!(permutations, permutations_row)
    end

    return permutations, features, k
end

permutations, features, k = read_permutations(ARGS[1])


function decode_permutation(p, k, n)
    decoded = Array(Uint64, n)
    for i in 1:n
        decoded[n-i+1] = p % k
        p = div(p, k)
    end
    return decoded
end


function interesting_transcripts(permutations, features)
    num_samples = length(permutations[1])
    interesting = Any[]
    cutoff = 0.5
    tally = DefaultDict(Uint64, Uint64, () -> 0)
    for i in 1:size(permutations, 1)
        for key in keys(tally)
            tally[key] = 0
        end

        for j in 1:num_samples
            tally[permutations[i][j]] += 1
        end

        for (k, v) in tally
            freq = v / num_samples
            if freq > cutoff
                push!(interesting, (features[i], k, freq))
            end
        end
    end

    return interesting
end



interesting = interesting_transcripts(permutations, features)

println("gene_names\tgene_ids\tincluding_transcript_ids\texcluding_transcript_ids\tlocus\tfeature_type\tprobability")
for feature in interesting
    println(join(feature[1], '\t'), '\t', feature[3])
end


