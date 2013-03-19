#!/usr/bin/env julia

# Subway Plots: a visualization for isoform abundance.
#
#

## along an axis
function amap(f::Function, A::StridedArray, axis::Integer)
    dimsA = size(A)
    ndimsA = ndims(A)
    axis_size = dimsA[axis]

    if axis_size == 0
        return f(A)
    end

    idx = ntuple(ndimsA, j -> j == axis ? 1 : 1:dimsA[j])
    r = f(sub(A, idx))
    R = Array(typeof(r), axis_size)
    R[1] = r

    for i = 2:axis_size
        idx = ntuple(ndimsA, j -> j == axis ? i : 1:dimsA[j])
        R[i] = f(sub(A, idx))
    end

    return R
end

using Compose
using Color
using Heap

import Base.isless

# Very simple parsing of GTF files.

module GTF

type Entry
    seqname::String
    name::String
    feature::String
    startpos::Int
    endpos::Int
    score::String
    strand::String
    frame::Union(Int, Nothing)
    attributes::Dict{String,String}
end

const attrib_pat = r"\s*(\S+)\s+(\"([^\"]*)\"|([^\"\s]*))[\s*;]?"

function parse_entry(line)
    line = strip(line)
    if isempty(line)
        return nothing
    end

    fields = split(line, "\t")
    if length(fields) != 9
        error("Malformed GTF")
    end

    attributes = Dict{String,String}()
    for mat in each_match(attrib_pat, fields[9])
        attributes[mat.captures[1]] = mat.captures[3]
    end

    Entry(fields[1], fields[2], fields[3], int(fields[4]), int(fields[5]),
          fields[6], fields[7], fields[8] == "." ? nothing : int(fields[8]),
          attributes)
end

end # module GTF

import GTF

if length(ARGS) < 1
    println("Usage: isolator-subway-plot genes.gtf samples.db gene_id")
    exit(1)
end

genes_fn, samples_fn, gene_id = ARGS[1:3]

print("Reading samples ... ");
cmd = `./isolator-dump-samples.py --gene_id=$gene_id $samples_fn`
transcript_ids = String[]
S = {}
for line in each_line(cmd)
    _, transcript_id, samples = split(strip(line), "\t")
    push(transcript_ids, transcript_id)
    push(S, map(float, split(samples, ",")))
end
S = hcat(S...)
println("done.")

# compute posterior means
posterior_means = amap(sum, S, 2)
posterior_means = posterior_means / sum(posterior_means)
posterior_means = [id => mu for (id, mu) in zip(transcript_ids, posterior_means)]

type Exon
    startpos::Int
    endpos::Int
    transcript_id::String
end


function isless(a::Exon, b::Exon)
    if a.startpos == b.startpos
        a.endpos < b.endpos
    else
        a.startpos < b.startpos
    end
end


exons = Exon[]

print("Reading GTF ... ");
ids = Set{String}()

# This is a hack, since parsing is pretty slow right row. (Grep is much faster.)
pat = @sprintf("gene_id\\s+\"%s\"", gene_id)
genes_f = `grep -P $pat $genes_fn`
for line in each_line(genes_f)
    entry = GTF.parse_entry(line)
    if entry === nothing
        break
    end
    if entry.feature == "exon" &&
       has(entry.attributes, "gene_id") &&
       entry.attributes["gene_id"] == gene_id
        add(ids, entry.attributes["transcript_id"])
        hpush(exons, Exon(entry.startpos - 1, entry.endpos - 1,
                          entry.attributes["transcript_id"]))
    end
end
@printf("done (%d exons)\n", length(exons))

S = Set{String}[]
exon_group = Exon[]
while !isempty(exons)
    exon = hpop(exons)
    startpos = exon.startpos
    push!(exon_group, exon)

    while !isempty(exons) && exons[1].startpos == startpos
        push(exon_group, hpop(exons))
    end

    endpos = min([exon.endpos for exon in exon_group])
    if !isempty(exons) && exons[1].startpos - 1 < endpos
        endpos = exons[1].startpos - 1
    end

    group_ids = Set{String}([exon.transcript_id for exon in exon_group]...)
    while !isempty(exon_group)
        exon = pop!(exon_group)
        if exon.endpos > endpos
            hpush(exons, Exon(endpos + 1, exon.endpos, exon.transcript_id))
        end
    end

    push!(S, group_ids)
end


function lab_rainbow(l, c, h0, n)
    ColorValue[LCHab(l, c, h0 + 360.0 * (i - 1) / n) for i in 1:n]
end


n = length(S)
m = length(ids)

colors = [id => c for (id, c) in zip(ids, lab_rainbow(90, 54, 0, m))]

function avg_colors(cs::LCHab...)
    LCHab(sum([c.l for c in cs]) / length(cs),
          sum([c.c for c in cs]) / length(cs),
          sum([c.h for c in cs]) / length(cs))
end


function weighted_avg_color(ws::Vector{Float64}, cs::LCHab...)
    sumws = sum(ws)
    LCHab(sum([w * c.l for (w, c) in zip(ws, cs)]) / sumws,
          sum([w * c.c for (w, c) in zip(ws, cs)]) / sumws,
          sum([w * c.h for (w, c) in zip(ws, cs)]) / sumws)
end


base_canvas = canvas(Units(-1, -1, 2, 2))

base_canvas <<= text(0, 0, gene_id, hcenter, vcenter) <<
                font("PT Sans") << fontsize(10pt) << fill("grey90") <<
                stroke(nothing)


# Position of the ith node
function nodeangle(i)
    p = 0.9
    p * 2pi * i/n + pi/2 + (1.0 - p)/2 * 2pi
end

function nodepos(i)
    theta = nodeangle(i)
    cos(theta), sin(theta)
end

# draw nodes
for (i, s) in enumerate(S)
    node_abundances = Float64[posterior_means[id] for id in s]
    abundance = sum(node_abundances)
    x, y = nodepos(i)
    c = weighted_avg_color(node_abundances, [colors[id] for id in s]...)
    base_canvas <<= circle(x, y, 1mm * abundance + 0.2mm) << fill(c)
end
base_canvas <<= stroke(nothing)

# draw edges
for i in 1:n
    ids_i = S[i]
    for j in (i+1):n
        if isempty(ids_i)
            break
        end

        ids = intersect(ids_i, S[j])
        ids_i -= S[j]
        if !isempty(ids)
            x0, y0 = nodepos(i)
            x1, y1 = nodepos(j)

            if i + 1 == j
                r = 1.0
            else
                r = 1.0 - 2*(j - i)/n
            end

            theta = (nodeangle(i) + nodeangle(j)) / 2
            ctrl = (r * cos(theta), r * sin(theta))

            line_abundances = Float64[posterior_means[id] for id in ids]
            abundance = sum(line_abundances)
            base_canvas <<=
                curve((x0, y0), ctrl, ctrl, (x1, y1))  <<
                linewidth(1.5mm * abundance) <<
                stroke(weighted_avg_color(line_abundances,
                                          [colors[id] for id in ids]...))
        end
    end
end

img = SVG("try.svg", 4inch, 4inch)
draw(img, canvas() << (rectangle() << fill("grey20")) << pad(base_canvas, 5mm))
finish(img)


