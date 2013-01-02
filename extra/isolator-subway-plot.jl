#!/usr/bin/env julia

# Subway Plots: a visualization for isoform abundance.
#
#

require("Compose")
using Compose

require("Heap")
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

function readentry(io)
    line = strip(readline(io))
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

print("Reading samples ... "); flush(stdout_stream)
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

print("Reading GTF ... "); flush(stdout_stream)
ids = Set{String}()
genes_f = open(genes_fn, "r")
while !((entry = GTF.readentry(genes_f)) === nothing)
    if entry.feature == "exon" &&
       has(entry.attributes, "gene_id") &&
       entry.attributes["gene_id"] == gene_id
        add(ids, entry.attributes["transcript_id"])
        hpush(exons, Exon(entry.startpos - 1, entry.endpos - 1,
                          entry.attributes["transcript_id"]))
    end
end
println("done.")

S = Set{String}[]
exon_group = Exon[]
while !isempty(exons)
    exon = hpop(exons)
    startpos = exon.startpos
    push(exon_group, exon)

    while !isempty(exons) && exons[1].startpos == startpos
        push(exon_group, hpop(exons))
    end

    endpos = min([exon.endpos for exon in exon_group])
    if !isempty(exons) && exons[1].startpos - 1 < endpos
        endpos = exons[1].startpos - 1
    end

    group_ids = Set{String}([exon.transcript_id for exon in exon_group]...)
    while !isempty(exon_group)
        exon = pop(exon_group)
        if exon.endpos > endpos
            hpush(exons, Exon(endpos + 1, exon.endpos, exon.transcript_id))
        end
    end

    push(S, group_ids)
end


function lab_rainbow(l, c, h0, n)
    Color[LCHab(l, c, h0 + 360.0 * (i - 1) / n) for i in 1:n]
end


n = length(S)
m = length(ids)

colors = [id => c for (id, c) in zip(ids, lab_rainbow(70, 54, 0, m))]

function avg_colors(cs::LCHab...)
    LCHab(sum([c.l for c in cs]) / length(cs),
          sum([c.c for c in cs]) / length(cs),
          sum([c.h for c in cs]) / length(cs))
end

base_canvas = canvas(Units(-1, -1, 2, 2))

for (i, s) in enumerate(S)
    abundance = sum([posterior_means[id] for id in s])
    theta = 0.9 * 2pi * i/n
    x, y = cos(theta), sin(theta)
    c = avg_colors([colors[id] for id in s]...)
    base_canvas <<= circle(x, y, 2mm * abundance + 0.5mm) << fill(c)
end
base_canvas <<= stroke(nothing)

img = SVG("try.svg", 6inch, 6inch)
draw(img, pad(base_canvas, 5mm))
finish(img)



# Is this going to work? Are there nicer alternatives we could do?

# What if we do a matrix of dots with components on one axis and transcripts on
# another.

# This at least shows splicing clearly. Can we show uncertaintly? Plotting
# components just once has the advantage that showing uncertaintly is less
# important.



