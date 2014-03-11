#!/usr/bin/env julia

using Compose, Color, DataFrames, Base.Collections

import Base: isless

function weighted_color_mean{S <: Number}(
        cs::AbstractArray{LAB,1}, ws::AbstractArray{S,1})
    l = 0.0
    a = 0.0
    b = 0.0
    sumws = sum(ws)
    for (c, w) in zip(cs, ws)
        w /= sumws
        l += w * c.l
        a += w * c.a
        b += w * c.b
    end
    LAB(l, a, b)
end


# Then functions return functions suitable for ContinuousColorScales.
function lab_gradient(cs::ColorValue...)
    if length(cs) < 2
        error("Two or more colors are needed for gradients")
    end

    cs_lab = [convert(LAB, c) for c in cs]
    n = length(cs_lab)
    function f(p::Float64)
        @assert 0.0 <= p <= 1.0
        i = 1 + min(n - 2, max(0, int(floor(p*(n-1)))))
        w = p*(n-1) + 1 - i
        weighted_color_mean([cs_lab[i], cs_lab[i+1]], [1.0 - w, w])
    end
    f
end


#colorer = lab_gradient(color("#001F3F"), color("#FFDC01"))
colorer = lab_gradient(
    LCHab(84,85,278),
    LCHab(84,85,87))
#color("#001F3F"), color("#FFDC01"))
#colorer = lab_gradient(color("grey20"), color("white"))


immutable Exon
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


# Parse the exons from a GTF corresponding to the specified gene id
function parse_gtf_gene_exons(filename, transcript_ids)
    transcript_id_pat = r"transcript_id\s+\"?([^\s]+)\"?;"

    exons = Exon[]

    for line in eachline(open(filename))
        row = split(line, '\t')

        if row[3] != "exon"
            continue
        end

        m = match(transcript_id_pat, row[9])
        @assert m != nothing
        row_transcript_id = m.captures[1]

        if !(row_transcript_id in transcript_ids)
            continue
        end

        push!(exons, Exon(parseint(row[4]), parseint(row[5]), row_transcript_id))
    end

    return exons
end


immutable Feature
    ids::Set{String}
    startpos::Int
    endpos::Int
end


# Figure out nodes and edges of the graph from a sorted array of exons.
function exon_graph_features(exon_array)
    i = 1
    exons = copy(exon_array)
    heapify!(exons)

    features = Feature[]
    exon_group = Exon[]
    while !isempty(exons)
        exon = heappop!(exons)
        startpos = exon.startpos
        push!(exon_group, exon)

        while !isempty(exons) && exons[1].startpos == startpos
            push!(exon_group, heappop!(exons))
        end

        endpos = minimum([exon.endpos for exon in exon_group])
        if !isempty(exons) && exons[1].startpos - 1 < endpos
            endpos = exons[1].startpos - 1
        end

        group_ids = Set{String}([exon.transcript_id for exon in exon_group]...)
        while !isempty(exon_group)
            exon = pop!(exon_group)
            if exon.endpos > endpos
                heappush!(exons, Exon(endpos + 1, exon.endpos, exon.transcript_id))
            end
        end

        push!(features, Feature(group_ids, startpos, endpos))
    end

    return features
end


function draw_graph(features, transcript_abundance)
    n = length(features)

    function nodeangle(i)
        p = 0.9
        p * 2pi * i/n + pi/2 + (1.0 - p)/2 * 2pi
    end

    function nodepos(i)
        theta = nodeangle(i)
        cos(theta), sin(theta)
    end

    # draw nodes
    node_canvases = Array(Compose.Canvas, n)
    nodesize = 4mm
    for (i, feature) in enumerate(features)
        theta = nodeangle(i)
        feature_abundance = 0.0
        for tid in feature.ids
            feature_abundance += transcript_abundance[tid]
        end
        x, y = nodepos(i)
        #node_forms[i] = compose(circle(x, y, 2mm), fill(colorer(feature_abundance)))
        node_canvases[i] = compose(
            canvas(0w, 0h, 1w, 1h, units_inherited=true,
                   rotation=Rotation(theta, x, y)),
            rectangle(x*cx - nodesize/2, y*cy - nodesize/2,
                      nodesize, nodesize),
            stroke("grey30"),
            fill(colorer(feature_abundance)))
    end

    # draw edges
    edge_forms = Compose.Form[]
    for i in 1:n
        ids_i = copy(features[i].ids)
        for j in (i+1):n
            if isempty(ids_i)
                break
            end
            ids_j = features[j].ids

            ids = intersect(ids_i, ids_j)
            setdiff!(ids_i, ids_j)
            if isempty(ids)
                continue
            end

            x0, y0 = nodepos(i)
            x1, y1 = nodepos(j)

            if i + 1 == j
                r = 1.0
            else
                r = 1.0 - 2*(j - i)/n
            end

            theta = (nodeangle(i) + nodeangle(j)) / 2
            ctrl = (r * cos(theta), r * sin(theta))

            line_abundance = sum([transcript_abundance[id] for id in ids])
            isintron = features[i].endpos + 1 < features[j].startpos

            push!(edge_forms,
                  compose(curve((x0, y0), ctrl, ctrl, (x1, y1)),
                          isintron ? linewidth(0.5mm) : linewidth(1.5mm),
                          stroke(colorer(line_abundance))))
        end
    end
    #form = Compose.combine(
        #Compose.combine(node_forms...),
        #Compose.combine(edge_forms...))
    return set_box(
        pad_inner(compose(
            canvas(unit_box=UnitBox(-1, -1, 2, 2)),
            Compose.combine(edge_forms...),
            node_canvases...),

        5mm),
        BoundingBox(0, 0, 1w, 1h))
end


function draw_scale(width)
    step = 0.01
    return compose(
        canvas(0, 0, width, 1h),
        compose(canvas(1cm, 2cm, 1cm, 10cm),
                Compose.combine([compose(rectangle(0, k, 1, step), fill(colorer(k)))
                                 for k in 0:step:1]...),
                svgattribute("shape-rendering", "crispEdges")))
end


function main()
    if length(ARGS) < 3
        println(STDERR, "Usage: isolator-subway-plot.jl genes.gtf transcript_id condition-splicing.tsv")
        return
    end

    gtf_filename = ARGS[1]
    transcript_id = ARGS[2]
    condition_splicing_filename = ARGS[3]

    splicing = readtable(condition_splicing_filename, separator='\t', header=true)
    splicing_tid_row = splicing[splicing[:transcript_id] .== transcript_id, :]
    if size(splicing_tid_row, 1) == 0
        error("Transcript $(transcript_id) not found in the isolator output.")
    end

    transcript_ids = splicing_tid_row[:transcription_group][1]
    println(STDERR, "transcript_ids: ", transcript_ids)

    transcript_id_set = Set(split(transcript_ids, ',')...)

    print(STDERR, "Parsing GTF file...")
    exons = parse_gtf_gene_exons(gtf_filename, transcript_id_set)
    println(STDERR, "done. (", length(exons), " exons)")

    features = exon_graph_features(exons)

    # TODO: we should plot multiple conditions side-by-side
    # This is awkward: one tgroup doesn't correspond to one gene_id.

    splicing = splicing[splicing[:transcription_group] .== transcript_ids, :]
    graphs = {}
    for i in 5:size(splicing, 2)
        splicing[i] ./= sum(splicing[i])
        transcript_abundance =
            [tid => prop for (tid, prop) in zip(splicing[:transcript_id], splicing[i])]
        push!(graphs, draw_graph(features, transcript_abundance))
    end

    scale = draw_scale((1/6)w)

    draw(SVG("subway-plot.svg", 6inch * length(graphs) + (1/6)inch, 6inch),
         compose(canvas(), (rectangle(), fill("grey20")), hstack(graphs..., scale)))
end


main()



