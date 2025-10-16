

using Pkg
using CalibrationTemplateFits
using BAT
using ArgParse
using PropDicts
using YAML
using Glob

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--config", "-c"
        help = "Path to YAML configuration file"
        arg_type = String
        required = true

        "--vary-fccd", "-f"
        help = "Vary the FCCD"
        action = :store_true

        "--vary-dlf", "-t"
        help = "Vary the DLF parameter"
        action = :store_true
    end
    return parse_args(s)
end

function main()
    @info "running using ", Base.Threads.nthreads(), " threads"

    args = parse_commandline()
    cfg = readprops(args["config"])
    s = YAML.write(Dict(cfg))

    if (args["vary-fccd"]|args["vary-dlf"])
        throw(NotImplementedError("vary fccd or dlf is not implemented"))
    end

    @info "Using config \n$s"

    # get the channel map
    @info "... reading channel map"
    chmap = readprops(cfg.channel_map)
    pos = cfg.pos
    dets = [ch for ch in keys(chmap) if chmap[ch].system=="geds"]

    # read data
    @info "... read the data"
    binning = parse_binning(cfg.binning)
    rawid_map = Dict(det=>chmap[det].daq.rawid for det in dets)

    # and the histograms
    data_hists = read_data_histograms(
        cfg.data_path,
        readprops(cfg.file_list)[pos],
        rawid_map,
        binning,
    )

    @info "... read mc"
    models = read_models(
        dets,
        glob(cfg.mc_label*"_"*"$pos*", cfg.mc_path),
        binning,
        r".*z-offset_([-\d.]+)_phi-offset_([-\d.]+)",
    )

    @info "... make likelihood"
    likelihood = build_likelihood(data_hists, models)

    # this can be in config but its hard to keep type stability
    prior = distprod(A = 0 .. 3000, z = -40.0 .. -80.0, φ = -6.0 .. 6.0)

    # sample
    @info "... start sampling"
    samples = bat_sample(
        posterior,
        MCMCSampling(mcalg = MetropolisHastings()),
        nsteps = 10^6,
        nchains = 4,
    ).result

    #@info "... make some summary plots"
    #make_summary_plots(samples)

    # save
    @info "... now save samples"
    bat_write(cfg.output, samples)

end


# Run the main function if this file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
