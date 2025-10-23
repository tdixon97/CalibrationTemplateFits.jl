using Pkg
using CalibrationTemplateFits
using BAT
using ArgParse
using PropDicts
using YAML
using Glob
using IntervalSets

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--config", "-c"
        help = "Path to YAML configuration file"
        arg_type = String
        required = true

        "--output", "-o"
        help = "Path to output folder"
        arg_type = String
        required = true

        "--pos", "-p"
        help = "Position to use"
        arg_type = Int
        required = true

        "--binning", "-b"
        help = "Binning to use, either range or list"
        arg_type = String
        required = false
        default = "2605:20:2625"

        "--vary-fccd", "-f"
        help = "Vary the FCCD"
        action = :store_true

        "--read-samples", "-s"
        help = "Only read the samples (to make summary plots)"
        action = :store_true

    end
    return parse_args(s)
end

function main()
    @info "running using ", Base.Threads.nthreads(), " threads"

    args = parse_commandline()
    cfg = readprops(args["config"])
    s = YAML.write(Dict(cfg))

    @info "Using config \n$s"

    pos = args["pos"]
    dets = YAML.load_file(cfg.det_list)
    @info "... using detectors $dets"

    spec = cfg.spec
    # read data
    @info "... read the data"
    binning = parse_binning(args["binning"])
    data_hists = read_data_histograms(cfg.data_path, "hist/$spec", dets, binning)


    list =
        haskey(cfg, "file-list") ? cfg["file-list"] :
        basename.(glob("*.lh5", cfg.data_path))

    @info "... read mc"
    models = read_models_hist(
        dets,
        glob(cfg.mc_label*"_"*"$pos*", cfg.mc_path),
        binning,
        r".*z-offset_([-\d.]+)_phi-offset_([-\d.]+)",
        "hist/$spec",
    )


    @info "... make likelihood"
    likelihood =
        build_likelihood(data_hists, models, n_sim = cfg.n_sim, livetime = cfg.livetime)

    # this can be in config but its hard to keep type stability
    prior = build_prior(dets, vary_fccd = args["vary-fccd"])

    posterior = PosteriorMeasure(likelihood, prior)

    # sample
    @info "... start sampling"

    dir = args["output"]
    isdir(dir) || mkdir(dir)

    if (!args["read-samples"])
        samples = bat_sample(
            posterior,
            TransformedMCMC(proposal = RandomWalk(), nsteps = 10^6, nchains = 4),
        ).result
        @info "... now save samples"

        bat_write(dir*"/samples.h5", samples)
    else
        samples = bat_read(dir*"/samples.h5").result
    end

    @info "... make some summary plots"
    plot_posteriors(dir*"/plots.pdf", samples, dets, args["vary-fccd"])
    plot_reconstruction_makie(
        data_hists,
        models,
        dir*"/best_fit.pdf",
        BAT.mode(samples),
        cfg.livetime/cfg.n_sim,
    )


end


# Run the main function if this file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
