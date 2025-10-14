

using Pkg
Pkg.activate(".") # activate the environment
Pkg.instantiate()
using CalibrationTemplateFits
using BAT
using ArgParse
using PropDicts
using YAML

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--config", "-c"
        help = "Path to YAML configuration file"
        arg_type = String
        required = true
    end
    return parse_args(s)
end

function main()
    @info "running using ", Base.Threads.nthreads(), " threads"

    args = parse_commandline()
    cfg = readprops(args["config"])
    s = YAML.write(cfg)

    @info "Using config ($cfg)\n$s"

    # get the channel map
    @info "... reading channel map"
    chmap = readprobs(cfg.channel_map)

    # read data
    @info "... read the data"
    file_list = readprops(cfg.file_list)
    data_table = io.read_data(cfg.data_path, file_list)

    # cut some events
    @info "... build data histograms"
    data_sel = data[(.!data.coincident.puls) .& (.!data.trigger.is_forced)]
    binning = io.parse_binning(cfg.binning)

    data_hists = Dict(
        det => io.get_data_histogram(det, data_sel, binning) for
        det in keys(chmap) if chmap[det].system=="geds"
    )

    @info "... read mc"
    models = io.read_models(cfg.mc_path, cfg.mc_label, binning)

    @info "... make likelihood"
    likelihood = likelihood.build_likelihood(data, models)

    # this can be in config but its hard to keep type stability
    prior = distprod(A = 0 .. 3000, z = -40 .. -80, phi = -6 .. 6)

    # sample
    @info "... start sampling"
    samples =
        bat_sample(posterior, MCMCSampling(mcalg = MetropolisHastings(); fit_kwargs...)).result

    # save
    @info "... now save samples"
    bat_write(cfg.output, samples)
end


# Run the main function if this file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
