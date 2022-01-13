# This is the run function used for commandline access
# for the DECaPS2 survey error bar correction. It can be
# used as a model for other survey specific run structures.

push!(LOAD_PATH, "/n/home12/saydjari/finksage/ExtSoftware/cloudCovErr.jl/src/")
using cloudCovErr
using ArgParse

"""
    parse_commandline()

A simple commandline parser that determines which files are run and makes most
code parameters user-facing.
"""
function parse_commandline()
    s=ArgParseSettings()
    @add_arg_table s begin
        # required file/dir name parameters
        "base"
            required=true
            arg_type=String
            help="base directory containing raw exposure image files"
        "date"
            required=true
            arg_type=String
            help="file name date (exposure datetime)"
        "filt"
            required=true
            arg_type=String
            help="file name filter"
        "vers"
            required=true
            arg_type=String
            help="file name version"
        "basecat"
            required=true
            arg_type=String
            help="base directory containing crowdsource cat files"
        # method parameter values
        "--thr", "-t"
            help = "threshold for psf-based masking of the residuals (int for ease)"
            arg_type = Int
            default = 20
        "--outthr"
            help = "threshold for residual-based masking (int for ease)"
            arg_type = Int
            default = 20000
        "--Np", "-p"
            help = "sidelength size of spatial covariance matrix (must be int, odd)"
            arg_type = Int
            default = 33
        "--wx", "-w"
            help = "sidelength size of sample averaging (must be int, odd)"
            arg_type = Int
            default = 129
        "--wy"
            help = "sidelength size of sample averaging (must be int, odd)"
            arg_type = Int
            default = -1
        "--tilex"
            help = "number of tiles to divide image into along x"
            arg_type = Int
            default = 1
        "--tiley"
            help = "number of tiles to divide image into along y"
            arg_type = Int
            default = -1
        "--cS7"
            help = "when set attempts to correct the amplifier offset in S7 fields"
            action = :store_true
        "--ftype", "-f"
            help = "use float n (32 or 64)"
            arg_type = Int
            default = 32
        # run settings
        "--resume", "-r"
            help = "resume partially completed file where leftoff"
            action = :store_true
        # limit run scope/specific subsets for testing
        "--ccdlist", "-c"
            help = "list of ccd_CAT names to be run"
            nargs = '+'
            default = String[]
    end
    return parse_args(s)
end


"""
    run_wrapper()

Reads commandline arguments and passes them to `proc_all` or start
processing all or a subset of a given exposure.
"""
function run_wrapper()
    parg = parse_commandline()

    # handle default x==y conditions
    if parg["wy"] == -1
        wy = parg["wx"]
    else
        wy = parg["wy"]
    end
    if parg["tiley"] == -1
        tiley = parg["tilex"]
    else
        tiley = parg["tiley"]
    end

    cloudCovErr.proc_all(parg["base"],parg["date"],parg["filt"],parg["vers"],
        parg["basecat"],ccdlist=parg["ccdlist"],resume=parg["resume"],corrects7=parg["cS7"],
        thr=parg["thr"],outthr=parg["outthr"],Np=parg["Np"],widx=parg["wx"],widy=wy,tilex=parg["tilex"],
        tiley=tiley,ftype=parg["ftype"])
end

run_wrapper()
