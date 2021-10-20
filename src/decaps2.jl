push!(LOAD_PATH, "/n/home12/saydjari/finksage/ExtSoftware/disCovErr.jl/src/")
using disCovErr
using ArgParse

function parse_commandline()
    s=ArgParseSettings()
    @add_arg_table s begin
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
        "--thr", "-t"
            help = "threshold for psf-based masking of the residuals (int for ease)"
            arg_type = Int
            default = 20
        "--Np", "-p"
            help = "sidelength size of spatial covariance matrix (must be int, odd)"
            arg_type = Int
            default = 33
        "--ccdlist", "-c"
            help = "list of ccd_CAT names to be run"
            nargs = '+'
            default = String[]
    end
    return parse_args(s)
end

function run_wrapper()
    parg = parse_commandline()
    proc_all(parg["base"],parg["date"],parg["filt"],parg["vers"],
        parg["basecat"],ccdlist=parg["ccdlist"],thr=parg["thr"],Np=parg["Np"])
end

run_wrapper()
