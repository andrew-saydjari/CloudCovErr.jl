# plotting functions to accompany the repo
module plotting

using Conda
using PyCall, LaTeXStrings, Formatting
import PyPlot; const plt = PyPlot

export cov_as_img
export plot_cov_compare

"""
    __int__()

Builds the required python plotting dependency.
"""
function __init__()
    if !haskey(Conda._installed_packages_dict(),"colorcet")
        Conda.add("colorcet",channel="conda-forge")
    end
    plt.matplotlib.style.use("dark_background")
    cc=pyimport("colorcet")
end

function cov_as_img(cov_out)
    (sx, sy) = size(cov_out)
    Np = sqrt(sx)
    Ns = (Np - 1) ÷ 2

    fig = plt.figure(figsize=(16,8), dpi=150)
    plt.subplots_adjust(wspace=0.4,hspace=0.1)
    ax = fig.add_subplot(1,1,1)
    input = reshape(cov_out[Ns+Np*Ns,:],Np,Np)
    sc = ax.imshow(
        input,
        origin="lower",
        interpolation="nearest",
        cmap="cet_fire",
        aspect="equal",
    )
    ax.set_title("Subimage of Cov")
    ax.axes.xaxis.set_ticks([])
    ax.axes.yaxis.set_ticks([])

    cax = fig.add_axes([ax.get_position().x1,ax.get_position().y0,0.02,ax.get_position().height])
    plt.colorbar(sc, cax=cax)
end

function plot_cov_compare(C, cov_out)
    (sx, sy) = size(cov_out)
    Np = sqrt(sx)
    Ns = (Np - 1) ÷ 2

    # Plotting WU
    fig = plt.figure(figsize=(16,8), dpi=150)
    plt.subplots_adjust(wspace=0.4,hspace=0.1)
    ax = fig.add_subplot(2,4,1)
    input = C
    sc = ax.imshow(
        input,
        origin="lower",
        interpolation="nearest",
        cmap="cet_fire",
        aspect="equal",
    )
    ax.set_title("Truth")
    ax.axes.xaxis.set_ticks([])
    ax.axes.yaxis.set_ticks([])

    cax = fig.add_axes([ax.get_position().x1,ax.get_position().y0,0.02,ax.get_position().height])
    plt.colorbar(sc, cax=cax)

    ax = fig.add_subplot(2,4,2)
    input = cov_out
    sc = ax.imshow(
        input,
        origin="lower",
        interpolation="nearest",
        cmap="cet_fire",
        aspect="equal",
    )
    ax.set_title("Estimate")
    ax.axes.xaxis.set_ticks([])
    ax.axes.yaxis.set_ticks([])

    cax = fig.add_axes([ax.get_position().x1,ax.get_position().y0,0.02,ax.get_position().height])
    plt.colorbar(sc, cax=cax)

    ax = fig.add_subplot(2,4,3)
    input = cov_out .- C
    sc = ax.imshow(
        input,
        origin="lower",
        interpolation="nearest",
        cmap="cet_bkr",
        aspect="equal",
        vmin=-0.11,
        vmax=0.11
    )
    ax.set_title("Difference")
    ax.axes.xaxis.set_ticks([])
    ax.axes.yaxis.set_ticks([])

    cax = fig.add_axes([ax.get_position().x1,ax.get_position().y0,0.02,ax.get_position().height])
    plt.colorbar(sc, cax=cax)

    ax = fig.add_subplot(2,4,4)
    input = (cov_out .- C)./C
    sc = ax.imshow(
        input,
        origin="lower",
        interpolation="nearest",
        cmap="cet_bkr",
        aspect="equal",
        vmin=-5,
        vmax=5
    )
    ax.set_title("Percent Difference")
    ax.axes.xaxis.set_ticks([])
    ax.axes.yaxis.set_ticks([])

    cax = fig.add_axes([ax.get_position().x1,ax.get_position().y0,0.02,ax.get_position().height])
    plt.colorbar(sc, cax=cax)

    ax = fig.add_subplot(2,4,5)
    input = reshape(C[Ns+Np*Ns,:],Np,Np)
    sc = ax.imshow(
        input,
        origin="lower",
        interpolation="nearest",
        cmap="cet_fire",
        aspect="equal",
    )
    ax.set_title("Truth")
    ax.axes.xaxis.set_ticks([])
    ax.axes.yaxis.set_ticks([])

    cax = fig.add_axes([ax.get_position().x1,ax.get_position().y0,0.02,ax.get_position().height])
    plt.colorbar(sc, cax=cax)

    ax = fig.add_subplot(2,4,6)
    input = reshape(cov_out[Ns+Np*Ns,:],Np,Np)
    sc = ax.imshow(
        input,
        origin="lower",
        interpolation="nearest",
        cmap="cet_fire",
        aspect="equal",
    )
    ax.set_title("Estimate")
    ax.axes.xaxis.set_ticks([])
    ax.axes.yaxis.set_ticks([])

    cax = fig.add_axes([ax.get_position().x1,ax.get_position().y0,0.02,ax.get_position().height])
    plt.colorbar(sc, cax=cax)

    ax = fig.add_subplot(2,4,7)
    input = reshape(cov_out[Ns+Np*Ns,:],Np,Np) .- reshape(C[Ns+Np*Ns,:],Np,Np)
    sc = ax.imshow(
        input,
        origin="lower",
        interpolation="nearest",
        cmap="cet_bkr",
        aspect="equal",
        vmin=-0.09,
        vmax=0.09
    )
    ax.set_title("Difference")
    ax.axes.xaxis.set_ticks([])
    ax.axes.yaxis.set_ticks([])

    cax = fig.add_axes([ax.get_position().x1,ax.get_position().y0,0.02,ax.get_position().height])
    plt.colorbar(sc, cax=cax)

    ax = fig.add_subplot(2,4,8)
    input = (reshape(cov_out[Ns+Np*Ns,:],Np,Np) .- reshape(C[Ns+Np*Ns,:],Np,Np))./reshape(C[Ns+Np*Ns,:],Np,Np)
    sc = ax.imshow(
        input,
        origin="lower",
        interpolation="nearest",
        cmap="cet_bkr",
        aspect="equal",
        vmin=-0.7,
        vmax=0.7
    )
    ax.set_title("Difference")
    ax.axes.xaxis.set_ticks([])
    ax.axes.yaxis.set_ticks([])

    cax = fig.add_axes([ax.get_position().x1,ax.get_position().y0,0.02,ax.get_position().height])
    plt.colorbar(sc, cax=cax)
end
end
