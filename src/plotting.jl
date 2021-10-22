# plotting functions to accompany the repo

function cov_as_img()
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
        vmin=-0.09,
        vmax=0.09
    )
    ax.set_title("Subimage of Cov")
    ax.axes.xaxis.set_ticks([])
    ax.axes.yaxis.set_ticks([])

    cax = fig.add_axes([ax.get_position().x1,ax.get_position().y0,0.02,ax.get_position().height])
    plt.colorbar(sc, cax=cax)
end
