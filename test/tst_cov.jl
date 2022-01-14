module tst_cov
    using Test
    using cloudCovErr

    @testset "CovEst" begin
        @test outest_bounds([-1,100],100) == 2
        @test outest_bounds([-1,105],100) == 5

        ref = [
        117  216  315  414  513  612  711  810  909;
        126  225  324  423  522  621  720  819  918;
        135  234  333  432  531  630  729  828  927;
        144  243  342  441  540  639  738  837  936;
        153  252  351  450  549  648  747  846  945;
        162  261  360  459  558  657  756  855  954;
        171  270  369  468  567  666  765  864  963;
        180  279  378  477  576  675  774  873  972;
        189  288  387  486  585  684  783  882  981;
        ]

        arr = reshape(1:121,11,11)
        out = zeros(Int64,9,9)
        tot = zeros(size(arr)[1]);
        boxsmooth!(out, arr, tot, 3, 3)
        @test out == ref

        Np = 5
        widx = 3
        widy = 3
        px0 = 1
        py0 = 1
        sx0 = 8
        sy0 = 8
        tilex = 2
        tiley = 2
        T = Float32
        Δx = (widx-1)÷2
        Δy = (widy-1)÷2
        padx = Np+Δx+px0
        pady = Np+Δy+py0
        stepx = (sx0+2) ÷ tilex
        stepy = (sy0+2) ÷ tiley
        in_subimage = zeros(T,stepx+2*padx,stepy+2*pady)
        in_subimage = convert(Array{Float32,2},reshape(1:size(in_subimage)[1]*size(in_subimage)[2],size(in_subimage)))
        ism = zeros(T,stepx+2*padx,stepy+2*pady)
        bimage = zeros(T,stepx+2*padx-2*Δx,stepy+2*pady-2*Δy)
        bism = zeros(T,stepx+2*padx-2*Δx,stepy+2*pady-2*Δy,2*Np-1, Np);
        cov_avg!(bimage, ism, bism, in_subimage, Np=Np, widx=widx, widy=widy)
        @test bimage[end,end] == 3069
        @test ism[1,end] == 21266.0
        @test bism[:,:,1:4,1] == zeros(stepx+2*padx-2*Δx,stepy+2*pady-2*Δy,4)
        @test bism[end,1,1,2] == 19488.0

        T = Float64
        in_subimage = zeros(T,stepx+2*padx,stepy+2*pady)
        in_subimage = convert(Array{Float32,2},reshape(1:size(in_subimage)[1]*size(in_subimage)[2],size(in_subimage)))
        ism = zeros(T,stepx+2*padx,stepy+2*pady)
        bimage = zeros(T,stepx+2*padx-2*Δx,stepy+2*pady-2*Δy)
        bism = zeros(T,stepx+2*padx-2*Δx,stepy+2*pady-2*Δy,2*Np-1, Np);
        cov_avg!(bimage,ism,bism,in_subimage,Np=Np,widx=widx,widy=widy,ftype=64)
        @test bimage[end,end] == 3069
        @test ism[1,end] == 21266.0
        @test bism[:,:,1:4,1] == zeros(stepx+2*padx-2*Δx,stepy+2*pady-2*Δy,4)
        @test bism[end,1,1,2] == 19488.0
    end
end
