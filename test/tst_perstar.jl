module tst_cov
    using Test
    using cloudCovErr

    @testset "PerStar" begin
        ttt_residimIn = ones(51,51)
        ttt_maskim = zeros(Bool,51,51)
        ttt_maskim[26,26] = true
        ttt_starim = ones(51,51)
        data_in, stars_in, kmasked2d = stamp_cutter(26,26,ttt_residimIn,ttt_starim,ttt_maskim)
        @test kmasked2d[17,17]

        function psfmodel_test(x,y,Np)
            psftest = [
             0.02  0.02  0.03  0.03  0.04  0.04  0.04  0.03  0.03  0.02  0.02;
             0.02  0.03  0.04  0.05  0.06  0.06  0.06  0.05  0.04  0.03  0.02;
             0.03  0.04  0.05  0.07  0.09  0.1   0.09  0.07  0.05  0.04  0.03;
             0.03  0.05  0.07  0.11  0.17  0.2   0.17  0.11  0.07  0.05  0.03;
             0.04  0.06  0.09  0.17  0.33  0.5   0.33  0.17  0.09  0.06  0.04;
             0.04  0.06  0.1   0.2   0.5   1.0   0.5   0.2   0.1   0.06  0.04;
             0.04  0.06  0.09  0.17  0.33  0.5   0.33  0.17  0.09  0.06  0.04;
             0.03  0.05  0.07  0.11  0.17  0.2   0.17  0.11  0.07  0.05  0.03;
             0.03  0.04  0.05  0.07  0.09  0.1   0.09  0.07  0.05  0.04  0.03;
             0.02  0.03  0.04  0.05  0.06  0.06  0.06  0.05  0.04  0.03  0.02;
             0.02  0.02  0.03  0.03  0.04  0.04  0.04  0.03  0.03  0.02  0.02;
            ]
            return psftest
        end
        kmasked2d = zeros(Bool,11,11)
        kmasked2d[6,end] = true
        kmasked2d[6,6] = true
        circmask = kstar_circle_mask(11;rlim=25)
        x_star = 10
        y_star = 10
        flux_star = 100
        psft, kstar, kpsf2d, kcond0, kcond, kpred, dnt = gen_pix_mask(kmasked2d,psfmodel_test,circmask,x_star,y_star,flux_star,Np=11,thr=2000)
        @test reshape(kstar,11,11)[6,end]
        @test !reshape(kstar,11,11)[6,1]
        @test (kcond0, kcond, kpred, dnt) == (71, 71, 9, 0)
        psft, kstar, kpsf2d, kcond0, kcond, kpred, dnt = gen_pix_mask(kmasked2d,psfmodel_test,circmask,x_star,y_star,flux_star,Np=11,thr=200)
        @test !reshape(kstar,11,11)[1,1]
        @test !reshape(kpsf2d,11,11)[1,1]
        function psfmodel_test(x,y,Np)
            psftest = [
             0.02  0.02  0.03  0.03  0.04  0.04  0.04  0.03  0.03  0.02  0.02;
             0.02  0.03  0.04  0.05  0.06  0.06  0.06  0.05  0.04  0.03  0.02;
             0.03  0.04  0.05  0.07  0.09  0.1   0.09  0.07  0.05  0.04  0.03;
             0.03  0.05  0.07  0.11  0.17  0.2   0.17  0.11  0.07  0.05  0.03;
             0.04  0.06  0.09  0.17  0.33  0.5   0.33  0.17  0.09  0.06  0.04;
             0.04  0.06  0.1   0.2   0.5   1.0   0.5   0.2   0.1   0.06  0.04;
             0.04  0.06  0.09  0.17  0.33  0.5   0.33  0.17  0.09  0.06  0.04;
             0.03  0.05  0.07  0.11  0.17  0.2   0.17  0.11  0.07  0.05  0.03;
             0.03  0.04  0.05  0.07  0.09  0.1   0.09  0.07  0.05  0.04  0.03;
             0.02  0.03  0.04  0.05  0.06  0.06  0.06  0.05  0.04  0.03  0.02;
             -0.2  0.02  0.03  0.03  0.04  0.04  0.04  0.03  0.03  0.02  0.02;
            ]
            return psftest
        end
        flux_star = 1e5
        psft, kstar, kpsf2d, kcond0, kcond, kpred, dnt = gen_pix_mask(kmasked2d,psfmodel_test,circmask,x_star,y_star,flux_star,Np=11,thr=200)
        @test (dnt .& 2^5) != 0

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
        cov1 = zeros(T,Np*Np,Np*Np)
        μ1 = zeros(T,Np*Np)
        cov_avg!(bimage, ism, bism, in_subimage, Np=Np, widx=widx, widy=widy)
        build_cov!(cov1,μ1,6,6,bimage,bism,Np,widx,widy)
        @test abs(cov1[1,1]-271.500183105)<1e-7
        @test abs(cov1[1,end]-271.499633789)<1e-7
        @test abs(cov1[end,end]-271.500732421875)<1e-7

        kpsf2d=zeros(Bool,3,3)
        kpsf2d[2,2] = true
        kpsf2d[2,3] = true
        km = kpsf2d[:]
        km[1] = true
        km[end] = true
        data_in = [
          0.516654    0.128205   1.94308;
         -0.0805772  -0.920908   0.212453;
         -0.774471    0.165229  -0.363535;
        ]
        stars_in = 2*ones(3,3)
        cov_loc = [
        47.2245   1.33946   0.628186   0.841306   0.288469  -0.437706   0.754434  -0.245601  -0.150857;
         0.0     47.2602    1.39354    0.570245   0.871201   0.243023  -0.457453   0.806491  -0.209147;
         0.0      0.0      47.166      1.40991    0.471935   0.767516   0.251835  -0.398556   0.772153;
         0.0      0.0       0.0       47.1832     1.44052    0.461676   0.757649   0.282371  -0.314066;
         0.0      0.0       0.0        0.0       47.2298     1.34296    0.471258   0.766785   0.218381;
         0.0      0.0       0.0        0.0        0.0       47.2845     1.36993    0.492573   0.715498;
         0.0      0.0       0.0        0.0        0.0        0.0       47.1459     1.3015     0.425096;
         0.0      0.0       0.0        0.0        0.0        0.0        0.0       47.1451     1.30432;
         0.0      0.0       0.0        0.0        0.0        0.0        0.0        0.0       47.2068;
        ]
        μ = [0.005037597380578518
             0.003629490500316024
             0.0026571564376354218
             0.0032599247060716152
             0.002772299339994788
             0.0034678715746849775
             0.0051109688356518745
             0.0037797277327626944
             0.003628455102443695]
        psft = [
        0.0274769  0.0360258  0.032605;
        0.0330318  0.0403823  0.0339776;
        0.0301219  0.033388   0.0255433;
        ]
        out = condCovEst_wdiag(cov_loc,μ,km,kpsf2d,data_in,stars_in,psft;Np=3,export_mean=false,n_draw=0,diag_on=true)
        @test length(out) == 1
        @test length(out[1]) == 6
        @test abs.(out[1][1] .- 133.81487945671913) .< 1e-7
        out = condCovEst_wdiag(cov_loc,μ,km,kpsf2d,data_in,stars_in,psft;Np=3,export_mean=true,n_draw=0,diag_on=true)
        @test length(out) == 2
        @test size(out[2]) == (3,3)
        out = condCovEst_wdiag(cov_loc,μ,km,kpsf2d,data_in,stars_in,psft;Np=3,export_mean=true,n_draw=2,diag_on=true)
        @test length(out) == 3
        @test size(out[3]) == (9,2)
    end
end
