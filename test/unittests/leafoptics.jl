using Test
import Emerald.EmeraldLand.LeafOptics as LO
import Emerald.EmeraldLand.Namespace as NS


@testset verbose = true "New Leaf Optics Model" begin
    @testset "Direct radiation ρ and τ of an interface" begin
        n₁ = 1.0
        n₂ = 1.33
        for θ₁ in collect(Float64, 0:1:89)
            ρ,τ = LO.interface_ρ_τ(n₁, n₂, θ₁);
            @test ρ + τ == 1.0;
            @test 0 <= ρ <= 1;
            @test 0 <= τ <= 1;
            ρ,τ = LO.interface_ρ_τ(n₂, n₁, θ₁);
            @test ρ + τ == 1.0;
            @test 0 <= ρ <= 1;
            @test 0 <= τ <= 1;
        end;
    end;

    @testset "Isotropic radiation τ of an interface" begin
        n₁ = 1.0
        n₂ = 1.33
        for θ₁ in collect(Float64, 1:1:90)
            τ = LO.interface_isotropic_τ(n₁, n₂, θ₁);
            @test 0 <= τ <= 1;
            τ = LO.interface_isotropic_τ(n₂, n₁, θ₁);
            @test 0 <= τ <= 1;
        end;
    end;

    @testset "Isotropic radiation τ of leaf sublayer" begin
        config = NS.SPACConfiguration{Float64}(DATASET = NS.LAND_2021_1NM);
        lha = config.LHA;
        bio = NS.HyperspectralLeafBiophysics(config);
        τ₁,f_cab_1,f_car_1 = LO.sublayer_τ(lha, bio, 5.0, 1/bio.MESOPHYLL_N, 10);
        τ₂,f_cab_2,f_car_2 = LO.sublayer_τ(lha, bio, 5.0, 1 - 1/bio.MESOPHYLL_N, 10);

        @test all(0 .< τ₁ .< 1);
        @test all(0 .< τ₂ .< 1);
        @test all(0 .<= f_cab_1 .<= 1);
        @test all(0 .<= f_cab_2 .<= 1);
        @test all(0 .<= f_car_1 .<= 1);
        @test all(0 .<= f_car_2 .<= 1);
    end;

    @testset "Isotropic radiation ρ and τ of leaf layer" begin
        config = NS.SPACConfiguration{Float64}(DATASET = NS.LAND_2021_1NM);
        lha = config.LHA;
        bio = NS.HyperspectralLeafBiophysics(config);
        ρ₁,τ₁ = LO.layer_ρ_τ(lha, bio, 5.0, 1/bio.MESOPHYLL_N, 90.0);
        ρ₂,τ₂ = LO.layer_ρ_τ(lha, bio, 5.0, 1 - 1/bio.MESOPHYLL_N, 90.0);
        @test all(0 .< ρ₁ .< 1);
        @test all(0 .< ρ₂ .< 1);
        @test all(0 .< τ₁ .< 1);
        @test all(0 .< τ₂ .< 1);
        @test all(ρ₁ .+ τ₁ .< 1);
        @test all(ρ₂ .+ τ₂ .< 1);
        ρ₁,τ₁,ρ₂,τ₂ = LO.layer_ρ_τ_diffuse(lha, bio, 5.0);
        @test all(0 .< ρ₁ .< 1);
        @test all(0 .< ρ₂ .< 1);
        @test all(0 .< τ₁ .< 1);
        @test all(0 .< τ₂ .< 1);
        @test all(ρ₁ .+ τ₁ .< 1);
        @test all(ρ₂ .+ τ₂ .< 1);
    end;

    @testset "Isotropic radiation ρ and τ of the leaf" begin
        config = NS.SPACConfiguration{Float64}(DATASET = NS.LAND_2021_1NM);
        lha = config.LHA;
        bio = NS.HyperspectralLeafBiophysics(config);
        ρ₁,τ₁ = LO.leaf_spectra(lha, bio, 5.0, 40.0);
        ρ₂,τ₂ = LO.leaf_spectra(lha, bio, 5.0, 59.0);
        @test all(0 .< ρ₁ .< 1);
        @test all(0 .< ρ₂ .< 1);
        @test all(0 .< τ₁ .< 1);
        @test all(0 .< τ₂ .< 1);
        @test all(ρ₁ .+ τ₁ .< 1);
        @test all(ρ₂ .+ τ₂ .< 1);
    end;

    @testset "Raw SIF excitation coefficient of leaf layer" begin
        config = NS.SPACConfiguration{Float64}(DATASET = NS.LAND_2021_1NM);
        lha = config.LHA;
        wls = config.WLSET;
        bio = NS.HyperspectralLeafBiophysics(config);
        ρ₁,τ₁ = LO.layer_ρ_τ(lha, bio, 5.0, 1/bio.MESOPHYLL_N, 40.0);
        ρ₂,τ₂ = LO.layer_ρ_τ(lha, bio, 5.0, 1 - 1/bio.MESOPHYLL_N, 90.0);
        α₁_sife = (1 .- ρ₁ .- τ₁)[wls.IΛ_SIFE];
        α₂_sife = (1 .- ρ₂ .- τ₂)[wls.IΛ_SIFE];
        sife₁ = LO.layer_raw_sife(lha, wls, bio, 5.0, 1/bio.MESOPHYLL_N, 40.0);
        sife₂ = LO.layer_raw_sife(lha, wls, bio, 5.0, 1 - 1/bio.MESOPHYLL_N, 90.0);
        @test all(0 .< sife₁ .< 1);
        @test all(0 .< sife₂ .< 1);
        @test all(α₁_sife .≈ sife₁);
        @test all(α₂_sife .≈ sife₂);
    end;

    @testset "Raw SIF emission matrices of the leaf" begin
        config = NS.SPACConfiguration{Float64}(DATASET = NS.LAND_2021_1NM);
        lha = config.LHA;
        wls = config.WLSET;
        bio = NS.HyperspectralLeafBiophysics(config);
        ρ,τ = LO.leaf_spectra(lha, bio, 5.0, 40.0);
        α_sife = (1 .- ρ .- τ)[wls.IΛ_SIFE];
        _,f_cab,f_car = LO.sublayer_τ(lha, bio, 5.0, 1/bio.MESOPHYLL_N, 10);
        ϕ_sife = f_cab[wls.IΛ_SIFE];
        mat_b,mat_f = LO.leaf_raw_sif_matrices(lha, wls, bio, 5.0, 40.0);
        rad = ones(size(mat_b,1));
        sif_b = mat_b' * rad;
        sif_f = mat_f' * rad;
        @test all(0 .<= mat_b .< 1);
        @test all(0 .<= mat_f .< 1);
        @test sum(sif_b .+ sif_f) ≈ (ϕ_sife .* α_sife)' * rad;
    end;

    @testset "SIF emission matrices of the leaf" begin
        config = NS.SPACConfiguration{Float64}(DATASET = NS.LAND_2021_1NM);
        lha = config.LHA;
        wls = config.WLSET;
        bio = NS.HyperspectralLeafBiophysics(config);
        ρ,τ = LO.leaf_spectra(lha, bio, 5.0, 40.0);
        α_sife = (1 .- ρ .- τ)[wls.IΛ_SIFE];
        _,f_cab,f_car = LO.sublayer_τ(lha, bio, 5.0, 1/bio.MESOPHYLL_N, 10);
        ϕ_sife = f_cab[wls.IΛ_SIFE];
        mat_b,mat_f = LO.leaf_sif_matrices(lha, wls, bio, 5.0, 40.0);
        rad = ones(size(mat_b,1));
        sif_b = mat_b' * rad;
        sif_f = mat_f' * rad;
        @test all(0 .<= mat_b .< 1);
        @test all(0 .<= mat_f .< 1);
        @test sum(sif_b .+ sif_f) .< (ϕ_sife .* α_sife)' * rad;
    end;
end;
