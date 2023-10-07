using Test
import Emerald.EmeraldLand.LeafOptics as LO
import Emerald.EmeraldLand.Namespace as NS


@testset verbose = true "Leaf Optics Model" begin
    @testset "Interface ρ and τ for direct radiation" begin
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

    @testset "Interface τ for isotropic radiation" begin
        n₁ = 1.0
        n₂ = 1.33
        for θ₁ in collect(Float64, 1:1:90)
            τ = LO.interface_isotropic_τ(n₁, n₂, θ₁);
            @test 0 <= τ <= 1;
            τ = LO.interface_isotropic_τ(n₂, n₁, θ₁);
            @test 0 <= τ <= 1;
        end;
    end;

    @testset "Leaf interface ρ and τ" begin
        config = NS.SPACConfiguration{Float64}(DATASET = NS.LAND_2021_1NM);
        bio = NS.LeafBio(config);
        LO.leaf_interface_ρ_τ!(config, bio, 40.0);

        @test all(0 .< bio.auxil.ρ_interface_θ  .< 1);
        @test all(0 .< bio.auxil.τ_interface_θ  .< 1);
        @test all(0 .< bio.auxil.ρ_interface_12 .< 1);
        @test all(0 .< bio.auxil.τ_interface_12 .< 1);
        @test all(0 .< bio.auxil.ρ_interface_21 .< 1);
        @test all(0 .< bio.auxil.τ_interface_21 .< 1);
    end;

    @testset "Leaf sublayer τ and f" begin
        config = NS.SPACConfiguration{Float64}(DATASET = NS.LAND_2021_1NM);
        bio = NS.LeafBio(config);
        LO.leaf_interface_ρ_τ!(config, bio, 40.0);
        LO.leaf_sublayer_f_τ!(config, bio, 5.0, 10);

        @test all(0 .<= bio.auxil.f_cab   .<= 1);
        @test all(0 .<= bio.auxil.f_car   .<= 1);
        @test all(0 .<= bio.auxil.τ_sub_1 .<= 1);
        @test all(0 .<= bio.auxil.τ_sub_2 .<= 1);
    end;

    @testset "Leaf layer ρ and τ" begin
        config = NS.SPACConfiguration{Float64}(DATASET = NS.LAND_2021_1NM);
        bio = NS.LeafBio(config);
        LO.leaf_interface_ρ_τ!(config, bio, 40.0);
        LO.leaf_sublayer_f_τ!(config, bio, 5.0, 10);
        LO.leaf_layer_ρ_τ!(bio, 10);

        @test all(0 .< bio.auxil.ρ_layer_θ .< 1);
        @test all(0 .< bio.auxil.τ_layer_θ .< 1);
        @test all(0 .< bio.auxil.ρ_layer_1 .< 1);
        @test all(0 .< bio.auxil.τ_layer_1 .< 1);
        @test all(0 .< bio.auxil.ρ_layer_2 .< 1);
        @test all(0 .< bio.auxil.τ_layer_2 .< 1);
        @test all(bio.auxil.ρ_layer_θ .+ bio.auxil.τ_layer_θ .< 1);
        @test all(bio.auxil.ρ_layer_1 .+ bio.auxil.τ_layer_1 .< 1);
        @test all(bio.auxil.ρ_layer_2 .+ bio.auxil.τ_layer_2 .< 1);
    end;

    @testset "Leaf layer effective interface ρ" begin
        config = NS.SPACConfiguration{Float64}(DATASET = NS.LAND_2021_1NM);
        bio = NS.LeafBio(config);
        LO.leaf_interface_ρ_τ!(config, bio, 40.0);
        LO.leaf_sublayer_f_τ!(config, bio, 5.0, 10);
        LO.leaf_layer_ρ_τ!(bio, 10);

        # test if the effective ρ_12 and ρ_21 are same as the modeled using the layer 1 data
        ρ_12 = LO.effective_ρ_12.(bio.auxil.ρ_layer_1, bio.auxil.τ_layer_1, bio.auxil.τ_sub_1.^10);
        ρ_21 = LO.effective_ρ_21.(bio.auxil.ρ_layer_1, bio.auxil.τ_layer_1, bio.auxil.τ_sub_1.^10);
        @test all(bio.auxil.ρ_interface_12 .≈ ρ_12);
        @test all(bio.auxil.ρ_interface_21 .≈ ρ_21);

        # test the effective ρ_12 and ρ_21 converge back to the layer 2 data
        ρ_2 = bio.auxil.ρ_layer_2;
        τ_2 = bio.auxil.τ_layer_2;
        τ_all = bio.auxil.τ_sub_2 .^ 10
        ρ_12 = LO.effective_ρ_12.(ρ_2, τ_2, τ_all);
        ρ_21 = LO.effective_ρ_21.(ρ_2, τ_2, τ_all);
        τ_re =                          (1 .- ρ_12) .* τ_all .* (1 .- ρ_21) ./ (1 .- τ_all .^ 2 .* ρ_21 .^ 2);
        ρ_re = ρ_12 .+ τ_all .* ρ_21 .* (1 .- ρ_12) .* τ_all .* (1 .- ρ_21) ./ (1 .- τ_all .^ 2 .* ρ_21 .^ 2);
        @test all(ρ_re .≈ ρ_2);
        @test all(τ_re .≈ τ_2);
    end;

    @testset "Leaf ρ and τ" begin
        config = NS.SPACConfiguration{Float64}(DATASET = NS.LAND_2021_1NM);
        bio = NS.LeafBio(config);

        LO.leaf_spectra!(config, bio, 5.0, 40.0; N = 10);
        @test all(0 .< bio.auxil.ρ_leaf .< 1);
        @test all(0 .< bio.auxil.τ_leaf .< 1);
        @test all(bio.auxil.ρ_leaf .+ bio.auxil.τ_leaf .< 1);

        LO.leaf_spectra!(config, bio, 5.0, 59.0; N = 10);
        @test all(0 .< bio.auxil.ρ_leaf .< 1);
        @test all(0 .< bio.auxil.τ_leaf .< 1);
        @test all(bio.auxil.ρ_leaf .+ bio.auxil.τ_leaf .< 1);
    end;

    @testset "Leaf SIF backward and forward matrices" begin
        config = NS.SPACConfiguration{Float64}(DATASET = NS.LAND_2021_1NM);
        bio = NS.LeafBio(config);
        LO.leaf_spectra!(config, bio, 5.0, 40.0; N = 10);
        rad = ones(Float64, size(bio.auxil.mat_b,2));
        sif_b = bio.auxil.mat_b * rad;
        sif_f = bio.auxil.mat_f * rad;

        @test all(0 .<= bio.auxil.mat_b .< 1);
        @test all(0 .<= bio.auxil.mat_f .< 1);
        @test sum(sif_b .+ sif_f) .< sum(bio.auxil.α_leaf[config.SPECTRA.IΛ_SIFE] .* bio.auxil.f_sife[config.SPECTRA.IΛ_SIFE] .* rad);
    end;
end;
