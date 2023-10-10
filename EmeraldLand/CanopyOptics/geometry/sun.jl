# This file contains function to compute the sun geometry related parameters (apart from sensor geometry)

function sun_geometry!(can::MultiLayerCanopy{FT}) where {FT}
    # extinction coefficients for the solar radiation
    for i in 1:9
        Cs = cosd(lia) * cosd(can.sun_geometry.state.sza);
        Ss = sind(lia) * sind(can.sun_geometry.state.sza);
        βs = (Cs >= Ss ? FT(π) : acos(-Cs/Ss));
        ks = 2 / FT(π) / cosd(can.sun_geometry.state.sza) * (Cs * (βs - FT(π)/2) + Ss * sin(βs));
    end;

    return nothing
end;
