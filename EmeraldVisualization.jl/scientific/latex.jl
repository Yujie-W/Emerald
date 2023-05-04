#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Dec-08: refactor the functions that set text font and render
#
#######################################################################################################################################################################################################
"""

    latex_symbol(
                mid::String;
                option::String = "mathrm",
                presub::Union{Nothing, String} = nothing,
                presup::Union{Nothing, String} = nothing,
                sub::Union{Nothing, String} = nothing,
                sup::Union{Nothing, String} = nothing)

Return a string for symbol with LaTeX syntax, given
- `mid` Symbol in the middle
- `option` Option to remove text format, such as `\\text` or `\\mathrm`
- `presub` Subscript before `mid`
- `presup` Superscript before `mid`
- `sub` Subscript after `mid`
- `sup` Superscript after `mid`

"""
function latex_symbol(
            mid::String;
            option::String = "mathrm",
            presub::Union{Nothing, String} = nothing,
            presup::Union{Nothing, String} = nothing,
            sub::Union{Nothing, String} = nothing,
            sup::Union{Nothing, String} = nothing)
    @assert option in ["mathrm", "text"] "Option `$(option)` not supported!";

    _str = "\$";
    _str *= (isnothing(presub) || length(presub) == 0) ? "" : "_\\$(option){$(presub)}";
    _str *= (isnothing(presup) || length(presup) == 0) ? "" : "^\\$(option){$(presup)}";
    _str *= length(mid) <= 1 ? mid : "\\$(option){$(mid)}";
    _str *= (isnothing(sub) || length(sub) == 0) ? "" : "_\\$(option){$(sub)}";
    _str *= (isnothing(sup) || length(sup) == 0) ? "" : "^\\$(option){$(sup)}";
    _str *= "\$";

    return _str
end


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Dec-08: refactor the functions that set text font and render
#
#######################################################################################################################################################################################################
"""

    text_render!(; font::String = "sans", use_latex::Bool = false)

Set text font and render for PyPlot, given
- `font` Font family (either sans or serif)
- `use_latex` If true, use LaTeX to render the text

"""
function set_text_render!(; font::String = "sans", use_latex::Bool = false)
    # switch between sans and serif
    @assert font in ["sans", "serif"] "Font must be sans or serif!";
    if font == "sans"
        rc("font", family="sans-serif", serif=["sans-serif"]);
        rc("mathtext", fontset="dejavusans");
    else
        rc("font", family="serif", serif=["Palatino"]);
        rc("mathtext", fontset="dejavuserif");
    end;

    # switch between LaTeX and normal renders
    rc("text", usetex = use_latex);
    if use_latex
        if font == "sans"
            rc("text.latex", preamble = "\\usepackage{amsmath,DejaVuSans,sansmath,upgreek} \\sansmath");
        else
            rc("text.latex", preamble = "\\usepackage{amsmath,newpxmath,newpxtext}");
        end;
    end;

    return nothing
end
