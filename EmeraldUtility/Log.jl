module Log

using Dates: format, now


#######################################################################################################################################################################################################
#
# Changes to these functions
# General
#     2022-Aug-24: move function outside of the folder
#     2023-Jan-19: make these functions to macros
#     2025-Nov-17: add support for multi-line messages using *_pre, *_mid, and *_end macros
#
#######################################################################################################################################################################################################
"""

    macro terror(exps...)

Add a time tag to @error expression, and display the message.

"""
macro terror(exps...)
    quote
        @error "$(format(now(),"yyyy-mm-dd HH:MM:SS"))\n       $($(esc(exps[1])))" $(exps[2:end]...)
    end;
end;


"""

    macro tinfo(exps...)

Add a time tag to @info expression, and display the message

"""
macro tinfo(exps...)
    quote
        @info "$(format(now(),"yyyy-mm-dd HH:MM:SS"))\n      $($(esc(exps[1])))" $(exps[2:end]...)
    end;
end;


"""

    macro twarn(exps...)

Add a time tag to @warn expression, and display the message

"""
macro twarn(exps...)
    quote
        @warn "$(format(now(),"yyyy-mm-dd HH:MM:SS"))\n         $($(esc(exps[1])))" $(exps[2:end]...)
    end;
end;


##########################################################################################################################################################################################################
#
# Changes to this function
# General
#     2025-Nov-17: add a general call to display messages based on message level
#
##########################################################################################################################################################################################################
"""

    display_message!(msg::String, msg_level::String = "println")

Display the message based on the message level, given
- `msg` message to display
- `msg_level` message level for displaying

"""
function display_message!(msg::String, msg_level::String = "println")
    if msg_level == "info"
        @info msg;
    elseif msg_level == "tinfo"
        @tinfo msg;
    elseif msg_level == "tinfo_pre"
        display_timed_info!(msg, "pre");
    elseif msg_level == "tinfo_mid"
        display_timed_info!(msg, "mid");
    elseif msg_level == "tinfo_end"
        display_timed_info!(msg, "end");

    elseif msg_level == "warn"
        @warn msg;
    elseif msg_level == "twarn"
        @twarn msg;
    elseif msg_level == "twarn_pre"
        display_timed_warning!(msg, "pre");
    elseif msg_level == "twarn_mid"
        display_timed_warning!(msg, "mid");
    elseif msg_level == "twarn_end"
        display_timed_warning!(msg, "end");

    elseif msg_level == "error"
        @error msg;
    elseif msg_level == "terror"
        @terror msg;
    elseif msg_level == "terror_pre"
        display_timed_error!(msg, "pre");
    elseif msg_level == "terror_mid"
        display_timed_error!(msg, "mid");
    elseif msg_level == "terror_end"
        display_timed_error!(msg, "end");

    elseif msg_level == "println"
        println(msg);
    elseif msg_level == "print"
        print(msg);
    end;

    return nothing
end;

function display_timed_info!(msg::String, msg_level::String)
    if msg_level == "pre"
        printstyled("┌ Info: "; bold = true, color = :cyan);
        printstyled("$(format(now(),"yyyy-mm-dd HH:MM:SS"))\n"; bold = false, color = :white);
        printstyled("│       "; bold = true, color = :cyan);
        printstyled(msg; bold = false, color = :white);
        print("\n");
    end;

    if msg_level == "mid"
        printstyled("│       "; bold = true, color = :cyan);
        printstyled(msg; bold = false, color = :white);
        print("\n");
    end;

    if msg_level == "end"
        printstyled("└       "; bold = true, color = :cyan);
        printstyled(msg; bold = false, color = :white);
        print("\n");
    end;

    return nothing
end;

function display_timed_warning!(msg::String, msg_level::String)
    if msg_level == "pre"
        printstyled("┌ Warning: "; bold = true, color = :yellow);
        printstyled("$(format(now(),"yyyy-mm-dd HH:MM:SS"))\n"; bold = false, color = :white);
        printstyled("│          "; bold = true, color = :yellow);
        printstyled(msg; bold = false, color = :white);
        print("\n");
    end;

    if msg_level == "mid"
        printstyled("│          "; bold = true, color = :yellow);
        printstyled(msg; bold = false, color = :white);
        print("\n");
    end;

    if msg_level == "end"
        printstyled("└          "; bold = true, color = :yellow);
        printstyled(msg; bold = false, color = :white);
        print("\n");
    end;

    return nothing
end;

function display_timed_error!(msg::String, msg_level::String)
    if msg_level == "pre"
        printstyled("┌ Error: "; bold = true, color = :red);
        printstyled("$(format(now(),"yyyy-mm-dd HH:MM:SS"))\n"; bold = false, color = :white);
        printstyled("│        "; bold = true, color = :red);
        printstyled(msg; bold = false, color = :white);
        print("\n");
    end;

    if msg_level == "mid"
        printstyled("│        "; bold = true, color = :red);
        printstyled(msg; bold = false, color = :white);
        print("\n");
    end;

    if msg_level == "end"
        printstyled("└        "; bold = true, color = :red);
        printstyled(msg; bold = false, color = :white);
        print("\n");
    end;

    return nothing
end;


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Aug-24: move function outside of the folder
#     2022-Aug-24: add support to the case if the value is a vector but not a vector of pairs
#
#######################################################################################################################################################################################################
"""

    pretty_display!(pvec::Union{Vector{Pair{String,String}}, Vector{Pair{String,Any}}, Vector{Pair{Any,String}}, Vector{Pair{Any,Any}}}, spaces::String = "    ")

Display the pairs in a pretty way, given
- `pvec` Vector of pairs to display
- `spaces` Leading spaces before displaying the pair key

---
Examples
```julia
_pairs = ["A" => "b", "d" => "A", "rr" => ["ra" => "rB", "rD" => "ra"]];
pretty_display!(_pairs);
pretty_display!(_pairs, "  ");
```

"""
function pretty_display! end;

pretty_display!(pair::Pair, max_len::Int, spaces = "    ") = (
    # print leading spaces
    print(spaces);

    # print the key
    printstyled(string(pair[1]); color = :light_magenta);

    # print spaces after key and arrow
    print(repeat(" ", max_len - length(string(pair[1]))) * " ⇨ ");

    # if the value is a vector of pairs, recursive display
    if typeof(pair[2]) <: Vector && typeof(pair[2][1]) <: Pair
        # display [ and line break
        print("[\n");
        pretty_display!(pair[2], spaces * repeat(" ", max_len + 5));

        # display a ] and the next line
        print(spaces * repeat(" ", max_len + 4) * "],\n");

        return nothing;
    end;

    # if the value is not array of pair, display it
    printstyled(string(pair[2]); color = :cyan);
    print(",\n");

    return nothing;
);

pretty_display!(pvec::Union{Vector{Pair{String,String}}, Vector{Pair{String,Any}}, Vector{Pair{Any,String}}, Vector{Pair{Any,Any}}}, spaces::String = "    ") = (
    # determine the length of the keys
    _max_len = maximum(length.([string(_p[1]) for _p in pvec]));

    # display the elements
    pretty_display!.(pvec, _max_len, spaces);

    return nothing
);


end; # module
