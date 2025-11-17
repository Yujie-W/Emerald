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

Users may choose to use `@terror_pre`, `@terror_mid`, and `@terror_end` to display multi-line error messages
```julia
begin
    @terror_pre "Timed error message:";
    @terror_mid "This is the second line of the error message.";
    @terror_end "This is the last line of the error message.";
end;
```
"""
macro terror(exps...)
    quote
        @error "$(format(now(),"yyyy-mm-dd HH:MM:SS"))\n       $($(esc(exps[1])))" $(exps[2:end]...)
    end;
end;

macro terror_pre(msg::String)
    printstyled("┌ Error: "; bold = true, color = :red);
    printstyled("$(format(now(),"yyyy-mm-dd HH:MM:SS"))\n"; bold = false, color = :white);
    printstyled("│        "; bold = true, color = :red);
    printstyled(msg; bold = false, color = :white);
    print("\n");
end;

macro terror_mid(msg::String)
    printstyled("│        "; bold = true, color = :red);
    printstyled(msg; bold = false, color = :white);
    print("\n");
end;

macro terror_end(msg::String)
    printstyled("└        "; bold = true, color = :red);
    printstyled(msg; bold = false, color = :white);
    print("\n");
end;


"""

    macro tinfo(exps...)

Add a time tag to @info expression, and display the message

Users may choose to use `@tinfo_pre`, `@tinfo_mid`, and `@tinfo_end` to display multi-line info messages
```julia
begin
    @tinfo_pre "Timed info message:";
    @tinfo_mid "This is the second line of the info message.";
    @tinfo_end "This is the last line of the info message.";
end;
```

"""
macro tinfo(exps...)
    quote
        @info "$(format(now(),"yyyy-mm-dd HH:MM:SS"))\n      $($(esc(exps[1])))" $(exps[2:end]...)
    end;
end;

macro tinfo_pre(msg::String)
    printstyled("┌ Info: "; bold = true, color = :cyan);
    printstyled("$(format(now(),"yyyy-mm-dd HH:MM:SS"))\n"; bold = false, color = :white);
    printstyled("│       "; bold = true, color = :cyan);
    printstyled(msg; bold = false, color = :white);
    print("\n");
end;

macro tinfo_mid(msg::String)
    printstyled("│       "; bold = true, color = :cyan);
    printstyled(msg; bold = false, color = :white);
    print("\n");
end;

macro tinfo_end(msg::String)
    printstyled("└       "; bold = true, color = :cyan);
    printstyled(msg; bold = false, color = :white);
    print("\n");
end;


"""

    macro twarn(exps...)

Add a time tag to @warn expression, and display the message

Users may choose to use `@twarn_pre`, `@twarn_mid`, and `@twarn_end` to display multi-line warning messages
```julia
begin
    @twarn_pre "Timed warning message:";
    @twarn_mid "This is the second line of the warning message.";
    @twarn_end "This is the last line of the warning message.";
end;
```
"""
macro twarn(exps...)
    quote
        @warn "$(format(now(),"yyyy-mm-dd HH:MM:SS"))\n         $($(esc(exps[1])))" $(exps[2:end]...)
    end;
end;

macro twarn_pre(msg::String)
    printstyled("┌ Warning: "; bold = true, color = :yellow);
    printstyled("$(format(now(),"yyyy-mm-dd HH:MM:SS"))\n"; bold = false, color = :white);
    printstyled("│          "; bold = true, color = :yellow);
    printstyled(msg; bold = false, color = :white);
    print("\n");
end;

macro twarn_mid(msg::String)
    printstyled("│          "; bold = true, color = :yellow);
    printstyled(msg; bold = false, color = :white);
    print("\n");
end;

macro twarn_end(msg::String)
    printstyled("└          "; bold = true, color = :yellow);
    printstyled(msg; bold = false, color = :white);
    print("\n");
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
        @tinfo_pre msg;
    elseif msg_level == "tinfo_mid"
        @tinfo_mid msg;
    elseif msg_level == "tinfo_end"
        @tinfo_end msg;

    elseif msg_level == "warn"
        @warn msg;
    elseif msg_level == "twarn"
        @twarn msg;
    elseif msg_level == "twarn_pre"
        @twarn_pre msg;
    elseif msg_level == "twarn_mid"
        @twarn_mid msg;
    elseif msg_level == "twarn_end"
        @twarn_end msg;

    elseif msg_level == "error"
        @error msg;
    elseif msg_level == "terror"
        @terror msg;
    elseif msg_level == "terror_pre"
        @terror_pre msg;
    elseif msg_level == "terror_mid"
        @terror_mid msg;
    elseif msg_level == "terror_end"
        @terror_end msg;

    elseif msg_level == "println"
        println(msg);
    elseif msg_level == "print"
        print(msg);
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
