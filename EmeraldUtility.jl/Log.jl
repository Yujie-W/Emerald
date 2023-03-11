module Log


#######################################################################################################################################################################################################
#
# Changes to these functions
# General
#     2022-Aug-24: move function outside of the folder
#     2023-Jan-19: make these functions to macros
#
#######################################################################################################################################################################################################
"""

    macro terror(exps...)

Add a time tag to @error expression, and display the message

"""
macro terror(exps...)
    quote
        @error "$(format(now(),"yyyy-mm-dd HH:MM:SS"))\n       $($(esc(exps[1])))" $(exps[2:end]...)
    end
end


"""

    macro tinfo(exps...)

Add a time tag to @info expression, and display the message

"""
macro tinfo(exps...)
    quote
        @info "$(format(now(),"yyyy-mm-dd HH:MM:SS"))\n      $($(esc(exps[1])))" $(exps[2:end]...)
    end
end


"""

    macro twarn(exps...)

Add a time tag to @warn expression, and display the message

"""
macro twarn(exps...)
    quote
        @warn "$(format(now(),"yyyy-mm-dd HH:MM:SS"))\n         $($(esc(exps[1])))" $(exps[2:end]...)
    end
end


end # module
