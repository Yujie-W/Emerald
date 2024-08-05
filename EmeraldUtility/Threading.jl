module Threading

using Distributed: addprocs, rmprocs, workers


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-11: add function to run the SPAC
#     2024-Feb-23: add exeflags option
#
#######################################################################################################################################################################################################
"""

    dynamic_workers!(nTH::Int; exeflags::String = "--project")

Add processors to run code in multiple threadings, given
- `threads` Number of threads
- `exeflags` Flags for the Julia executable

"""
function dynamic_workers!(threads::Int; exeflags::String = "--project")
    max_threads = Sys.CPU_THREADS;

    # if        no worker yet, and nTH <= MaxThreads, hire nTH
    # elseif    no worker yet, but nTH > MaxThreads, hire MaxThreads
    # elseif    some workers already, and nTH <= MaxThreads, hire more
    # elseif    some workers already, and nTH > MaxThreads, hire more
    # else      workers is more than expected, remove the extra
    if (length(workers()) == 1) && (workers()[1] == 1) && (threads <= max_threads)
        addprocs(threads; exeflags = exeflags);
    elseif (length(workers()) == 1) && (workers()[1] == 1) && (threads > max_threads)
        addprocs(max_threads; exeflags = exeflags);
    elseif length(workers()) < threads && (threads <= max_threads)
        addprocs(threads - length(workers()); exeflags = exeflags);
    elseif length(workers())<max_threads && (threads > max_threads)
        addprocs(max_threads-length(workers()); exeflags = exeflags);
    else
        to_remove = workers()[(threads + 1):end];
        rmprocs(to_remove...);
    end;

    return nothing
end;


end; # module
