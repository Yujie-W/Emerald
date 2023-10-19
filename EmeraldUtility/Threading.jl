module Threading

using Distributed: addprocs, rmprocs, workers


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2023-Mar-11: add function to run the SPAC
#
#######################################################################################################################################################################################################
"""

    dynamic_workers!(nTH::Int)

Add processors to run code in multiple threadings, given
- `threads` Number of threads

"""
function dynamic_workers!(threads::Int)
    _MaxThreads = Sys.CPU_THREADS;

    # if        no worker yet, and nTH <= MaxThreads, hire nTH
    # elseif    no worker yet, but nTH > MaxThreads, hire MaxThreads
    # elseif    some workers already, and nTH <= MaxThreads, hire more
    # elseif    some workers already, and nTH > MaxThreads, hire more
    # else      workers is more than expected, remove the extra
    if (length(workers()) == 1) && (workers()[1] == 1) && (threads <= _MaxThreads)
        addprocs(threads; exeflags="--project");
    elseif (length(workers()) == 1) && (workers()[1] == 1) && (threads > _MaxThreads)
        addprocs(_MaxThreads; exeflags="--project");
    elseif length(workers()) < threads && (threads <= _MaxThreads)
        addprocs(threads - length(workers()); exeflags="--project");
    elseif length(workers())<_MaxThreads && (threads > _MaxThreads)
        addprocs(_MaxThreads-length(workers()); exeflags = "--project");
    else
        _to_remove = workers()[(threads + 1):end];
        rmprocs(_to_remove...);
    end;

    return nothing
end;


end; # module
