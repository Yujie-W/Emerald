""" General structure to map parameter saving settings to functions """
mutable struct ParameterFunctionMapper
    "Parameter name"
    name::String
    "Whether to save the parameter"
    to_save::Bool
    "Function to compute the parameter"
    func::Function
    "Extra parameters to be passed to the function other than config and spac"
    params::Vector
end;
