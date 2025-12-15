module MyVVPackage
    import JuMP
    import PowerModelsDistribution
    import InfrastructureModels
    const IM = InfrastructureModels
    const PMD = PowerModelsDistribution

    abstract type AbstractExplicitNeutralIVRModelVoltVar <: PMD.AbstractExplicitNeutralIVRModel end

    include("opf_vv.jl")
    include("Data.jl")
    include("VoltVarconstraints.jl")
end 
