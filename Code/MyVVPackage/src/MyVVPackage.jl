module MyVVPackage
    import JuMP
    import PowerModelsDistribution
    import PowerModelsDistribution: @pmd_fields, nw_id_default, pmd_it_sym, WYE
    import InfrastructureModels
    const IM = InfrastructureModels
    const PMD = PowerModelsDistribution

    abstract type AbstractExplicitNeutralIVRModelVoltVar <: PMD.AbstractNLExplicitNeutralIVRModel end
    mutable struct IVRENPowerModelVoltVar <: AbstractExplicitNeutralIVRModelVoltVar @pmd_fields end

    include("Data.jl")
    include("VoltVarconstraints.jl")
    include("opf_vv.jl")
end 
