function build_mc_opf(pm::AbstractExplicitNeutralIVRModelVoltVar)
    # Variables
    PMD.variable_mc_bus_voltage(pm)
    PMD.variable_mc_branch_current(pm)
    PMD.variable_mc_load_current(pm) #creates dictionary position for load currents (Is this still needed if currents are defined as variables?)
    PMD.variable_mc_load_power(pm) #creates dictionary position for load powers (Is this still needed if powers are defined as variables?)
    PMD.variable_mc_generator_current(pm)
    PMD.variable_mc_generator_power(pm)
    PMD.variable_mc_transformer_current(pm)
    PMD.variable_mc_transformer_power(pm)
    PMD.variable_mc_switch_current(pm)
    variable_mc_load_power_VoltVar(pm) #Defines powers as variables for loads with VoltVAr control
    #variable_mc_load_current_VoltVar(pm) #Defines currents as variables for loads with VoltVAr control (if wanted)

    # Constraints
    for i in ids(pm, :bus)

        if i in ids(pm, :ref_buses)
            PMD.constraint_mc_voltage_reference(pm, i)
        end

        PMD.constraint_mc_voltage_absolute(pm, i)
        PMD.constraint_mc_voltage_pairwise(pm, i)
    end

    # components should be constrained before KCL, or the bus current variables might be undefined

    for id in ids(pm, :gen)
        PMD.constraint_mc_generator_power(pm, id)
        PMD.constraint_mc_generator_current(pm, id)
    end

    for id in ids(pm, :load) 
        PMD.constraint_mc_load_power(pm, id)
        constraint_mc_load_current(pm, id)
        #constraint_mc_load_current_vars(pm, id) #adds constraints for VoltVAr loads if the currents are defined as variables
    end

    for i in ids(pm, :transformer)
        PMD.constraint_mc_transformer_voltage(pm, i)
        PMD.constraint_mc_transformer_current(pm, i)
        PMD.constraint_mc_transformer_thermal_limit(pm, i)
    end

    for i in ids(pm, :branch)
        PMD.constraint_mc_current_from(pm, i)
        PMD.constraint_mc_current_to(pm, i)
        PMD.constraint_mc_bus_voltage_drop(pm, i)

        PMD.constraint_mc_branch_current_limit(pm, i)
        PMD.constraint_mc_thermal_limit_from(pm, i)
        PMD.constraint_mc_thermal_limit_to(pm, i)
    end

    for i in ids(pm, :switch)
        PMD.constraint_mc_switch_current(pm, i)
        PMD.constraint_mc_switch_state(pm, i)

        PMD.constraint_mc_switch_current_limit(pm, i)
        PMD.constraint_mc_switch_thermal_limit(pm, i)
    end

    for i in ids(pm, :bus)
        PMD.constraint_mc_current_balance(pm, i)
    end

    # Objective
    PMD.objective_mc_min_fuel_cost(pm)
end
