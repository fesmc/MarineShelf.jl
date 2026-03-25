# Physical constants for marine shelf / ice-ocean interaction.
#
# All values are Float64 literals; call `T(ρ_ice)` etc. at use-sites
# to convert to the working precision T.

"""Density of ice [kg/m³]"""
const ρ_ice = 918.0

"""Density of freshwater [kg/m³]"""
const ρ_w = 1000.0

"""Density of seawater [kg/m³]"""
const ρ_sw = 1028.0

"""Gravitational acceleration [m/s²]"""
const g = 9.81

"""Specific heat capacity of the ocean mixed layer [J/(kg·K)]"""
const cp_o = 3974.0

"""Latent heat of fusion of ice [J/kg]"""
const L_ice = 3.34e5

"""Reference temperature (0 °C in Kelvin) [K]"""
const T0 = 273.15

"""Ice-to-seawater density ratio [-]"""
const ρ_ratio = ρ_ice / ρ_sw

"""
    ω  [1/K]

Combined physical constant appearing in the linear and quadratic BMB formulae
(Favier et al. 2019):

    ω = (ρ_sw · cp_o) / (ρ_ice · L_ice)
"""
const ω = (ρ_sw * cp_o) / (ρ_ice * L_ice)
