
function energy = advance_E_compiled2(energy, d_energy, timestep)
            
energy = energy + d_energy .* timestep;