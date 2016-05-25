program cooling
     use types_and_parameters, only: energy_lev,                          &
                                     radiative_coeffs,                    &
                                     population

    
     use level_population, only: lev_pop
    
    type(energy_lev) :: energy
    type(radiative_coeffs) :: a21
    type(population) :: x
    
     call lev_pop(energy, a21, x)
     
     
end program cooling