program cooling
     use types_and_parameters, only: energy_lev,                            &
                                     radiative_coeffs, collisional_coeffs,  &
                                     population

    
     use level_population, only: lev_pop, tests
    
    type(energy_lev) :: energy
    type(radiative_coeffs) :: a21
    type(collisional_coeffs) :: rr    
    type(population) :: x
    
     call lev_pop(energy, a21, rr, x)
     
     call tests(rr)
     
     
end program cooling