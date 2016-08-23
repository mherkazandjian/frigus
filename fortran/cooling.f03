program cooling
     use types_and_parameters, only: it, energy_lev,                            &
                                     radiative_coeffs, collisional_coeffs,  &
                                     reaction_matrix, population

     use level_population, only: lev_pop, tests
    
    type(energy_lev) :: energy
    type(radiative_coeffs) :: a21, b21, r21, b12, r12, jnu
    type(collisional_coeffs) :: rr, rr21, rr12
    type(reaction_matrix)  :: coll_rad_matrix    
    type(population) :: x
    
     call lev_pop(energy, a21, b21, r21, b12, jnu, r12, rr, rr21, rr12, it, coll_rad_matrix, x)
      
     call tests(energy, rr, rr21, rr12, a21, b21, r21, b12, jnu, r12, coll_rad_matrix, it)

     
     
end program cooling