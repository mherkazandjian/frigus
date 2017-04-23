module read_data

    use types_and_parameters, only: nlev, energy_lev, vmax, jmax,         &
                                    ntrans, ntemp, i,                     &
                                    radiative_coeffs, collisional_coeffs, &
                                    reaction_matrix,                      &
                                    population,                           &
                                    Trad, ini, fin,                       &
                                    nlev_lique, vi, vf, ji, jf,           &
                                    id_temp, id_temp_test
                                    
    use energy_levels, only: reading_data_energies
    use radiation,     only: reading_data_radiative,                      &
                             radiative_downwards, radiative_upwards,      &
                             reading_data_radiative_lipovka
    use collisions,    only: reading_data_collisions

    use testing_data, only: tests, writing_files

    contains

    subroutine get_data(Trad, id_temp, energy, a21, b21, r21, b12, jnu, r12, rr, rr21, rr12)
         type(energy_lev) :: energy
         type(radiative_coeffs) :: a21, b21, b12, jnu, r21, r12, rad
         type(collisional_coeffs) :: rr, rr21, rr12
         type(reaction_matrix) :: coll_rad_matrix
         real*8 :: Trad

         call reading_data_energies(energy)

         !call reading_data_radiative(energy, a21)

         call reading_data_radiative_lipovka(energy, a21)         

         call reading_data_collisions(energy, rr, rr21, rr12)

         call radiative_downwards(energy, Trad, a21, b21, r21)

         call radiative_upwards(energy, Trad, b21, b12, jnu, r12)

         call tests(energy, rr, rr21, rr12, a21, b21, r21, b12, jnu, r12, coll_rad_matrix, id_temp)

         !call writing_files(a21, b21, b12, jnu, id_temp, rr, coll_rad_matrix)

        return
     end subroutine get_data

end module read_data