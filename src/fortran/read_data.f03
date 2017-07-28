module read_data

    use types_and_parameters, only: nlev, energy_lev, vmax, jmax,         &
                                    ntrans, ntemp, i,                     &
                                    radiative_coeffs, collisional_coeffs, &
                                    reaction_matrix,                      &
                                    population,                           &
                                    Trad, ini, fin,                       &
                                    nlev_lique, vi, vf, ji, jf,           &
                                    id_temp, id_temp_test, ilique_flag,   &
                                    iflower_flag, ilipovka_flag
                                    
    use energy_levels, only: reading_data_energies_stancil,               &
                             reading_data_energies_lique,                 &
                             reading_data_energies_flower,                &
                             reading_data_energies_lipovka
    
    
    use radiation,     only: reading_data_radiative,                      &
                             radiative_downwards, radiative_upwards,      &
                             reading_data_radiative_lipovka

    use collisions,    only: reading_data_collisions_lique,                 &
                             reading_data_collisions_flower,                &
                             reading_data_collisions_lipovka

    use testing_data, only: tests, writing_files

    contains

    subroutine get_data(ilique_flag, iflower_flag, ilipovka_flag, Trad, &
    id_temp, energy, a21, b21, r21, b12, jnu, r12, rr, rr21, rr12)
         type(energy_lev) :: energy
         type(radiative_coeffs) :: a21, b21, b12, jnu, r21, r12, rad
         type(collisional_coeffs) :: rr, rr21, rr12
         type(reaction_matrix) :: coll_rad_matrix
         real*8 :: Trad

         if(ilique_flag.eq.1) then
            call reading_data_energies_lique(energy)
         elseif(iflower_flag.eq.1) then
            call reading_data_energies_flower(energy)
         elseif(ilipovka_flag.eq.1) then
            call reading_data_energies_lipovka(energy)
         else 
            print*, 'Error raised: no flag selected for calling energy data reader'
         endif

        if(ilique_flag.eq.1.or.iflower_flag.eq.1) then
            call reading_data_radiative(energy, a21)
         elseif(ilipovka_flag.eq.1) then
            call reading_data_radiative_lipovka(energy, a21)
         else 
            print*,  'Error raised: no flag selected for calling radiative data reader'
        endif


         if(ilique_flag.eq.1) then
            call reading_data_collisions_lique(energy, rr, rr21, rr12)
         elseif(iflower_flag.eq.1) then
            call reading_data_collisions_flower(energy, rr, rr21, rr12)
         elseif(ilipovka_flag.eq.1) then
            call reading_data_collisions_lipovka(energy, rr, rr21, rr12)
         else
            print*, 'Error raised: no flag selected for calling collisional data reader'
         endif
         
         call radiative_downwards(energy, Trad, a21, b21, r21)

         call radiative_upwards(energy, Trad, b21, b12, jnu, r12)

         call tests(energy, rr, rr21, rr12, a21, b21, r21, b12, jnu, r12, coll_rad_matrix, id_temp)

         call writing_files(a21, b21, b12, jnu, id_temp, rr, coll_rad_matrix)

        return
     end subroutine get_data

end module read_data