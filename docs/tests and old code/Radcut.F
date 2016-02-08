c a cut on the invariant mass is only necessary for DSUB_I in radiation phase


        ! cut on the onvariant mass sij
        logical doradcut
        double precision theta
        double precision delta_cut
        ! invariants
        double precision sij, mij, mi, mj, mk
        ! indices
        integer i,j,k

        ! ==- cut on invariant mass when in radiation-phase -==
        ! We investigated that there are many negative events generated
        ! close to the resonant region. We introduced a cut on the 
        ! invariant mass Sqrt(sij) of the intermediate squark to avoid
        ! negative amplitudes squared. We don't have to do this
        ! in the DR-I case.
        doradcut = flg_rad ! cut only in radiation-phase
#ifdef DR_I
        doradcut = .false.
#endif

        if(doradcut) then

#ifdef DEBUG
          print*,"in Reals_osres.f. TODO: check rad_cut"
          stop
#endif
        
          ! reset the temp. variables
          theta = 0D0
          delta_cut = 1D0

          call increasecnt('total calls to real-routine in rad-gen')
          
          ! TODO: setze ichan je nachdem wie rflav gesetzt wurde.
          call set_channel(osresID(ichan),i,j,k,mi,mj,mk,mij)
          sij = momsum2sq(p(0:3,i),p(0:3,j)) ! invariant mass of particle i + j

          ! we need to cut only in case of on-shell-resonances!
          if( (s12.ge.(mij+mk)**2) .and. (mij.ge.(mi+mj)) ) then
            theta = dabs(mij-dsqrt(sij))-delta_cut
          endif  
        
          ! theta function
          if((theta.le.0D0) .and. (ampreal.le.0.d0)) then
            ampreal = 1D-20
            call increasecnt('radiation in on-shell-res. region vetoed')
          endif
        endif  
        ! ==- end cut on invariant mass when in radiation-phase -==