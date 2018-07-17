subroutine check_and_save(SNPP)
#include "cpp_options.h"
    use global, only :tt,saveFreq,diagFreq,pickupFreq,&
        rec_num,dt_reinit,iswitch,pickup
    implicit none
    INTEGER*8 :: t0,t1,IPP
    INTEGER*8, intent(in) :: SNPP

    !$OMP PARALLEL SECTIONS

    !$OMP SECTION
    if (mod(tt,saveFreq) .eq. 0.0) then
        t0=abs(iswitch-1)
        t1=iswitch
        print*, "write data to files at step ",rec_num,' tt=',tt
        do IPP=1,SNPP

#ifdef saveTSG
            call load_tsg(rec_num,iswitch)
            call interp_tracer(t0,t1,IPP)
#endif

#ifdef saveGradient
            call calc_gradient(t0,IPP)
#endif
            call save_data(IPP)
        enddo

    endif

    !$OMP SECTION
    !call count_stagnant()
#ifdef monitoring
    if (mod(tt,diagFreq) .eq. 0) then
        call diag()
    endif
#endif

!    !$OMP SECTION
!    if (mod(tt,pickupFreq) .eq. 0) then
!        print*, "==========================================="
!        print*, " Dump pickup data at record ", rec_num
!        print*, "==========================================="
!        call save_pickup()
!    endif

    !$OMP SECTION
    !reinitialize particles if invoked
    if (dt_reinit>0 .and. mod(tt,dt_reinit) .eq. 0) then
        pickup=0d0
        call init_particles()
    endif
!$OMP END PARALLEL SECTIONS
end subroutine check_and_save
