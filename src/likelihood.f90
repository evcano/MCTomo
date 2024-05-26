! module of likelihood
! 
module m_likelihood

    use m_exception,       only : exception_raiseError
    use m_utils,           only : ii10, vs2vp_3d, vp2rho_3d
    use m_settings,        only : T_MOD
    use run_info,          only : T_RUN_INFO
    use like_settings,     only : T_DATA, T_LIKE_BASE, T_LIKE_SET, likeBaseSetup
    use cgal_delaunay,     only : d3
    use m_body_likelihood, only : body_likelihood, body_likelihood_grads, body_noise_likelihood, body_location_likelihood
    use m_surf_likelihood, only : surf_likelihood, surf_noise_likelihood, surf_likelihood_grads, surf_noise_likelihood_2

    implicit none

    private

    public :: T_LIKE, like_setup
    public :: likelihood, noise_likelihood, location_likelihood
    public :: update_lgP_grads
    
    ! define a likelihood type based on the base type
    type T_LIKE
        real(kind=ii10) :: like, misfit, unweighted_misfit
        real(kind=ii10), dimension(:), allocatable :: grads
        type(T_LIKE_BASE), dimension(:), allocatable :: likelihoods
    end type

contains

    subroutine like_setup(like, dat, like_set, ncell_max)
        implicit none
        type(T_LIKE), intent(out) :: like
        type(T_DATA), dimension(:), intent(in) :: dat
        type(T_LIKE_SET), intent(in) :: like_set
        integer, intent(in) :: ncell_max

        integer i
        
        like%like = huge(like%like)
        like%misfit = huge(like%misfit)
        like%unweighted_misfit = huge(like%unweighted_misfit)
        
        allocate(like%grads(2*ncell_max+4*dat(1)%nsrc))
        ! evcano:
        !  if datatype is 4, the dimension of likelihoods is 2.
        !  likelihoods(1) is for rgrp and rpha information
        !  likelihoods(2) is for lgrp  and lpha information
        if (like_set%datatype==3 .or. like_set%datatype==4) then
            allocate(like%likelihoods(2))
        else
            allocate(like%likelihoods(1))
        endif

        ! evcano:
        !  If datatype is 4, dat(1) holds rgrp data, dat(2) rpha data, dat(3) lgrp data, dat(4) lpha data.
        !  We setup likelihoods(1) using rgrp data and likelihoods(2) using lgrp data.
        !  likeBaseSetup does not use the data itself, it only uses the number of sources, receivers,
        !  frequencies and modes of each dataset.
        !  rgrp and rpha share the same number of sources, receivers, etc. The same goes for lgrp and lpha.
        !  Thus, it doesnt matter if we use rgrp or rpha to initialize likelihoods(1). The same applies for
        !  likelihoods(2).
        if (like_set%datatype==4) then
            call likeBaseSetup(like%likelihoods(1),dat(1),like_set,ncell_max)
            call likeBaseSetup(like%likelihoods(2),dat(3),like_set,ncell_max)
        else
            do i = 1, size(dat)
                call likeBaseSetup(like%likelihoods(i),dat(i),like_set,ncell_max)
            enddo
        endif

    end subroutine like_setup

    subroutine likelihood(dat,model,RTI,perturbed_box,like_set,like)
        implicit none
        type(T_DATA), dimension(:), intent(in)      :: dat
        type(T_MOD), intent(inout)                  :: model
        type(T_RUN_INFO), intent(inout)             :: RTI
        type(d3), dimension(2), intent(in)          :: perturbed_box
        type(T_LIKE_SET), intent(in)                :: like_set
        type(T_LIKE), intent(inout)                 :: like

        select case (like_set%datatype)
        case (0,1)
            call body_likelihood(dat(1),model,RTI,perturbed_box,like_set,like%likelihoods(1))
            like%like = like%likelihoods(1)%like
            like%misfit = like%likelihoods(1)%misfit
            like%unweighted_misfit = like%likelihoods(1)%unweighted_misfit
        case (2)
            ! only surface waves are used, set the vp according to vs and
            ! density according to vp
            call vs2vp_3d(model%vs,model%vp)
            call vp2rho_3d(model%vp,model%rho)
            call surf_likelihood(dat(1),model,RTI,perturbed_box,like_set,like%likelihoods(1))
            like%like = like%likelihoods(1)%like
            like%misfit = like%likelihoods(1)%misfit
            like%unweighted_misfit = like%likelihoods(1)%unweighted_misfit
        case (3)
            call body_likelihood(dat(1),model,RTI,perturbed_box,like_set,like%likelihoods(1))
            call surf_likelihood(dat(2),model,RTI,perturbed_box,like_set,like%likelihoods(2))
            like%like = like%likelihoods(1)%like + like%likelihoods(2)%like
            like%misfit = like%likelihoods(1)%misfit + like%likelihoods(2)%misfit
            like%unweighted_misfit = like%likelihoods(1)%unweighted_misfit + like%likelihoods(2)%unweighted_misfit
        ! evcano: TODO add case 4
        ! probably good to create new surf_likelihood function that evaluates group and phase measurements simultaneously
        ! and then call it twice, one for love and other for rayleigh
        ! thus, i need to modify like%likelihoods so it holds two, one for rayleigh (group+phase) and other for love(group+phase)
        end select

        !write(*,*) '-loglikelihood: ', like%like
        !write(*,*) 'misfits: ', like%misfit
        !write(*,*) 'unweighted misfits: ', like%unweighted_misfit

    end subroutine likelihood

    subroutine update_lgP_grads(dat,RTI,like_set,like)
        implicit none
        type(T_DATA), dimension(:), intent(in) :: dat
        type(T_RUN_INFO), intent(inout) :: RTI
        type(T_LIKE_SET), intent(in) :: like_set
        type(T_LIKE), intent(inout) :: like

        select case (like_set%datatype)
        case (0,1)
            call body_likelihood_grads(dat(1),RTI,like_set,like%likelihoods(1))
            like%like = like%likelihoods(1)%like
            like%misfit = like%likelihoods(1)%misfit
            like%unweighted_misfit = like%likelihoods(1)%unweighted_misfit
            like%grads = like%likelihoods(1)%grads
        case (2)
            ! only surface waves are used, set the vp according to vs and
            ! density according to vp
            !call vs2vp_3d(model%vs,model%vp)
            !call vp2rho_3d(model%vp,model%rho)
            call surf_likelihood_grads(dat(1),RTI,like_set,like%likelihoods(1))
            like%like = like%likelihoods(1)%like
            like%misfit = like%likelihoods(1)%misfit
            like%unweighted_misfit = like%likelihoods(1)%unweighted_misfit
            like%grads(1:RTI%ncells) = like%likelihoods(1)%grads(1:RTI%ncells)
        case (3)
            call body_likelihood_grads(dat(1),RTI,like_set,like%likelihoods(1))
            call surf_likelihood_grads(dat(2),RTI,like_set,like%likelihoods(2))
            like%like = like%likelihoods(1)%like + like%likelihoods(2)%like
            like%misfit = like%likelihoods(1)%misfit + like%likelihoods(2)%misfit
            like%unweighted_misfit = like%likelihoods(1)%unweighted_misfit + like%likelihoods(2)%unweighted_misfit
            like%grads = like%likelihoods(1)%grads
            like%grads(RTI%ncells+1:2*RTI%ncells) = like%likelihoods(1)%grads(RTI%ncells+1:2*RTI%ncells) + like%likelihoods(2)%grads(1:RTI%ncells)
        ! evcano: it seems this subroutine is only called  if HMC is used
        case (4)
            call exception_raiseError('update_lgP_grads not implemented for LIKE_SET%DATATYPE=4')
        end select

    endsubroutine update_lgP_grads

    subroutine noise_likelihood(dat,RTI,like_set,like)
        implicit none
        type(T_DATA), dimension(:), intent(in)      :: dat
        type(T_RUN_INFO), intent(inout)             :: RTI
        type(T_LIKE_SET), intent(in)                :: like_set
        type(T_LIKE), intent(inout)                 :: like

        select case (like_set%datatype)
        case (0,1)
            call body_noise_likelihood(dat(1),RTI,like_set,like%likelihoods(1))
            like%like = like%likelihoods(1)%like
            like%misfit = like%likelihoods(1)%misfit
            like%unweighted_misfit = like%likelihoods(1)%unweighted_misfit
        case (2)
            call surf_noise_likelihood(dat(1),RTI,like%likelihoods(1))
            like%like = like%likelihoods(1)%like
            like%misfit = like%likelihoods(1)%misfit
            like%unweighted_misfit = like%likelihoods(1)%unweighted_misfit
        case (3)
            call body_noise_likelihood(dat(1),RTI,like_set,like%likelihoods(1))
            call surf_noise_likelihood(dat(2),RTI,like%likelihoods(2))
            like%like = like%likelihoods(1)%like + like%likelihoods(2)%like
            like%misfit = like%likelihoods(1)%misfit + like%likelihoods(2)%misfit
            like%unweighted_misfit = like%likelihoods(1)%unweighted_misfit + like%likelihoods(2)%unweighted_misfit
        ! evcano: update new noise parameters
        case (4)
            ! noise likelihood of rayleigh-group and rayleigh-phase
            call surf_noise_likelihood_2(dat(1),dat(2),RTI,like%likelihoods(1),1)
            ! noise likelihood of love-group and love-phase
            call surf_noise_likelihood_2(dat(3),dat(4),RTI,like%likelihoods(2),2)

            like%like = like%likelihoods(1)%likeGroup + like%likelihoods(1)%likePhase + &
                like%likelihoods(2)%likeGroup + like%likelihoods(2)%likePhase

            like%misfit = like%likelihoods(1)%misfitGroup + like%likelihoods(1)%misfitPhase + &
                like%likelihoods(2)%misfitGroup + like%likelihoods(2)%misfitPhase

            like%unweighted_misfit = like%likelihoods(1)%unweighted_misfitGroup + like%likelihoods(1)%unweighted_misfitPhase + &
                like%likelihoods(2)%unweighted_misfitGroup + like%likelihoods(2)%unweighted_misfitPhase
        end select

        !write(*,*) '-loglikelihood: ', like%like
        !write(*,*) 'misfits: ', like%misfit
        !write(*,*) 'unweighted misfits: ', like%unweighted_misfit
    end subroutine noise_likelihood

    
    subroutine location_likelihood(dat,model,RTI,perturbed_box,like_set,like)
        implicit none
        type(T_DATA), dimension(:), intent(in)      :: dat
        type(T_MOD), intent(inout)                  :: model
        type(T_RUN_INFO), intent(inout)             :: RTI
        type(d3), dimension(2), intent(in)          :: perturbed_box
        type(T_LIKE_SET), intent(in)                :: like_set
        type(T_LIKE), intent(inout)                 :: like

        if(like_set%datatype == 2) &
            call exception_raiseError('Error. Surface wave tomography cannot change source locations!')

        call body_location_likelihood(dat(1),model,RTI,perturbed_box,like_set,like%likelihoods(1))

        select case (like_set%datatype)
        case(0,1)
            like%like = like%likelihoods(1)%like
            like%misfit = like%likelihoods(1)%misfit
            like%unweighted_misfit = like%likelihoods(1)%unweighted_misfit
        case(3)
            like%like = like%likelihoods(1)%like + like%likelihoods(2)%like
            like%misfit = like%likelihoods(1)%misfit + like%likelihoods(2)%misfit
            like%unweighted_misfit = like%likelihoods(1)%unweighted_misfit + like%likelihoods(2)%unweighted_misfit
        case (4)
            call exception_raiseError('location_likelihood not implemented for LIKE_SET%DATATYPE=4')
        end select

        !write(*,*) '-loglikelihood: ', like%like
        !write(*,*) 'misfits: ', like%misfit
        !write(*,*) 'unweighted misfits: ', like%unweighted_misfit
    end subroutine location_likelihood

end module m_likelihood
