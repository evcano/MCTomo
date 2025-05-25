module m_writesample
    use iso_c_binding

    type T_SAMPLE
        integer(c_int)  step    ! propose type
        logical(c_bool) accepted  ! accepted or not
        integer(c_int)  vindex  ! index of vertex
        integer(c_size_t)  ncells  ! number of cells
        real( kind=c_double ) misfit
        real( kind=c_double ) unweighted_misfit
        real( kind=c_double ) like
        real( kind=c_double ), dimension(3) :: coord
        real( kind=c_double ), dimension(3) :: values
        real( kind=c_double )  :: noise0
        real( kind=c_double )  :: noise1
    endtype

    ! static values
    integer, parameter, public :: sample_unit = 11
    integer, parameter, public :: outfile_unit = 12

contains

    subroutine read_samples(samples,nsamples)
        implicit none

        type(T_SAMPLE), dimension(:), intent(inout) :: samples
        integer, intent(out) :: nsamples

        integer i
        integer stat

        nsamples = 0
        do i = 1, size(samples)
            read(sample_unit,iostat=stat) samples(i)%step, samples(i)%accepted,&
            samples(i)%vindex, samples(i)%ncells, samples(i)%misfit,&
            samples(i)%unweighted_misfit,samples(i)%like,  &
            samples(i)%coord, samples(i)%values, samples(i)%noise0, &
            samples(i)%noise1
            if(stat /= 0) exit
            nsamples = nsamples + 1
        enddo

    end subroutine
    
    subroutine write_samples_txt(samples,nsamples)
        implicit none

        type(T_SAMPLE), dimension(:), intent(in) :: samples
        integer, intent(in) :: nsamples
        integer i
    
        if(nsamples==0) return
    
        ! write samples to the file sample by sample
        do i = 1, nsamples
            write(outfile_unit,*) samples(i)%step, samples(i)%accepted,&
            samples(i)%vindex, samples(i)%ncells, samples(i)%misfit,&
            samples(i)%unweighted_misfit,samples(i)%like,  &
            samples(i)%coord, samples(i)%values,&
            samples(i)%noise0, samples(i)%noise1
        enddo

    end subroutine

end module m_writesample

program write2txt
    use m_writesample

    implicit none

    ! static value
    integer, parameter :: N = 1000000

    logical lexist

    character(len=100) :: filename_in
    character(len=100) :: filename_out

    type(T_SAMPLE), dimension(:), allocatable  ::  samples
    integer nsamples

    filename_in = './samples_1.out'
    filename_out = './samples_1.txt'

    allocate(samples(N))

    ! open file to read
    inquire(file=filename_in,exist=lexist)
    if(.not.lexist) write(*,*) 'samples file does not exist'
    open(unit=sample_unit,file=filename_in,status='old',access='stream',action='read')

    ! open file to write
    open(unit=outfile_unit,file=filename_out,status='replace',action='write')

    do
        call read_samples(samples,nsamples)
        call write_samples_txt(samples,nsamples)
        if(nsamples<N) exit
    enddo

    ! close files
    close(sample_unit)
    close(outfile_unit)

end program
