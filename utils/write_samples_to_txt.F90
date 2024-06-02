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

contains

    subroutine read_samples(samples,filename)
        character( len=* ), intent(in) :: filename
        type(T_SAMPLE), dimension(:), intent(inout) :: samples
    
        integer i
        integer stat
        logical lexist
    
        inquire(file=filename,exist=lexist)
        if(.not.lexist) write(*,*) 'error while reading file'
        open(unit=11,file=filename,status='old',access='stream',action='read')
            do i = 1, size(samples)
                read(11,iostat=stat) samples(i)%step, samples(i)%accepted,&
                samples(i)%vindex, samples(i)%ncells, samples(i)%misfit,&
                samples(i)%unweighted_misfit,samples(i)%like,  &
                samples(i)%coord, samples(i)%values, samples(i)%noise0, &
                samples(i)%noise1
                if(stat /= 0)then
                    write(*,*) 'error while reading file'
                endif
            enddo
        close(11)
    end subroutine
    
    subroutine write_samples_txt(samples,filename)
        implicit none
        character( len=* ), intent(in) :: filename
        type(T_SAMPLE), dimension(:), intent(in) :: samples
    
        integer i
        integer nsamples
        logical lexist
    
        if(size(samples)==0) return
    
        ! > open file for writing
        open(unit=12,file=filename,status='replace',action='write')
    
        ! > write samples to the file sample by sample
        nsamples = size(samples)
        do i = 1, nsamples
            write(12,*) samples(i)%step, samples(i)%accepted,&
            samples(i)%vindex, samples(i)%ncells, samples(i)%misfit,&
            samples(i)%unweighted_misfit,samples(i)%like,  &
            samples(i)%coord, samples(i)%values,&
            samples(i)%noise0, samples(i)%noise1
        enddo
        close(12)
    end subroutine

end module m_writesample

program write2txt
    use m_writesample

    implicit none

    character(len=100) :: filename_in
    character(len=100) :: filename_out

    type(T_SAMPLE), dimension(:), allocatable  ::  samples
    integer nsamples

    filename_in = './samples_1.out'
    filename_out = './sample_1.txt'
    nsamples = 100000

    allocate(samples(nsamples))
    call read_samples(samples,filename_in)
    call write_samples_txt(samples,filename_out)
end program
