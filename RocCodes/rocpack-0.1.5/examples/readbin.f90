        program readbin
        implicit none
        integer :: i, j, k, ni, nj, nk
        byte, allocatable :: packing(:,:,:)

        open(1, file='pack.bin', form='unformatted', access='stream')
        
        read(1) ni, nj, nk

        print *, ni, nj, nk

        allocate (packing(ni, nj, nk))

        read(1) packing(:,:,:)
        
        ! packing(i,j,k) now contains the tag number of voxel(i,j,k)

        ! print a slice through the middle of the pack on the screen

        k = nk/2
        do j = 1, nj
                do i = 1, ni
                        write (*, '(I1)', advance='no') packing(i,j,k)
                enddo
                print *, ''
        enddo

        print *, ''

        close(1)

        end program readbin
