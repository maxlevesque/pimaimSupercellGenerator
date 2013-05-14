module derived_types
    use precision_kinds
    implicit none
    private
    type atom
        real(dp), dimension(3) :: pos
        character(2) :: name
    end type
    type molecule
        real(dp) :: moFrac ! mole Fraction
        real(dp), dimension(:), allocatable :: pos
        type (atom), dimension(:), allocatable :: at
    end type
    type(molecule), dimension(:), allocatable, public :: mo

    type supercell
        real(dp) :: len ! length of the supercell. For now cubic cell, i.e. constant length along each direction
        integer(i2b), dimension(:,:,:), allocatable :: nod
    end type
    type (supercell), public :: cell
    real(dp), public, dimension(:,:,:), allocatable :: coordinates

end module

module input
    use precision_kinds
    implicit none
    private
    public :: init
    contains
        subroutine init
            call checkInputFileExists
            call readNumberOfDifferentMoleculesInMixture
            call readMoleFractionOfEachDifferentMolecule
            call readDescriptionsOfMolecules
            call readCellLength
            call readNumberOfMoleculesInEachDirectionOfSupercell
            call closeInput
        end subroutine
        
        subroutine checkInputFileExists
            logical :: isOK
            inquire (file="in", exist=isOK)
            if (.not. isOK) stop "input file does not exists"
        end subroutine
        
        subroutine readNumberOfDifferentMoleculesInMixture
            use derived_types, only: mo
            integer(i2b) :: nMo
            open (11, file="in")
            read(11,*) nMo
            allocate (mo(1:nMo))
            print*,'There are ',size(mo),' molecules.'
        end subroutine
        
        subroutine readMoleFractionOfEachDifferentMolecule
            use derived_types, only: mo
            integer(i2b) :: i
            do i= 1, size(mo)
                read(11,*) mo(i)%moFrac
            end do
            if (any(mo%moFrac<0._dp)) stop "STOP. No mole fraction should be negative"
            if (any(mo%moFrac>1._dp)) stop "STOP. No mole fraction should be more than 1"
        end subroutine
        
        subroutine readDescriptionsOfMolecules
            use derived_types, only: mo
            integer(i2b) :: m, i
            do m= 1, size(mo) ! size(mo) is the total number of molecules, the first line of the input file
                read(11,*) i ! number of atoms in molecule m
                allocate (mo(m)%at(i))
                do i= 1, size(mo(m)%at)
                    read(11,'(A)') mo(m)%at(i)%name
                end do
                print*,'molecule ',m,' has ',size(mo(m)%at),' atoms: ',mo(m)%at(:)%name,'. Mole fraction ',mo(m)%moFrac
            end do
        end subroutine
        
        subroutine readCellLength
            use derived_types, only: cell
            read(11,*) cell%len
        end subroutine
        
        subroutine readNumberOfMoleculesInEachDirectionOfSupercell
            use derived_types, only: cell
            integer(i2b) :: i
            read(11,*) i
            allocate (cell%nod(i,i,i)) ! the supercell can only be cubic for now
        end subroutine
        
        subroutine closeInput
            close(11)
        end subroutine
        
end module

module supercell
    use precision_kinds
    use derived_types, only: cell
    implicit none
    private
    public :: init
    
    contains
    
        subroutine init
            call toEachNodeAllocateAMolecule
            call transformEachMoleculeIntoAtomicPositions
        end subroutine
        
        subroutine toEachNodeAllocateAMolecule
            use derived_types, only: mo
            integer(i2b) :: i
            real(dp), dimension(3) :: rn
            integer(i2b), dimension(3) :: irn
            integer(i2b) :: tar, acc ! target, accepted
            cell%nod = 0 ! initialization to something one will never want
            ! The number of molecules of type 1 is mo(1)%moFrac * 5*5*5 if cell%nod allocated as (5:5:5)
            ! note that all nodes may not be fulfilled with a molecule because of the real number given in the line above to integer transforms.
            do i= 1, size(mo) ! for each molecule
                tar = int (mo(i)%moFrac * size(cell%nod))
                if (tar == 0) tar = 1
                acc = 0
                do while (acc < tar)
                    call random_number (rn(1:3))
                    irn = int (rn*real(size(cell%nod,1))) + 1
                    if (cell%nod(irn(1),irn(2),irn(3)) == 0) then
                        cell%nod(irn(1),irn(2),irn(3)) = i
                        acc = acc +1
                    end if
                end do
            end do
        end subroutine
        
        subroutine transformEachMoleculeIntoAtomicPositions
            use derived_types, only: mo
            integer(i2b) :: m, nAtPerMo
            do m = 1, size(mo)
                nAtPerMo = size(mo(m)%at)
                if( nAtPerMo == 0 ) then
                    stop "STOP. Impossible. This problem should have been detected sooner"
                else if( nAtPerMo == 1 ) then
                    mo(m)%at(1)%pos(:) = real([0, 0, 0],dp)
                else if( nAtPerMo == 2 ) then
                    mo(m)%at(1)%pos(:) = 2._dp/sqrt(3._dp)*real([-1, -1, -1],dp)
                    mo(m)%at(2)%pos(:) = real([0, 0, 0],dp)
                else if( nAtPerMo == 3 ) then
                    mo(m)%at(1)%pos(:) = 2._dp/sqrt(3._dp)*real([1, 1, 1],dp)
                    mo(m)%at(2)%pos(:) = 2._dp/sqrt(3._dp)*real([-1, -1, -1],dp)
                    mo(m)%at(3)%pos(:) = real([0, 0, 0],dp)
                else if( nAtPerMo == 4 ) then
                    mo(m)%at(1)%pos(:) = 2._dp*real([1, 0, 0],dp)
                    mo(m)%at(2)%pos(:) = real([-1, 1, 1],dp)
                    mo(m)%at(3)%pos(:) = real([-1, -1, -1],dp)
                    mo(m)%at(4)%pos(:) = real([0, 0, 0],dp)
                else if( nAtPerMo == 5 ) then
                    mo(m)%at(1)%pos(:) = 2._dp/sqrt(3._dp)*real([1,-1,-1],dp)
                    mo(m)%at(2)%pos(:) = 2._dp/sqrt(3._dp)*real([-1,1,-1],dp)
                    mo(m)%at(3)%pos(:) = 2._dp/sqrt(3._dp)*real([-1,-1,1],dp)
                    mo(m)%at(4)%pos(:) = 2._dp/sqrt(3._dp)*real([1,1,1],dp)
                    mo(m)%at(5)%pos(:) = real([0, 0, 0],dp)
                else
                    stop "STOP. More than 5 sites per molecule is not yet supported. Implement it please."
                end if
            end do
        end subroutine

end module

module final
    use precision_kinds
    use derived_types, only: mo, cell
    implicit none
    private
    public init
    type atoms
        real(dp), dimension(3) :: pos ! coordinate
        character(2) :: name
    end type
    type (atoms), dimension(:), allocatable :: at
    integer(i2b) :: nbOfAtomicTypes
    character(2), dimension(:), allocatable :: atomicTypes
    
    contains
    
        subroutine init
            call listAllPositionsInCellOrder
            call determineNumberOfAtomTypes
            call writeXSF
            call writeRestartDotDat
            call printInfoOut
        end subroutine
        
        subroutine printInfoOut
            integer(i2b) :: i
            do i= 1, nbOfAtomicTypes
                print*,'nb of type ',i,count (at%name == atomicTypes(i))
            end do
        end subroutine
        
        subroutine listAllPositionsInCellOrder
            use derived_types, only: cell, mo
            integer :: i, j, k, l, totalNbOfAtomsInSupercell, moInThisNod, a
            real(dp), dimension(3) :: nod
            totalNbOfAtomsInSupercell = 0
            do concurrent (i=1:size(cell%nod,1) , j=1:size(cell%nod,2) , k=1:size(cell%nod,3), cell%nod(i,j,k)/=0)
                totalNbOfAtomsInSupercell = totalNbOfAtomsInSupercell + size(mo(cell%nod(i,j,k))%at)
            end do
            allocate (at(totalNbOfAtomsInSupercell))
            l = 0
            do i= 1, size(cell%nod,1)
                do j= 1, size(cell%nod,2)
                    do k= 1, size(cell%nod,3)
                        call coordinatesOfCenterOfNode (nod,i,j,k)
                        moInThisNod = cell%nod(i,j,k)
                        if( moInThisNod == 0) cycle ! for the few nodes not containing molecule
                        do a = 1, size(mo(moInThisNod)%at)
                            l = l + 1
                            at(l)%pos(:) = mo(moInThisNod)%at(a)%pos(:) + nod(:)
                            at(l)%name = mo(moInThisNod)%at(a)%name
                        end do
                    end do
                end do
            end do
        end subroutine

        subroutine coordinatesOfCenterOfNode (nod,i,j,k)
            use derived_types, only: cell
            integer(i2b), intent(in) :: i, j, k
            real(dp), dimension(3), intent(out) :: nod
            integer(i2b) :: imax, jmax, kmax
            imax = size(cell%nod,1)
            jmax = size(cell%nod,2)
            kmax = size(cell%nod,3)
            nod = cell%len * [real(i-1,dp)/real(imax,dp), real(j-1,dp)/real(jmax,dp), real(k-1,dp)/real(kmax,dp)]
        end subroutine

        subroutine determineNumberOfAtomTypes
            integer(i2b) :: i, j
            j = 0
            do i = 1, size(mo)
                j = j + size(mo(i)%at)
            end do
            allocate( atomicTypes (j) )
            atomicTypes = "0" ! zero is not an atomic symbol
            nbOfAtomicTypes = 0
            do i = 1, size(mo)
                do j = 1, size(mo(i)%at)
                    if (any(atomicTypes==mo(i)%at(j)%name)) then ! it's a known name
                        cycle
                    else
                        nbOfAtomicTypes = nbOfAtomicTypes + 1
                        atomicTypes (nbOfAtomicTypes) = mo(i)%at(j)%name
                        print*,'type ',nbOfAtomicTypes,' is ',atomicTypes (nbOfAtomicTypes)
                    end if
                end do
            end do
            block ! reshape atomicTypes (an intrinsic function exists that does that but the semantics are cryptic)
                character(2), dimension(size(atomicTypes)) :: tmp
                tmp = atomicTypes
                deallocate (atomicTypes)
                allocate (atomicTypes(nbOfAtomicTypes))
                atomicTypes(1:nbOfAtomicTypes) = tmp(1:nbOfAtomicTypes)
            end block
        end subroutine
       
        subroutine writeXSF
            integer(i2b) :: i, j
            open(10, file='restart.dat.xyz')
            write(10,*) size(at) ! first line should be the number of positions to be written in cartesian units
            write(10,*) ! second line should a comment line (blank here)
            do i= 1, nbOfAtomicTypes ! then all positions
                do j= 1, size(at)
                    if (at(j)%name == atomicTypes(i)) then
                        write(10,*) at(j)%name, modulo( at(j)%pos ,cell%len )
                    else
                        cycle
                    end if
                end do
            end do
            close(10)
        end subroutine
        
        subroutine writeRestartDotDat
            integer(i2b) j
            open(10, file='restart.dat')
            write(10,*)'T'
            write(10,*)'F'
            write(10,*)'F'
            write(10,*)'F'
            do j= 1, size(at)
                write(10,*) modulo( at(j)%pos ,cell%len )
            end do
            write(10,*) '1.0 0.0 0.0'
            write(10,*) '0.0 1.0 0.0'
            write(10,*) '0.0 0.0 1.0'
            write(10,*) cell%len
            write(10,*) cell%len
            write(10,*) cell%len
            close(10)
        end subroutine

end module

program pimaimSupercellGenerator

    use precision_kinds
    use input, only: initInput => init
    use supercell, only: initSupercell => init
    use final, only: initFinal => init
    implicit none

    call initInput
    call initSupercell
    call initFinal

end program
