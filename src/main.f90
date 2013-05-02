module typeDefinitions
    implicit none
    public
    type geometry ! object that may have a length, a volume or whatever related to a geometry
        real :: length
    end type
    type type_species
        character(len=2) :: name
        integer :: nbInCell
        real :: charge
        real :: compo
        real :: coord ! coordination
    end type
    type type_liq
        type(type_species), dimension(:), allocatable :: species
        integer :: nTot, nSpecies
    end type
end module

program pimaimSupercellGenerator
    use typeDefinitions
    implicit none
    integer, parameter :: Zr=1, Na=2, F=3, Y=4
    integer, allocatable, dimension(:,:,:) :: inside
    integer :: nBinPerDir, countPrintedZr, countPrintedNa, countPrintedF
    integer :: i, j, k
    character(len=*), parameter :: filenamedat='restart.inpt'
    character(len=*), parameter :: filenamexyz='restart.inpt.xyz'
    character(len=*), parameter :: inputFileName='in'
    type (type_liq) :: liq
    type (geometry) :: cell
    type myneighbors
        real, dimension(3) :: vecOne = 2.*[1,0,0]
        real, dimension(3) :: vecTwo1 = 2.*[1,0,0],&
                              vecTwo2 = 2.*[-1,0,0]
        real, dimension(3) :: vecThree1 = 2.*[1,0,0],&
                              vecThree2 = 2.*[-1,0,0],&
                              vecThree3 = 2.*[0,1,0]
        real, dimension(3) :: vecTetra1 = 4./sqrt(3.)*0.5*[-1,1,-1] ,&
                                         vecTetra2 = 4./sqrt(3.)*0.5*[-1,-1,1] ,&
                                         vecTetra3 = 4./sqrt(3.)*0.5*[1,1,1]   ,&
                                         vecTetra4 = 4./sqrt(3.)*0.5*[1,-1,-1]
    end type
    type(myneighbors) :: neighbor
    type myRandomNumbers
        real, dimension(3) :: dim3Real
        integer, dimension(3) :: dim3Int
    end type
    type(myRandomNumbers) :: harvest
    real :: myx, myy, myz, distanceBetweenTwoCations
    if( .not. isPresent(inputFileName) ) then ! check presence of input file
        stop 'STOP. The input file called in is not present'
    else
        open(21,file=inputFileName)
    end if
    read(21,*) liq%nSpecies ! number of species in the liquid, e.g. 3 in NaFZrF4
    if( liq%nSpecies<1 ) stop 'STOP. too few species'
    allocate( liq%Species(liq%nSpecies) )
    if(      liq%nSpecies==1 ) then ! read names then charge horizontaly. Not elegant but don't know how to do better.
        read(21,*) liq%species(1)%name
        read(21,*) liq%species(1)%charge
        liq%species(1)%compo = 1.0d0
    else if( liq%nSpecies==2 ) then
        read(21,*) liq%species(1)%name, liq%species(2)%name
        read(21,*) liq%species(1)%charge, liq%species(2)%charge
        read(21,*) liq%species(1)%compo, liq%species(2)%compo
    else if( liq%nSpecies==3 ) then
        read(21,*) liq%Species(1)%name, liq%species(2)%name, liq%species(3)%name
        read(21,*) liq%Species(1)%charge, liq%species(2)%charge, liq%species(3)%charge
        read(21,*) liq%Species(1)%compo, liq%species(2)%compo, liq%species(3)%compo
    else if( liq%nSpecies==4 ) then
        read(21,*) liq%Species(1)%name, liq%species(2)%name, liq%species(3)%name, liq%species(4)%name
        read(21,*) liq%Species(1)%charge, liq%species(2)%charge, liq%species(3)%charge, liq%species(4)%charge
        read(21,*) liq%Species(1)%compo, liq%species(2)%compo, liq%species(3)%compo, liq%species(4)%compo
    else if( liq%nSpecies>=5 ) then
        stop 'STOP. 5 or more species are not allowed for now. Implement it!'
    end if
    do i=1,liq%nSpecies ! all species names should be different
        do j=i,liq%nSpecies
            if (liq%species(i)%name == liq%species(j)%name .and. i /= j) stop 'STOP. Two species have same name.'
        end do
    end do
    if( count(liq%species%charge<0.)==0 .or. count(liq%species%charge>0.)==0 ) &
            stop 'STOP. There should be positive **and** negative charges in the supercell'
    liq%species%coord = abs(liq%species%charge)
    if( sum(liq%species%compo)<0.9999 .or. sum(liq%species%compo)>1.0001 ) stop 'STOP. Sum of compositions should be 1'
    if ( any(liq%species%compo<0.) ) stop 'STOP. No molar fraction should be negative.'
    if ( any(liq%species%compo>1.d0) )
        .or. any(liq%species%compo>1.)
    read(21,*) cell%length ! length of the orthorombic supercell
    if( cell%length < 1. ) stop 'STOP. length is too small'
    read(21,*)liq%nTot ! approximate total number of atoms in the supercell
    if( liq%nTot < 1 ) stop 'STOP. too few atoms. Check input.'
    
    ! now we compute the number of atoms of each species
    myz = real(liq%nTot)
    myx = real(liq%species(1)%coord+1)
    myy = real(liq%species(1)%coord-liq%species(2)%coord)
    liq%species(2)%nbInCell = liq%species(2)%compo*myz/myx/(1.-liq%species(2)%compo*myy/myx)
    
    
    
    liq%species(2)%nbInCell = int( real(liq%nTot)/(5.+2.*(1.-liq%species(2)%compo)/liq%species(2)%compo) )
    print*,'There will be ',liq%species(2)%nbInCell,' Zr in the supercell'
    liq%species(3)%nbInCell = int( (1.-liq%species(2)%compo)/liq%species(2)%compo * liq%species(2)%nbInCell )
    print*,'There will be ',liq%species(3)%nbInCell,' Na in the supercell'
    liq%species(1)%nbInCell = liq%species(3)%nbInCell + 4*liq%species(2)%nbInCell
    print*,'There will be ',liq%species(1)%nbInCell,' F in the supercell'
    liq%nTot = liq%species(1)%nbInCell+liq%species(3)%nbInCell+liq%species(2)%nbInCell
    print*,'There will be ',liq%nTot,' atoms in total (exactly) in the supercell.'
    ! the total number is not exactly the one wanted because of xZrF4 being a real, not an integer. And x*int is not an int. So we approximate it to a real.
    print*,'volume per atome = ',cell%length**3/liq%nTot
    distanceBetweenTwoCations = (cell%length**3/liq%nTot)**(1./3.)
    print*,'distance between two cations = ',distanceBetweenTwoCations
    if( distanceBetweenTwoCations < 4. ) stop 'the supercell length is too small'
    nBinPerDir = int(cell%length/ ((cell%length**3/(liq%species(2)%nbInCell+liq%species(3)%nbInCell))**(1./3.)))+1
    print*,'number of bins in each direction = ', nBinPerDir
    if (nBinPerDir**3 < (liq%species(2)%nbInCell+liq%species(3)%nbInCell) ) stop 'STOP. There should be more bins. pb in algo :('
    if( nBinPerDir < 2 ) stop 'STOP. too many atom per cubic cell.'


    allocate( inside(nBinPerDir,nBinPerDir,nBinPerDir), source=0)
    do while( count(inside==Na) < (liq%species(3)%nbInCell+liq%species(2)%nbInCell) ) ! pick up liq%species(2)%nbInCell+liq%species(3)%nbInCell sites randomly and attribute it to Na
        call random_number(harvest%dim3Real)
        harvest%dim3Int = nint( harvest%dim3Real*(nBinPerDir-1) +1 )
        if( maxval(harvest%dim3Int) > nBinPerDir ) stop 'STOP. random number too big ! :('
        if( minval(harvest%dim3Int) < 1 ) stop 'STOP. random number too small ! :('
        inside(harvest%dim3Int(1),harvest%dim3Int(2),harvest%dim3Int(3)) = Na
    end do
    if( count(inside==Na)/=liq%species(3)%nbInCell+liq%species(2)%nbInCell ) stop 'STOP. not the good number of Na in inside'
    do while( count(inside==Zr) < liq%species(2)%nbInCell ) ! replace default Na by Zr as long as there is not enough Zr
        call random_number(harvest%dim3Real)
        harvest%dim3Int = nint( harvest%dim3Real*(nBinPerDir-1) +1 )
        if( maxval(harvest%dim3Int) > nBinPerDir ) stop 'STOP. random number too big ! :('
        if( minval(harvest%dim3Int) < 1 ) stop 'STOP. random number too small ! :('
        if( inside(harvest%dim3Int(1),harvest%dim3Int(2),harvest%dim3Int(3)) == Na ) then
            inside(harvest%dim3Int(1),harvest%dim3Int(2),harvest%dim3Int(3)) = Zr
        end if
    end do
    if( count(inside==Zr)/=liq%species(2)%nbInCell ) stop 'STOP. not the good number of Zr' ! all Zr and Na should now be in the cell
    if( count(inside==Na)/=liq%species(3)%nbInCell ) stop 'STOP. not the good number of Na'
    ! print the positions
    open(10, file=filenamexyz)
    open(11, file=filenamedat)
    write(10,*)liq%nTot
    write(10,*)
    write(11,*)'T'
    write(11,*)'F'
    write(11,*)'F'
    write(11,*)'F'
    ! print all F
    countPrintedF = 0
    do i=1,nBinPerDir
        do j=1,nBinPerDir
            do k=1,nBinPerDir
                myx = real(i-1)*cell%length/real(nBinPerDir)
                myy = real(j-1)*cell%length/real(nBinPerDir)
                myz = real(k-1)*cell%length/real(nBinPerDir)
                if( inside(i,j,k)==Zr ) then
                    write(10,*) 'F',modulo( [myx,myy,myz] + neighbor%vecTetra1 ,cell%length )
                    write(10,*) 'F',modulo( [myx,myy,myz] + neighbor%vecTetra2 ,cell%length )
                    write(10,*) 'F',modulo( [myx,myy,myz] + neighbor%vecTetra3 ,cell%length )
                    write(10,*) 'F',modulo( [myx,myy,myz] + neighbor%vecTetra4 ,cell%length )
                    write(11,*) modulo( [myx,myy,myz] + neighbor%vecTetra1 ,cell%length )
                    write(11,*) modulo( [myx,myy,myz] + neighbor%vecTetra2 ,cell%length )
                    write(11,*) modulo( [myx,myy,myz] + neighbor%vecTetra3 ,cell%length )
                    write(11,*) modulo( [myx,myy,myz] + neighbor%vecTetra4 ,cell%length )
                    countPrintedF = countPrintedF +4
                else if( inside(i,j,k)==Na ) then
                    write(10,*) 'F',modulo( [myx,myy,myz] + neighbor%vecOne ,cell%length )
                    write(11,*) modulo( [myx,myy,myz] + 2.*[1,0,0] ,cell%length )
                    countPrintedF = countPrintedF +1
                end if
            end do
        end do
    end do
    if( countPrintedF /= liq%species(1)%nbInCell ) then
        print*,
        print*,'STOP. bad number of F :('
        print*,'countPrintedF = ',countPrintedF,', instead of ',liq%species(1)%nbInCell
        stop
    end if
    ! print all Zr
    countPrintedZr = 0
    do i=1,nBinPerDir
        do j=1,nBinPerDir
            do k=1,nBinPerDir
                myx = real(i-1)*cell%length/real(nBinPerDir)
                myy = real(j-1)*cell%length/real(nBinPerDir)
                myz = real(k-1)*cell%length/real(nBinPerDir)
                if( inside(i,j,k)==Zr ) then
                    write(10,*) 'Zr',modulo( [myx,myy,myz] ,cell%length )
                    write(11,*) modulo( [myx,myy,myz] ,cell%length )
                    countPrintedZr = countPrintedZr +1
                end if
            end do
        end do
    end do
    if( countPrintedZr /= liq%species(2)%nbInCell ) then
        print*,
        print*,'STOP. bad number of Zr :('
        print*,'countPrintedF = ',countPrintedZr,', instead of ',liq%species(2)%nbInCell
        stop
    end if
    ! print all Na
    countPrintedNa = 0
    do i=1,nBinPerDir
        do j=1,nBinPerDir
            do k=1,nBinPerDir
                myx = real(i-1)*cell%length/real(nBinPerDir)
                myy = real(j-1)*cell%length/real(nBinPerDir)
                myz = real(k-1)*cell%length/real(nBinPerDir)
                if( inside(i,j,k)==Na ) then
                    write(10,*) 'Na',modulo( [myx,myy,myz] ,cell%length )
                    write(11,*) modulo( [myx,myy,myz] ,cell%length )
                    countPrintedNa = countPrintedNa +1
                end if
            end do
        end do
    end do
    if( countPrintedNa /= liq%species(3)%nbInCell ) then ! Check that "Na" was printed the attended number of times.
        print*,
        print*,'STOP. bad number of Na :('
        print*,'countPrintedNa = ',countPrintedNa,', instead of ',liq%species(3)%nbInCell
        stop
    end if
    close(10)
    write(11,*) '1.0 0.0 0.0'
    write(11,*) '0.0 1.0 0.0'
    write(11,*) '0.0 0.0 1.0'
    write(11,*) cell%length
    write(11,*) cell%length
    write(11,*) cell%length
    close(11)
    if( .not. printRuntimeInpt() ) stop 'Problem writing runtime.inpt'
    print*,
    print*,
    print*,'OK. Everything is fine. Look at files ',filenamedat,' and ', filenamexyz
    contains
        logical function isPresent(filename) ! Check that file called filename is present. returns true if present, false if not found
            character(*), intent(in) :: filename
            inquire(file=filename, exist=isPresent)
        end function
        logical function printRuntimeInpt() ! Print runtime.inpt needed by pimaim for a first short unpolar (RIM) annealing of the box to randomize a little.
            character(*), parameter :: runtimeInptFilename = 'runtime.inpt'
            open(99, file=runtimeInptFilename)
            write(99,*) '10000      Number of steps in the run'
            write(99,*) '1399       Translational temperature'
            write(99,*) '1          Number of ionic molecular species'
            write(99,*) '3          Number of ionic species'
            write(99,"(I5.3,A,I5.3,A,I5.3,A)") liq%species(1)%nbInCell,', ',liq%species(2)%nbInCell,', ',liq%species(3)%nbInCell,&
                        '  Number of species of type F, Zr, Na'
            write(99,*) '-1.0,4.0,1.0                  Permanent charge on ions.'
            write(99,*) '18.9984,91.224,22.9898          Atomic masses'
            write(99,*) '.true.,.true.,.true.                Polarizable ?'
            write(99,*) '.false.,.false.,.false.               Deformable ?'
            write(99,*) '20.671     Timestep (a.u.).'
            write(99,*) 'rim   Type of run (rim,dippim,quadpim,cim,aim)'
            write(99,*) 'epp'
            write(99,*) '.false.     AIM effects on anion-anion interactions?'
            write(99,*) '.false.     Like-like multipole damping?'
            write(99,*) '.false.     Isolated cluster?'
            write(99,*) '.false.     Environmental effects on PIM?'
            write(99,*) '.false.     Environmental effects on AIM?'
            write(99,*) '.true.      Conjugate gradient minimisation? (PIM) 1.0d-08 1.0d-08'
            write(99,*) '1.0d-08'
            write(99,*) '1.0d-08'
            write(99,*) '.false.     Conjugate gradient minimisation? (AIM)'
            write(99,*) '.true.       Restart?'
            write(99,*) '.true.  Set up velocities?'
            write(99,*) '.false.     Velocity rescale prior to main run ?'
            write(99,*) '.false.     Random displacement of ions.'
            write(99,*) '.true.     Move ions?'
            write(99,*) '.false.     Do a dynamical matrix calculation?'
            write(99,*) '.false.      Relax input structure?'
            write(99,*) '100          Number of steps inbetween periodic output (energies).'
            write(99,*) '100          Number of steps inbetween periodic output (velocities etc).'
            write(99,*) '100          Number of steps inbetween periodic output (frictions etc).'
            write(99,*) '100          Number of steps inbetween periodic output (pressure etc).'
            write(99,*) '500           Number of steps inbetween rdf call in main loop.'
            write(99,*) '1           Number of ions to monitor.'
            write(99,*) '1           Ion number to monitor. (1)'
            write(99,*) '5.60d0      eta = <x>/boxlen.'
            write(99,*) '16.1d0     rcut (au).'
            write(99,*) '1.0d-7      convergence parameter'
            write(99,*) '0.1d0       convergence factor.'
            write(99,*) '16.1d0     rcut (au) short range.'
            write(99,*) '.true.     Nose-Hoover thermostat? (if true then enter a relaxation time)'
            write(99,*) '10000'
            write(99,*) '.false.     Periodic rescale of temperature?  1000'
            write(99,*) '.true.     Isotropic barostat? 10000 1.0d-8'
            write(99,*) '10000'
            write(99,*) '1.0d-8'
            write(99,*) '.false.     Anisotropic barostat?'
            write(99,*) '.false.     Orthorhombic cell?'
            write(99,*) '.false.     VBC ?'
            write(99,*) '.false.      verbose'
            write(99,*) '.false.      verbose'
            write(99,*) '.false.      verbose'
            write(99,*) '.false.      verbose'
            write(99,*) '.false.      verbose'
            write(99,*) '.false.      verbose'
            write(99,*) '.false.      verbose'
            close(99)
            printRuntimeInpt = .true.
        end function
end program
