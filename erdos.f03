!##################################################################################
!   ASSIGNMENT 1 - COMPUTATIONAL METHODS FOR COMPLEX SYSTEMS
!##################################################################################

program ludecomposition

  ! Use Marsenne Twist RNG
    use mtmod
    
    integer :: time, counter,ctemp
    real    :: ratio
    integer :: i,j,t0,t1
    integer, parameter  :: n = 10, nr=1000000
    real, parameter     :: pmin = 0.05, pmax = 0.9, dp = 0.05 
    
    real*8 :: a(n,n)

    call sgrnd(time())
    
    open (11, file = "n10.res", status = "unknown")
    write(11,'(2f15.8)'),0.0,1.0
    
    do j=0,nint((pmax-pmin)/dp)
    
    t0=time()
    counter = 0
    
        do i=1,nr
            ctemp = counter
            call erdos(a,n,pmin+j*dp,counter)
            if (counter > ctemp) cycle 
 
            call ludecpiv(a,n,counter)
        enddo
    
    ratio = float(counter)/nr
    
    write(11,'(2f15.8)'),pmin+j*dp,ratio
    t1=time()
    print*,pmin+j*dp,ratio,float((t1-t0))
    enddo
    
    close(11)
    
end program

!##################################################################################

! Subroutine to generate erdos matrices
subroutine erdos(a,n,p,c)

    use mtmod
    
    integer, intent(in)    :: n            ! Number of Nodes
    integer, intent(inout) :: c            ! Counter
    real, intent(in)       :: p            ! Probability of Link
    real*8, intent(inout)  :: a(n,n)       ! Adjacency Matrix
    
    
    integer :: sing(n)                     ! To check isolated nodes
    integer :: i,j
    
    ! Initialisation
    a = 0.d0
    sing = 1
    
    ! The sing vector is initialised to one. Every time node i gets a link, sing(i) is decreased by 1.
    ! If in the end any element of sing is still 1, it means that node has no link.
    
    do j=1,n-1
        do i=j+1,n
            if (grnd()<p) then
                a(i,j)=1.d0
                a(j,i)=1.d0
                sing(i) = sing(i)-1
                sing(j) = sing(j)-1
            endif
        enddo
    enddo
    
    if (any(sing==1)) c=c+1
    return
    
end subroutine

!##################################################################################

!Just a subroutine to swap two vectors
subroutine swap(a,b,n)
    
    integer, intent(in) :: n
    real*8, dimension(n), intent(inout) :: a, b
    real*8, dimension(n) :: work
    
    work = a
    a = b
    b = work
    return
    
end subroutine swap

!##################################################################################

!Subroutine for the LU decomposition of an n*n matrix a.
!c is the counter, increased if the matrix is singular
!The DOOLITTLE's algorithm is implemented
subroutine ludecpiv(a,n,c)

    integer, intent(in) :: n
    real*8, intent(inout)  :: a(n,n)
    integer, intent(inout) :: c
    
    real*8  :: eps    
    integer :: i,j,k
    
    eps = 10*epsilon(0.d0)  
    
    do i=1,n
        
        ! First we calculate the i-th diagonal element of U and
        ! the elements of the i-th column of L
        
        do j=i,n
            do k=1,i-1
                a(j,i) = a(j,i)-a(j,k)*a(k,i)
            enddo
        enddo
        
!        if (i==n) then
!        if (abs(a(i,i)) < eps) then  ! This means that we have a 0 on the diagonal of U, so the determinant is 0, so singular
!            c = c+1
!            return
!        else
!        return
!        endif
!        endif    ! There are no off diagonal elements for i=n
        
        !START OF THE PIVOTING SECTION!
        ! Then we want to make sure that the diagonal element of U is the biggest in the i-th column
        ! Find the maximum in the i-th column of U.
        k = maxval(maxloc(abs(a(i:n,i))))+i-1    
               
        if (abs(a(k,i)) < eps) then  ! This means that we have a 0 on the diagonal of U, so the determinant is 0, so singular
            c = c+1
            return
        endif
       
        if (k /= i) call swap(a(i,:),a(k,:),n)  ! Pivoting. No need to keep track because not solving system. Swapping columns
        !END OF THE PIVOTING SECTION!
            
        ! Then we calculate the other elements of the i-th row of U   
           
        do j=i+1,n
            do k=1,i-1
                a(i,j) = a(i,j)-a(i,k)*a(k,j)
            enddo
            a(j,i) = a(j,i)/a(i,i)
        enddo
        
    enddo
    
    return
    
end subroutine

!##################################################################################

!Just a subroutine to print matrices
subroutine printmatrix(a,n)
    
    integer :: n,i,j
    real*8  :: a(n,n)
    
    do i=1,n
        print '(100f5.2)',(a(i,j),j=1,n)
    enddo
end subroutine
    
