module m_motore
implicit none
contains

!////////////////////////////////
!//// SUBROUTINE DE CALCUL //////
!////////////////////////////////

!construction de v_cst en fonction de v_relief
subroutine c_relief(v_relief,v_cst,c,m1,m2,motif,epaiss)
	real,dimension(:),intent(in)::v_relief
	real,dimension(:),intent(out)::v_cst
	real,intent(in)::c,m1,m2,motif,epaiss
	integer::i,n
	n=size(v_relief)
	do i=1,n
	if(mod(abs(v_relief(i)),motif)<epaiss) then
		v_cst(i)=c*m1
	else
		v_cst(i)=c*m2
	end if
	end do
end subroutine c_relief

!construction de la matrice
subroutine m_matrice(mat,di,dj,saut,dt,dxy,v_cst,nx,sig,dcxy)
	real,dimension(:,:),intent(inout)::mat
	real,dimension(:),intent(in)::v_cst,dcxy
	real,intent(in)::dt,sig,dxy
	integer,intent(in)::di,dj,saut,nx
	integer::c,k,nxy
	nxy=size(mat,1)
	c=0
	do k=1,nxy-(di+dj)
		c=c+1
		if (c/=saut) then
			mat(k+di,k+dj)=-dt/(dxy**2.)*(v_cst(k+dj)+sig*dcxy(k+dj)/4.)
		else
			c=0	
		end if
	end do
end subroutine m_matrice

!Calcul de la derivee selon x v_dcx de v_cst
subroutine constru_dcx(v_dcx,v_cst,nx,ny,nxy)
	real,dimension(:),intent(inout)::v_cst,v_dcx
	integer,intent(in)::nx,ny,nxy
	integer::i,j,k
	do i=2,nx-1 ; do j=1,ny
		k=nx*(i-1)+j
		v_dcx(k)=v_cst(k-nx)-v_cst(k+nx)
	end do ; end do
	v_dcx(1:ny)=0.
	v_dcx(nxy-ny+1:nxy)=0.
end subroutine constru_dcx

!Calcul de la derivee selon y v_dcy de v_cst
subroutine constru_dcy(v_dcy,v_cst,nx,ny,nxy)
	real,dimension(:),intent(inout)::v_cst,v_dcy
	integer,intent(in)::nx,ny,nxy
	integer::i,j,k
	do i=1,nx ; do j=2,ny-1
		k=nx*(i-1)+j	
		v_dcy(k)=v_cst(k-1)-v_cst(k+1)
	end do ; end do
	v_dcy(1:nxy-ny+1:ny)=0.!c-v_cst(2:nxy-ny+2:ny)
	v_dcy(ny:nxy:ny)=0.!v_cst(ny-1:nxy-1:ny)-c
end subroutine constru_dcy

!construction de la topo
subroutine m_topo(ao,nx,ny,ampx,ampy,xp,yp,ptp,sig,div)
	real,dimension(:),intent(inout)::ao
	integer,intent(in)::nx,ny
	real,intent(in)::ampx,ampy,ptp,sig,xp,yp,div
	integer::i,j,k
	real::x,y,xxp,yyp
	do i=1,ny
	do j=1,nx
	k=nx*(i-1)+j ; x=1./div*i ; y=1./div*j ; xxp=xp/div ; yyp=yp/div
	ao(k)=ao(k)+sig*exp(-ampx*(x-xxp)**2)*exp(-ampy*(y-yyp)**2)*ptp
	end do
	end do
end subroutine m_topo

!construction de la topo utilisant la subroutine m_topo
subroutine m_topomulti(v_relief,nx,ny,div)
	real,dimension(:),intent(inout)::v_relief
	integer,intent(in)::nx,ny
	real,intent(in)::div

	!integer::i,j,k
	!real::x,y

	!do i=1,ny
	!do j=1,nx
	!k=nx*(i-1)+j ; x=1.*i ; y=1.*j 
	!v_relief(k)=v_relief(k)+(x+y)/25.
	!end do
	!end do

	call m_topo(v_relief,nx,ny,0.,0.3,0.,2.,4.3,1.,div) !chaine

	call m_topo(v_relief,nx,ny,0.08,0.,1.,0.,7.,1.,div) !chaine2

	call m_topo(v_relief,nx,ny,0.3,0.08,20.,9.,6.8,1.,div) !montOK
	call m_topo(v_relief,nx,ny,0.2,0.2,25.,25.,5.,1.,div) !montOK
	call m_topo(v_relief,nx,ny,0.18,0.18,2.,3.,10.,-1.,div) !montOK

	!call m_topo(v_relief,nx,ny,0.,0.5,0.,6.,6.9,1.,div) !chaine
	!call m_topo(v_relief,nx,ny,0.4,0.,5.,0.,4.7,1.,div) !chaine2
	!call m_topo(v_relief,nx,ny,0.4,0.2,14.,12.,7.,1.,div) !montOK
	!call m_topo(v_relief,nx,ny,0.2,0.2,22.,18.,5.,1.,div) !montOK
	!call m_topo(v_relief,nx,ny,0.5,0.5,5.,6.,4.,-1.,div) !montOK
end subroutine m_topomulti

!////////////////////////////////
!//// SUBROUTINE INVERSION //////
!////////////////////////////////
SUBROUTINE FINDInv(matrix, inverse, n, errorflag)
        IMPLICIT NONE
        !Declarations
        INTEGER, INTENT(IN) :: n
        INTEGER, INTENT(OUT) :: errorflag  !Return error status. -1 for error, 0 for normal
        REAL, INTENT(IN), DIMENSION(n,n) :: matrix  !Input matrix
        REAL, INTENT(OUT), DIMENSION(n,n) :: inverse !Inverted matrix
        
        LOGICAL :: FLAG = .TRUE.
        INTEGER :: i, j, k!, l
        REAL :: m
        REAL, DIMENSION(n,2*n) :: augmatrix !augmented matrix
        
        !Augment input matrix with an identity matrix
        !write(*,*)'1 sur 6'
        DO i = 1, n
                DO j = 1, 2*n
                        IF (j <= n ) THEN
                                augmatrix(i,j) = matrix(i,j)
                        ELSE IF ((i+n) == j) THEN
                                augmatrix(i,j) = 1
                        Else
                                augmatrix(i,j) = 0
                        ENDIF
                END DO
        END DO
        
        !Reduce augmented matrix to upper traingular form
        !write(*,*)'2 sur 6'
        DO k =1, n-1
                IF (augmatrix(k,k) == 0) THEN
                        FLAG = .FALSE.
                        DO i = k+1, n
                                IF (augmatrix(i,k) /= 0) THEN
                                        DO j = 1,2*n
                                                augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
                                        END DO
                                        FLAG = .TRUE.
                                        EXIT
                                ENDIF
                                IF (FLAG .EQV. .FALSE.) THEN
                                        PRINT*, "Matrix is non - invertible"
                                        inverse = 0
                                        errorflag = -1
                                        return
                                ENDIF
                        END DO
                ENDIF
                DO j = k+1, n                        
                        m = augmatrix(j,k)/augmatrix(k,k)
                        DO i = k, 2*n
                                augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
                        END DO
                END DO
        END DO
        
        !write(*,*)'3 sur 6'
        !Test for invertibility
        DO i = 1, n
                IF (augmatrix(i,i) == 0) THEN
                        PRINT*, "Matrix is non - invertible"
                        inverse = 0
                        errorflag = -1
                        return
                ENDIF
        END DO
        
        !write(*,*)'4 sur 6'
        !Make diagonal elements as 1
        DO i = 1 , n
                m = augmatrix(i,i)
                DO j = i , (2 * n)                                
                           augmatrix(i,j) = (augmatrix(i,j) / m)
                END DO
        END DO
        
        !write(*,*)'5 sur 6'
        !Reduced right side half of augmented matrix to identity matrix
        DO k = n-1, 1, -1
                DO i =1, k
                m = augmatrix(i,k+1)
                        DO j = k, (2*n)
                                augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
                        END DO
                END DO
        END DO                                
        
        !write(*,*)'6 sur 6'
        !store answer
        DO i =1, n
                DO j = 1, n
                        inverse(i,j) = augmatrix(i,j+n)
                END DO
        END DO
        errorflag = 0
END SUBROUTINE FINDinv  

!////////////////////////////////
!//// SUBROUTINE D AFFICHAGE ////
!////////////////////////////////

!AFfichage COupe REal
subroutine afcore(vect,nx,ny,n,nom,en)
	real,dimension(:),intent(in)::vect
	integer,intent(in)::nx,ny,n,en
	character(LEN=*),intent(in)::nom
	integer::j
	real::x
	open(UNIT=n,FILE=nom)
	do j=1,ny
		x=j*1.
		write(n,*)x,vect(nx*(j-1)+en)
	end do
		close(n)
end subroutine afcore

!AFfichage NAppe REal
subroutine afnare(vect,nx,ny,n,nom)
	real,dimension(:),intent(in)::vect
	integer,intent(in)::nx,ny,n
	character(LEN=*),intent(in)::nom
	integer::i,j,k
	real::x,y
	open(UNIT=n,FILE=nom)
	do i=1,nx ; do j=1,ny
		k=nx*(i-1)+j ; x=1.*i ; y=1.*j
		write(n,*)x,y,vect(k)
	end do ; end do
	close(n)
end subroutine afnare

!AFfiche Matrice REal... |A|
!subroutine afmare(tab,n,nom)
!	real,dimension(:,:),intent(in)::tab
!	integer,intent(in)::n
!	character(LEN=*),intent(in)::nom
!	integer::l,i
!	open(UNIT=n,FILE=nom)
!	l=size(tab,1)
!	do i=1,l ; write(n,*)tab(i,:) ; end do
!	close(n)
!end subroutine afmare

!AFfiche VEcteur REal... |A|
!subroutine afvere(vect,n,nom)
!	real,dimension(:),intent(in)::vect
!	integer,intent(in)::n
!	character(LEN=*),intent(in)::nom
!	integer::l,i
!	open(UNIT=n,FILE=nom)
!	l=size(vect)
!	do i=1,l ; write(n,*)vect(i) ; end do
!	close(n)
!end subroutine afvere

!affichage d informations diverses au lancement du programme
subroutine apropos(div,nx,c,motif,epaiss,m2,nite,ncapture)
	real,intent(in)::div,c,motif,epaiss,m2
	integer,intent(in)::nx,nite,ncapture
	write(*,*)'--------------------------------'
	write(*,*)''
	write(*,*)'MOuntain TOp REmoval : motore_r1'
	write(*,*)''
	write(*,*)'--------------------------------'
	write(*,*)'>>la topo est un carre de cote 1./div avec div=',div
	write(*,*)'le nb de subdivision d un cote de la topo est de nx=ny=',nx
	write(*,*)''
	write(*,*)'>>constante c=',c
	write(*,*)'l altitude varie de 0 à 12 environ'
	write(*,*)'si l altitude n est pas un multiple de motif=',motif
	write(*,*)'plus ou moins epaiss=',epaiss
	write(*,*)'alors la constante c est multiplie par m2=',m2
	write(*,*)''
	write(*,*)'>>nb iterations=',nite
	write(*,*)'enregistrement des donnes tout les ncapture=',ncapture
	write(*,*)''
	write(*,*)'--------------------------------'
	write(*,*)''
end subroutine apropos

end module m_motore
