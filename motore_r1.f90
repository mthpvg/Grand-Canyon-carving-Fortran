program motore
use m_motore
implicit none

	integer,parameter::nt=100,nx=26,ny=26
	real,parameter::div=2.
	real,dimension(nx*ny)::v_relief,v_cst,v_dcx,v_dcy,v_relief_p1
	real,dimension(nx*ny,nx*ny)::mat,invmat
	real::dx,dy,dt,c,m1,m2,motif,epaiss
	integer::nxy,k,error,l,compteur,i,nite,ncapture
	character(LEN=16)::nom,nom2,nom3
	
!VARIABLES
	!v_relief contient la hauteur de tous les points de la topo, la subroutine m_topomulti lui donne ses valeurs
	!v_cst contient la valeur du coefficient c de tous les points de la topo, il depend de la hauteur de ces derniers la subroutine c_relief lui donne ses valeurs
	!v_dcx est la dérivée de v_cst selon x
	!v_dcy est la dérivée de v_cst selon y
	!Le programme effectue des iterations sur l equation : v_relief_p1=invmat*v_relief

!Initialistion des variables
	nxy=nx*ny
	dx=1./div/nx ; dy=1./div/ny ; dt=1./nt
	v_relief=0. ; v_cst=0. ; v_dcx=0. ; v_dcy=0. ; mat=0. ; v_relief_p1=0.
	error=0 ; compteur=0
	c=0.0001 !4
	m1=3. ; m2=3. ; motif=2. ; epaiss=0.2!0.1
	nite=600 ; ncapture=3
	
	call apropos(div,nx,c,motif,epaiss,m2,nite,ncapture)

!Construction et ecriture de v_relief
	call m_topomulti(v_relief,nx,ny,div)
	call afnare(v_relief,nx,ny,1,'v_relief0.res')
	call afcore(v_relief,nx,ny,9,'v_1coupe0.res',9)
	call afcore(v_relief,nx,ny,9,'v_2coupe0.res',12)
	!call afvere(v_relief,14,'v_relief.res')
	write(*,*)'relief construit'

!Boucle iterative
	do l=1,nite
		compteur=compteur+1
	!Construction de v_cst
		call c_relief(v_relief,v_cst,c,m1,m2,motif,epaiss)
		!call afvere(v_cst,11,'v_cst.res')
		do i=1,nx
		!v_cst(nx*(9-1)+i)=m2*c*1.1
		end do
	!Construction de v_dcx
		call constru_dcx(v_dcx,v_cst,nx,ny,nxy)
		!v_dcx=0.
		!call afvere(v_dcx,12,'v_dcx.res')
	!Construction de v_dcy
		call constru_dcy(v_dcy,v_cst,nx,ny,nxy)
		!v_dcy=0.
		!call afvere(v_dcy,13,'v_dcy.res')
	!Construction de la matrice
		call m_matrice(mat,0,1,nx,dt,dy,v_cst,nx,-1.,v_dcy)
		call m_matrice(mat,0,nx,0,dt,dx,v_cst,nx,-1.,v_dcx)
		call m_matrice(mat,1,0,nx,dt,dy,v_cst,nx,1.,v_dcy)
		call m_matrice(mat,nx,0,0,dt,dx,v_cst,nx,1.,v_dcx)
		do k=1,nxy
			mat(k,k)=1.+2.*v_cst(k)*dt*(1./dx**2.+1./dy**2.)
		end do
		!call afmare(mat,10,'matrice.dat')
	!Inversion de la matrice
		call FINDInv(mat,invmat,nxy,error)
	!Calcul et ecriture du nouveau vecteur v_relief
		v_relief_p1=matmul(invmat,v_relief)
	!write(*,*)sum(abs(v_relief_p1-v_relief))
		v_relief=v_relief_p1
		if (compteur==ncapture) then
			write (nom(9:12),'(I4)')l+1000
			nom(1:8)='v_relief'
			nom(13:16)='.res'
			call afnare(v_relief,nx,ny,7,nom)

			write (nom2(9:12),'(I4)')l+1000
			nom2(1:8)='v_1coupe'
			nom2(13:16)='.res'
			call afcore(v_relief,nx,ny,8,nom2,9)

			write (nom3(9:12),'(I4)')l+1000
			nom3(1:8)='v_2coupe'
			nom3(13:16)='.res'
			call afcore(v_relief,nx,ny,9,nom3,12)
			compteur=0
		end if

		write(*,*)'iteration',l
	end do
end program motore
