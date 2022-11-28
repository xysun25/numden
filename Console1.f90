module atomdef
implicit none
type atom
character(len=10)::element
real x,y,z,m,e,sig,eps
end type

type mol
character(len=10)::numberofmol
real x,y,z
end type

end module

program main
use atomdef
implicit none
character(len=20) filename,mmol1,mmol2,mmol3,trjfile
integer num,numberofAIM
integer i,j,k,u,l,o,p,q,r,trj,ls,u1,u2,u3,ff1,ff2,lnn,np,nn,mm,up1,upan,uu1,uu2,uu3,nsite,tt
integer numberoftrj
integer num1,num2,num3,numberofAIM1,numberofAIM2,numberofAIM3
integer numberofatom1,numberofatom2,numberofatom3,numberofatoms_2
type(atom),allocatable::atoms_1(:),atoms_2(:),atoms(:,:)
real X,Y,Z,M,Xcom,Ycom,Zcom,x00,y00,dr,high
real a0,b0,c0,a,b,c,interv,ax,ay,az,bx,by,bz,cx,cy,cz,theta,bin,angle,c1,c2,average
integer numberofinterv
integer::numberofions,numberofatoms
real sumdensityofcation1,sumdensityofanion1,sumdensityofcharge1,sumdensityofcation2,sumdensityofanion2,sumdensityofcharge2
real,allocatable::numberofcation(:,:),numberofanion(:,:),numberofcharge(:,:)
real,allocatable::den_an(:,:),den_ca(:,:),angledis(:,:)



read (*,*) filename
write(*,*) 'the input file is:', filename
open(100,file=filename,status="old")
read(100,*) trjfile,numberoftrj
write(*,*) 'trjfile and number of trj',trjfile,numberoftrj
read (100,*) mmol1,mmol2,mmol3
write(*,*) 'the mmol files are:', mmol1,mmol2,mmol3
read (100,*) num1,num2,num3
write(*,*) 'the number of molecules are:', num1,num2,num3
read(100,*) a0,b0,c0,ax,by,cz
read(100,*) interv, numberofinterv



open(11,file=mmol1,status="old")
read(11,*)
read(11,*) numberofAIM1
!write(*,*) numberofAIM1
allocate(atoms_1(numberofAIM1))
do i=1,numberofAIM1
read(11,*)atoms_1(i)%element,atoms_1(i)%x,atoms_1(i)%y,atoms_1(i)%z,atoms_1(i)%m,atoms_1(i)%e,atoms_1(i)%sig,atoms_1(i)%eps
!write (*,*) atoms_1(i)%m
end do
rewind(11)

open(12,file=mmol2,status="old")
read(12,*)
read(12,*) numberofAIM2
!write(*,*) numberofAIM2
allocate (atoms_2(numberofAIM2))
do j=1,numberofAIM2
read(12,*)atoms_2(j)%element,atoms_2(j)%x,atoms_2(j)%y,atoms_2(j)%z,atoms_2(j)%m,atoms_2(j)%e,atoms_2(j)%sig,atoms_2(j)%eps
!write (*,*) atoms_2(j)%m
end do
rewind(12)

open(13,file=mmol3,status="old")
read(13,*) numberofAIM3
close(13)
numberofatom3=num3*numberofAIM3

num=num1*numberofAIM1+num2*numberofAIM2
numberofAIM=numberofAIM1+numberofAIM2


open(10,file=trjfile,status="old")
allocate(numberofcation(2,numberofinterv))
allocate(numberofanion(2,numberofinterv))
allocate(numberofcharge(2,numberofinterv))
!allocate(numberofcation2(numberoftrj,numberofinterv))
!allocate(numberofanion2(numberoftrj,numberofinterv))
!allocate(numberofcharge2(numberoftrj,numberofinterv))


ff1 = 0 
ff2 = 0
upan = 0
up1 = 0
np = 100


do trj=1,numberoftrj
allocate(atoms(num,numberofAIM))
read(10,*)
read(10,*)
!do k=1,num
!do u=1,numberofAIM
!read(10,*) atoms(k,u)%element,atoms(k,u)%x,atoms(k,u)%y,atoms(k,u)%z
!end do
!end do
!numberofdensity
    do k=1,num1
        do u=1,numberofAIM1
        read(10,*) atoms(k,u)%element,atoms(k,u)%x,atoms(k,u)%y,atoms(k,u)%z
        end do
    end do

    do k=num1+1,num2+num1
        do o=1,numberofAIM2
        read (10,*)atoms(k,o)%element,atoms(k,o)%x,atoms(k,o)%y,atoms(k,o)%z
        end do
    end do

!tune the Cation across the boundary
! transfer CA
   do k=1,num1
    tt = 0
!!! jude the x boundary
	x00 = atoms(k,1)%x
        do u=2,numberofAIM1
          dr = abs(atoms(k,u)%x-x00)
          if(dr >(0-a0)) then
           tt =1 !cross the bounder
          end if
		end do
        if (tt ==1) then
          do u=1,numberofAIM1
            if(atoms(k,u)%x > 0.0) then
              atoms(k,u)%x = atoms(k,u)%x + a0*2
            end if
          end do
        end if
!!! jude the y boundary
    y00 = atoms(k,1)%y
    tt = 0
        do u=2,numberofAIM1
          dr = abs(atoms(k,u)%y-y00)
          if(dr >(0-b0)) then
           tt =1 !cross the bounder
          end if
        end do
        if (tt == 1) then
          do u=1,numberofAIM1
            if(atoms(k,u)%y >0.0) then
              atoms(k,u)%y = atoms(k,u)%y + b0*2
            end if
          end do
        end if
  end do 


! transfer AN
   do k=num1+1,num2+num1
    tt = 0
!!! jude the x boundary
    x00 = atoms(k,1)%x
        do u=2,numberofAIM2
          dr = abs(atoms(k,u)%x-x00)
          if(dr >(0-a0)) then
           tt =1 !cross the bounder
          end if
        end do
        if (tt ==1) then
          do u=1,numberofAIM2
            if(atoms(k,u)%x > 0.0) then
              atoms(k,u)%x = atoms(k,u)%x + a0*2
            end if
          end do
        end if
!!! jude the y boundary
    y00 = atoms(k,1)%y
    tt = 0
        do u=2,numberofAIM2
          dr = abs(atoms(k,u)%y-y00)
          if(dr >(0-b0)) then
           tt =1 !cross the bounder
          end if
        end do
        if (tt == 1) then
          do u=1,numberofAIM2
            if(atoms(k,u)%y >0.0) then
              atoms(k,u)%y = atoms(k,u)%y + b0*2
            end if
          end do
        end if
   end do


    do r=1,numberofatom3
    read(10,*) 
	end do
	
do k=1,num1
         do u=1,numberofAIM1   ! find the cation com 
			 X=X+atoms(k,u)%x*atoms_1(u)%m
			 Y=Y+atoms(k,u)%y*atoms_1(u)%m
			 Z=Z+atoms(k,u)%z*atoms_1(u)%m
			 M=M+atoms_1(u)%m
			 do ls=1,numberofinterv
				 b=by+(ls-1)*interv
				 if(atoms(k,u)%y>=b.and.atoms(k,u)%y<(b+interv))then
					 !write(*,*) "com-trans completed"
					 if(atoms(k,u)%z<cz) then
						 numberofcharge(1,ls)=numberofcharge(1,ls)+atoms_1(u)%e/numberoftrj
					 elseif(atoms(k,u)%z>(0-cz))then
						 numberofcharge(2,ls)=numberofcharge(2,ls)+atoms_1(u)%e/numberoftrj
					 end if
				 end if
			 end do
		 end do 
     Xcom=X/M
	 Ycom=Y/M
	 Zcom=Z/M
     do ls=1,numberofinterv
		 b=by+(ls-1)*interv
        if(Ycom>=b.and.Ycom<(b+interv))then
			if(Zcom<cz)then
				numberofcation(1,ls)=numberofcation(1,ls)+1.0/numberoftrj
			elseif(Zcom>(0-cz))then
				numberofcation(2,ls)=numberofcation(2,ls)+1.0/numberoftrj
			end if
		end if
    end do
    X=0.0
    Y=0.0
    Z=0.0
    M=0.0
end do 

do k=num1+1,num2+num1
      do o=1,numberofAIM2
		  X=X+atoms(k,o)%x*atoms_2(o)%m
		  Y=Y+atoms(k,o)%y*atoms_2(o)%m
		  Z=Z+atoms(k,o)%z*atoms_2(o)%m
		  M=M+atoms_2(o)%m
		  do ls=1,numberofinterv
				 b=by+(ls-1)*interv
				 if(atoms(k,o)%y>=b.and.atoms(k,o)%y<(b+interv))then
					 if(atoms(k,o)%z<cz) then
						 numberofcharge(1,ls)=numberofcharge(1,ls)+atoms_2(o)%e/numberoftrj
					 elseif(atoms(k,o)%z>(0-cz))then
						 numberofcharge(2,ls)=numberofcharge(2,ls)+atoms_2(o)%e/numberoftrj
					 end if
				 end if
			 end do
	  end do 
	  Xcom=X/M
	  Ycom=Y/M
	  Zcom=Z/M
	  do ls=1,numberofinterv
		 b=by+(ls-1)*interv
		 if(Ycom>=b.and.Ycom<(b+interv))then
			if(Zcom<cz)then
				numberofanion(1,ls)=numberofanion(1,ls)+1.0/numberoftrj
			elseif(Zcom>(0-cz))then
				numberofanion(2,ls)=numberofanion(2,ls)+1.0/numberoftrj
			end if
		 end if
	 end do
    X=0.0
    Y=0.0
    Z=0.0
    M=0.0
end do 
deallocate(atoms)
end do


!average density
open(14,file="average",status="unknown")
do i=1,numberofinterv
	average=1000/interv/(ax-a0)/(cz-c0)
        sumdensityofcation1=numberofcation(1,i)*average
        sumdensityofanion1=numberofanion(1,i)*average
        sumdensityofcharge1=numberofcharge(1,i)*average
		sumdensityofcation2=numberofcation(2,i)*average
        sumdensityofanion2=numberofanion(2,i)*average
        sumdensityofcharge2=numberofcharge(2,i)*average
    write(14,"(I4,2x,6F10.4)")ls,sumdensityofcation1,sumdensityofanion1,sumdensityofcharge1,&
sumdensityofcation2,sumdensityofanion2,sumdensityofcharge2
	!sumdensityofcation1=0
 !   sumdensityofanion1=0
 !   sumdensityofcharge1=0
	!sumdensityofcation2=0
 !   sumdensityofanion2=0
 !   sumdensityofcharge2=0
end do
deallocate(numberofcation)
deallocate(numberofanion)
deallocate(numberofcharge)
!deallocate(numberofcation2)
!deallocate(numberofanion2)
!deallocate(numberofcharge2)

stop
end program
