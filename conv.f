        subroutine conv4(mat,mx,my,mz,mu,mw,filtermat,ksiz,smoothed) 

        integer mx,my,mz,mu,mw,ksiz,m1,m2,m3,m4,e,i,j,k,m,l
        double precision mat(mx,my,mz,mu,mw)
        double precision filtermat(ksiz,ksiz,ksiz,ksiz)
        double precision smoothed(mx,my,mz,mu,mw),total,f


        e=(ksiz+1)/2
        do i=1,mx
           do j=1,my
              do k=1,mz
                do l=1,mu
                  do m=1,mw
                       total=0.0
                       f=0.0

                       do m1=1,ksiz
                          do m2=1,ksiz
                             do m3=1,ksiz
                               do m4=1,ksiz

      if((0.lt.(i-e+m1)).and.((i-e+m1).le.mx)) then
      if((0.lt.(j-e+m2)).and.((j-e+m2).le.my)) then
      if((0.lt.(k-e+m3)).and.((k-e+m3).le.mz)) then
      if((0.lt.(l-e+m4)).and.((l-e+m4).le.mu)) then

        total=total+filtermat(m1,m2,m3,m4)*
     1  mat(i-e+m1,j-e+m2,k-e+m3,l-e+m4,m)

        f=f+filtermat(m1,m2,m3,m4)
				     end if
                                   end if
                                end if
			      end if
                             end do
                          end do
                       end do
                     end do
                       smoothed(i,j,k,l,m)=total/f

                    end do
                 end do
              end do
           end do
        end do

        return
        end

  	subroutine conv3(mat,mx,my,mz,mw,filtermat,ksiz,smoothed) 

        integer mx,my,mz,mw,ksiz,m1,m2,m3,e,i,j,k,m
        double precision mat(mx,my,mz,mw),filtermat(ksiz,ksiz,ksiz)
        double precision smoothed(mx,my,mz,mw),total,f


        e=(ksiz+1)/2
        do i=1,mx
           do j=1,my
              do k=1,mz
                 do m=1,mw
                       total=0.0
                       f=0.0

                       do m1=1,ksiz
                          do m2=1,ksiz
                             do m3=1,ksiz

       if((0.lt.(i-e+m1)).and.((i-e+m1).le.mx).and.(0.lt.(j-e+m2))) then
       if(((j-e+m2).le.my).and.(0.lt.(k-e+m3)).and.((k-e+m3).le.mz))then


        total=total+filtermat(m1,m2,m3)*mat(i-e+m1,j-e+m2,k-e+m3,m)

        f=f+filtermat(m1,m2,m3)
                                   end if
                                end if
                             end do
                          end do
                       end do
                     smoothed(i,j,k,m)=total/f
                 end do
              end do
           end do
        end do

        return
        end