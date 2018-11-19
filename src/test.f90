! subroutine test(tt0,tt1,delta,trialref,patientref,trt,nrow,som)
subroutine test(donnee,nrow,som)

implicit none

integer, intent(in):: nrow
double precision, dimension(nrow,6), intent(inout)::donnee
! double precision, dimension(3288), intent(in)::tt0,tt1,delta,trialref,patientref,trt
double precision, dimension(nrow), intent(out)::som
integer::i

do i=1,nrow
    ! som(i)=tt0(i)+tt1(i)+delta(i)+trialref(i)+patientref(i)+trt(i)
    som(i)=sum(donnee(i,:))
enddo

endsubroutine test